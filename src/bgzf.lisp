;;;
;;; Copyright (C) 2009-2010 Keith James. All rights reserved.
;;;
;;; This file is part of cl-sam.
;;;
;;; This program is free software: you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation, either version 3 of the License, or
;;; (at your option) any later version.
;;;
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with this program.  If not, see <http://www.gnu.org/licenses/>.
;;;

(in-package :sam)

(defconstant +xlen+ 6)
(defconstant +sf1+ (char-code #\B))
(defconstant +sf2+ (char-code #\C))
(defconstant +slen+ 2)

(defconstant +member-header-length+ 18)
(defconstant +member-footer-length+ 8)
(defconstant +bgz-max-size+ (expt 2 16))

(defvar *empty-bgzf-record*
  (make-array 28 :element-type 'octet
              :initial-contents
              '(#o037 #o213 #o010 #o004 #o000 #o000 #o000 #o000 #o000 #o377
                #o006 #o000 #o102 #o103 #o002 #o000 #o033 #o000 #o003 #o000
                #o000 #o000 #o000 #o000 #o000 #o000 #o000 #o000))
  "SAMtools version >0.1.5 appends an empty BGZF record to allow
  detection of truncated BAM files. These 28 bytes constitute such a
  record.")

(defvar *bgz-read-buffer* (make-array 8 :element-type 'octet :initial-element 0)
  "The buffer used by {defun read-bgz-member} for reading block gzip
data. Rebind this per thread to make {defun read-bgz-member}
re-entrant.")

(defvar *bgz-write-buffer* (make-array 8 :element-type 'octet
                                       :initial-element 0)
  "The buffer used by {defun write-bgz-member} for writing block gzip
data. Rebind this per thread to make {defun write-bgz-member}
re-entrant.")

(defstruct (bgz-member (:include gz:gz-member))
  "A gzip member with extensions, as defined by RFC1952. The
extensions are described in RFC1952 and the SAM format specification."
  (sf1 +sf1+ :type uint8)
  (sf2 +sf2+ :type uint8)
  (slen +slen+ :type uint16)
  (bsize 0 :type uint16)
  (udata (make-array 0 :element-type 'octet) :type simple-octet-vector))

(defun read-bgz-member (stream &optional (buffer *bgz-read-buffer*))
  (declare (optimize (speed 3)))
  (flet ((decode-bytes (s n)
           (when (plusp (read-sequence buffer s :end n))
             (ecase n
               (1 (decode-uint8le buffer))
               (2 (decode-uint16le buffer))
               (4 (decode-uint32le buffer))))))
    (let ((id1 (decode-bytes stream 1)))
      (when id1
        (let* ((id2 (decode-bytes stream 1))
               (cm (decode-bytes stream 1))
               (flg (decode-bytes stream 1))
               (mtime (decode-bytes stream 4))
               (xfl (decode-bytes stream 1))
               (os (decode-bytes stream 1))
               (xlen (decode-bytes stream 2))
               (sf1 (decode-bytes stream 1))
               (sf2 (decode-bytes stream 1))
               (slen (decode-bytes stream 2))
               (bsize (decode-bytes stream 2))) ; 18 header bytes
          (declare (ignore xfl))
          (assert (and (= gz:+id1+ id1)
                       (= gz:+id2+ id2)
                       (= gz:+cm-deflate+ cm)
                       (plusp (logand flg gz:+flag-extra+))
                       (= +xlen+ xlen)
                       (= +sf1+ sf1)
                       (= +sf2+ sf2)
                       (= +slen+ slen)) () "Invalid block gzip header")
          (let* ((deflated-size (1+ (logand bsize #xffff)))
                 (cdata (make-array (- deflated-size
                                       +member-header-length+
                                       +member-footer-length+)
                                    :element-type '(unsigned-byte 8))))
            (read-sequence cdata stream)
            (let* ((crc32 (decode-bytes stream 4))
                   (isize (decode-bytes stream 4))) ; 8 footer bytes
              (make-bgz-member :mtime mtime :os os :xlen xlen
                               :bsize deflated-size :cdata cdata
                               :crc32 crc32 :isize isize))))))))

(defun write-bgz-member (bgz stream &optional (buffer *bgz-write-buffer*))
  (declare (optimize (speed 3)))
  (flet ((encode-bytes (x n s)
           (ecase n
             (1 (encode-int8le x buffer))
             (2 (encode-int16le x buffer))
             (4 (encode-int32le x buffer)))
           (write-sequence buffer s :end n)))
    (encode-bytes (bgz-member-id1 bgz) 1 stream)
    (encode-bytes (bgz-member-id2 bgz) 1 stream)
    (encode-bytes (bgz-member-cm bgz) 1 stream)
    (encode-bytes (bgz-member-flg bgz) 1 stream)
    (encode-bytes (get-universal-time) 4 stream)
    (encode-bytes (bgz-member-xfl bgz) 1 stream)
    (encode-bytes gz:+os-unknown+ 1 stream)
    (encode-bytes (bgz-member-xlen bgz) 2 stream)
    (encode-bytes (bgz-member-sf1 bgz) 1 stream)
    (encode-bytes (bgz-member-sf2 bgz) 1 stream)
    (encode-bytes (bgz-member-slen bgz) 2 stream)
    (encode-bytes (1- (bgz-member-bsize bgz)) 2 stream) ; 18 header bytes
    (write-sequence (bgz-member-cdata bgz) stream) ; Payload
    (encode-bytes (bgz-member-crc32 bgz) 4 stream)
    (encode-bytes (bgz-member-isize bgz) 4 stream) ; 8 footer bytes
    bgz))

(defstruct bgzf
  (pathname nil :type t)
  (stream nil :type t)
  (buffer (make-array +bgz-max-size+ :element-type 'octet)
          :type simple-octet-vector)
  (position 0 :type (unsigned-byte 48))
  (offset 0 :type uint16)
  (pointer 0 :type uint16)
  (init nil :type t)
  (eof nil :type t))

(defmacro with-bgzf-file ((var filespec &key (direction :input)
                               (if-exists :supersede))
                          &body body)
  "Executes BODY with VAR bound to a BGZF handle structure created by
opening the file denoted by FILESPEC.

Arguments:

- var (symbol): The symbol to be bound.
- filespec (pathname designator): The file to open.

Key:

- direction (keyword): The direction, one of either :read or :write ."
  `(let ((,var (bgzf-open ,filespec :direction ,direction
                          :if-exists ,if-exists)))
    (unwind-protect
         (progn
           ,@body)
      (when ,var
        (bgzf-close ,var)))))

(defun bgzf-open (filespec &key (direction :input) (if-exists :supersede))
  "Opens a block gzip file for reading or writing.

Arguments:

- filespec (pathname designator): The file to open.

Key:

- direction (keyword): The direction, one of either :read or :write .

Returns:

- A BGZF structure."
  (let ((stream (open filespec :element-type 'octet :direction direction
                      :if-exists if-exists)))
    (make-bgzf :stream stream :pathname filespec)))

(defun bgzf-close (bgzf)
  "Closes an open block gzip handle.

Arguments:

- bgzf (bgzf structure): The handle to close.

Returns:

- T on success."
  (let ((stream (bgzf-stream bgzf)))
    (when (output-stream-p stream)
       (bgzf-flush bgzf))
    (close stream)))

(defun bgzf-open-p (bgzf)
  (open-stream-p (bgzf-stream bgzf)))

(defun bgzf-seek (bgzf position)
  "Seeks with the file encapsulated by a block gzip handle.

Arguments:

- bgzf (bgzf structure): The handle to seek.
- position (integer): The position to seek. Only values previously
  returned by {defun bgzf-tell} may be used.

Returns:

- The new position."
  (let ((stream (bgzf-stream bgzf))
        (new-position (ash position -16))
        (new-offset (logand position #xffff))
        (current-position (bgzf-position bgzf)))
    (cond ((= current-position new-position) ; Seek within current bgz member
           (setf (bgzf-offset bgzf) new-offset))
          ((file-position stream new-position)
           (setf (bgzf-buffer bgzf) (bgz-member-udata (read-bgz-member stream))
                 (bgzf-offset bgzf) new-offset
                 (bgzf-position bgzf) new-position))
          (t
           (error 'bgzf-io-error :text "failed to seek")))))

(defun bgzf-tell (bgzf)
  "Returns the current position in the encapsulated file of a block gzip
handle.

Arguments:

- bgzf (bgzf structure): The handle.

Returns:

- The file position."
  (logior (ash (bgzf-position bgzf) 16) (bgzf-offset bgzf)))

;; As yet untested
(defun bgzf-eof-p (bgzf)
  (let* ((stream (bgzf-stream bgzf))
         (pos (file-position stream)))
    (unwind-protect
         (cond ((< (file-length stream) (length *empty-bgzf-record*))
                (error 'bgzf-io-error :text "incomplete file"))
               ((file-position stream (- (file-length stream)
                                         (length *empty-bgzf-record*)))
                (loop
                   for byte across *empty-bgzf-record*
                   always (equal byte (read-byte stream))))
               (t
                nil))
      (file-position stream pos))))
