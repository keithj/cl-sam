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

(deftype bgz-payload-index ()
  '(integer 0 65536))

(defconstant +xlen+ 6
  "The value of the RCF1952 eXtra LENgth field used by BGZ members.")
(defconstant +sf1+ (char-code #\B)
  "The value of the RCF1952 first extra subfield used by BGZ members.")
(defconstant +sf2+ (char-code #\C)
  "The value of the RCF1952 first extra subfield used by BGZ members.")
(defconstant +slen+ 2
  "The value of the RFC1952 subfield length field used by BGZ members.")

(defconstant +member-header-length+ 18
  "The total number of bytes in the BGZ header.")
(defconstant +member-footer-length+ 8
  "The total number of bytes in the BGZ footer.")

(defconstant +bgz-max-payload-length+ (expt 2 16)
  "The maximium size in bytes of a BGZ member payload. This is
  dictated by the fact that the SAM spec makes 16 bits are avaliable
  for addressing positions with the BGZ member.")

(defvar *empty-bgz-record*
  (make-array 28 :element-type 'octet
              :initial-contents
              '(#o037 #o213 #o010 #o004 #o000 #o000 #o000 #o000 #o000 #o377
                #o006 #o000 #o102 #o103 #o002 #o000 #o033 #o000 #o003 #o000
                #o000 #o000 #o000 #o000 #o000 #o000 #o000 #o000))
  "SAMtools version >0.1.5 appends an empty BGZF record to allow
  detection of truncated BAM files. These 28 bytes constitute such a
  record.")

(defstruct (bgz-member (:include gz:gz-member))
  "A block-gzip data chunk; a gzip member with extensions, as defined
by RFC1952. The extensions are described in RFC1952 and the SAM format
specification.

- sf1: RCF1952 first extra SubField.
- sf2: RCF1952 second extra SubField.
- slen: RFC1952 Subfield LENgth.
- bsize: SAM spec total Block (member) SIZE. The serialized value is
  bsize -1.
- udata: Uncompressed DATA."
  (sf1 +sf1+ :type uint8 :read-only t)
  (sf2 +sf2+ :type uint8 :read-only t)
  (slen +slen+ :type uint16 :read-only t)
  (bsize 0 :type uint16)
  (udata (make-array 0 :element-type 'octet) :type simple-octet-vector))

(defstruct bgzf
  "A block gzip file.

- pathname: The file pathname.
- stream: The file stream.
- buffer: A simple-octet-vector used for buffering reads.
- position: The file position component of the BGZF virtual file
  offset (most significant 48 bits). This is the position in the
  stream at which this member starts.
- offset: The within-member offset component of the BGZF virtual file
  offset (least significant 16 bits). This is a position within the
  uncompressed data of the member.

- init: T if the struct has been initialized (used internally in
  decompression and reading).

- eof: T if the decompression process has reached EOF (used internally
  in decompression and reading)."
  (pathname nil :type t)
  (stream nil :type t)
  (compression 5 :type (integer 0 9))
  (buffer (make-array +bgz-max-payload-length+ :element-type 'octet)
          :type simple-octet-vector)
  (position 0 :type (unsigned-byte 48))
  (offset 0 :type uint16)
  (pointer 0 :type uint16)
  (init nil :type t)
  (eof nil :type t))

(defmacro with-bgzf ((var filespec &key (direction :input)
                          (if-exists :overwrite)) &body body)
  "Executes BODY with VAR bound to a BGZF handle structure created by
opening the file denoted by FILESPEC.

Arguments:

- var (symbol): The symbol to be bound.
- filespec (pathname designator): The file to open.

Key:

- direction (keyword): The direction, one of either :input or :output .
- if-exists (keyword): Behaviour with respoect to existing
  files. Defaults to :overwrite ."
  `(let ((,var (bgzf-open ,filespec :direction ,direction
                          :if-exists ,if-exists)))
    (unwind-protect
         (progn
           ,@body)
      (when ,var
        (bgzf-close ,var)))))

(defun bgzf-open (filespec &key (direction :input) compression
                  (if-exists :overwrite))
  "Opens a block gzip file for reading or writing.

Arguments:

- filespec (pathname designator): The file to open.

Key:

- direction (keyword): The direction, one of either :read or :write .
- if-exists (keyword): Behaviour with respoect to existing
  files. Defaults to :overwrite .

Returns:

- A BGZF structure."
  (let ((stream (open filespec :element-type 'octet :direction direction
                      :if-exists if-exists)))
    (if compression
        (make-bgzf :stream stream :pathname filespec :compression compression)
      (make-bgzf :stream stream :pathname filespec))))

(defun bgzf-close (bgzf)
  "Closes an open block gzip file.

Arguments:

- bgzf (bgzf structure): The file to close.

Returns:

- T on success."
  (let ((stream (bgzf-stream bgzf)))
    (when (output-stream-p stream)
       (bgzf-flush bgzf))
    (close stream)))

(defun bgzf-open-p (bgzf)
  (open-stream-p (bgzf-stream bgzf)))

(defun bgzf-seek (bgzf position)
  "Seeks with the file encapsulated by a block gzip file.

Arguments:

- bgzf (bgzf structure): The handle to seek.

- position (integer): The position to seek. Only values previously
  returned by {defun bgzf-tell} may be used. The most significant 48
  bits denote the real file position and the least significant 16 bits
  the offset within the uncompressed gzip member (see the SAM spec).

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
file.

Arguments:

- bgzf (bgzf structure): The file.

Returns:

- The file position."
  (logior (ash (bgzf-position bgzf) 16) (bgzf-offset bgzf)))

;; As yet untested
(defun bgzf-eof-p (bgzf)
  (let* ((stream (bgzf-stream bgzf))
         (pos (file-position stream)))
    (unwind-protect
         (cond ((< (file-length stream) (length *empty-bgz-record*))
                (error 'bgzf-io-error :text "incomplete file"))
               ((file-position stream (- (file-length stream)
                                         (length *empty-bgz-record*)))
                (loop
                   for byte across *empty-bgz-record*
                   always (equal byte (read-byte stream))))
               (t
                nil))
      (file-position stream pos))))
