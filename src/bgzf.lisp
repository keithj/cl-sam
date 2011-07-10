;;;
;;; Copyright (c) 2009-2011 Keith James. All rights reserved.
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

(defvar *default-compression* 5
  "The default bgz compression level (Zlib level 5).")

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
- bsize: SAM spec total Block (member) SIZE. The serialized value is (1- bsize).
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
- pointer: The buffer index after the last usable data element. i.e.
  the index at which new data may be written to the buffer prior to
  deflating. Calling (subseq buffer 0 bgzf-pointer) will extract the
  current usable bytes.
- util-buffer: a 4 byte buffer used internally to store an integer.
- loaded-p: T if the bgz data has been loaded into the buffer (used
  internally in decompression and reading).
- load-seek: On loading the bgz data, proceed immediately to this
  offset (used internally in decompression and reading).
- eof: T if the decompression process has reached EOF (used internally
  in decompression and reading)."
  (pathname nil :type t)
  (stream nil :type t)
  (compression *default-compression* :type (integer 0 9))
  (buffer (make-array +bgz-max-payload-length+ :element-type 'octet)
          :type simple-octet-vector)
  (position 0 :type (unsigned-byte 48))
  (offset 0 :type uint16)
  (pointer 0 :type (integer 0 65536))
  (util-buffer (make-array 4 :element-type 'octet :initial-element 0)
               :type (simple-array octet (4)))
  (loaded-p nil :type t)
  (load-seek 0 :type uint16)
  (eof nil :type t))

(defmacro with-bgzf ((var filespec &rest args) &body body)
  "Executes BODY with VAR bound to a BGZF handle structure created by
opening the file denoted by FILESPEC.

Arguments:

- var (symbol): The symbol to be bound.
- filespec (pathname designator): The file to open.

Rest:

- args: Arguments applicable to bgzf-open."
  `(let ((,var (bgzf-open ,filespec ,@args)))
    (unwind-protect
         (progn
           ,@body)
      (when ,var
        (bgzf-close ,var)))))

(defun bgzf-open (bgzfspec &rest args &key (compression *default-compression*)
                  &allow-other-keys)
  "Opens a block gzip file for reading or writing.

Arguments:

- bgzfspec (pathname designator, open stream or bgzf object): The file
  to open or an open stream or bgzf object (which is returned unmodified).
  The *standard-input* or *standard-ouput* streams may be used.

Key:

- compression (keyword): The zlib compression level (for writing).

Also accepts the keyword arguments applicable to CL:OPEN. However,
an :element-type argument will be ignored as the type is always octet.

Returns:

- A BGZF structure."
  (typecase bgzfspec
    (bgzf bgzfspec)
    (t (let ((stream (if (streamp bgzfspec)
                         bgzfspec
                         (apply #'open bgzfspec :element-type 'octet
                                (remove-key-values '(:element-type
                                                     :compression) args)))))
         (make-bgzf :stream stream :pathname bgzfspec
                    :compression compression)))))

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

(defun bgzf-coffset (position)
  "Returns the compressed offset (most significant bits) component of
BGZF virtual offset POSITION."
  (logand (ash position -16) #xffffffffffff))

(defun bgzf-uoffset (position)
  "Returns the uncompressed (least significant bits) component of BGZF
virtual offset POSITION."
  (logand position #xffff))

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
        (new-position (bgzf-coffset position))
        (new-offset (bgzf-uoffset position))
        (current-position (bgzf-position bgzf)))
    (cond ((and (bgzf-loaded-p bgzf)
                (= current-position new-position))
           ;; Pseudo-seek within the current initialised bgz member
           (setf (bgzf-offset bgzf) new-offset))
          ((= current-position new-position)
           ;; Pseudo-seek within the current uninitialised bgz member
           (setf (bgzf-load-seek bgzf) new-offset))
          ((file-position stream  new-position)
           ;; Real seek ensuring that the new bgz member will be
           ;; initialised
           (setf (bgzf-position bgzf) new-position
                 (bgzf-loaded-p bgzf) nil
                 (bgzf-load-seek bgzf) new-offset))
          (t
           (error 'bgzf-io-error :bgzf bgzf :text "failed to seek")))))

(defun bgzf-tell (bgzf)
  "Returns the current position in the encapsulated file of a block gzip
file.

Arguments:

- bgzf (bgzf structure): The file.

Returns:

- The file position."
  (logior (ash (bgzf-position bgzf) 16) (bgzf-offset bgzf)))

(defun bgzf-empty-p (bgzf)
  "Returns T if all decompressed bytes have been read from the current
decompressed BGZ block."
  (= (bgzf-offset bgzf) (1- (length (bgzf-buffer bgzf)))))

;; As yet untested
(defun bgzf-eof-p (bgzf)
  "Returns T if the BGZF stream is terminated by an empty record."
  (let* ((stream (bgzf-stream bgzf))
         (pos (file-position stream)))
    (unwind-protect
         (cond ((< (file-length stream) (length *empty-bgz-record*))
                (error 'bgzf-io-error :bgzf bgzf :text "incomplete file"))
               ((file-position stream (- (file-length stream)
                                         (length *empty-bgz-record*)))
                (loop
                   for byte across *empty-bgz-record*
                   always (equal byte (read-byte stream))))
               (t
                nil))
      (file-position stream pos))))
