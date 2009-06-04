;;;
;;; Copyright (C) 2009 Keith James. All rights reserved.
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

(defconstant +bam-magic+ #(66 65 77 1)
  "The BAM file magic header bytes.")
(defconstant +null-byte+ #x00
  "The termination byte for BAM strings.")
(defconstant +tag-size+ 2
  "The size of a BAM auxilliary tag in bytes.")

(defstruct bgzf
  (file "" :type t)
  (ptr nil :type t)
  (open-p nil :type t))

(defun read-bytes (bgzf n)
  "Reads N bytes from the handle BGZF and returns them as a Lisp array
of unsigned-byte 8. If fewer than N bytes are available an a
{define-condition bgzf-read-error} is raised."
  (declare (optimize (speed 3)))
  (declare (type fixnum n))
  (with-foreign-object (array-ptr :unsigned-char n)
    (let ((num-read (bgzf-ffi:bgzf-read (bgzf-ptr bgzf) array-ptr n)))
      (declare (type fixnum num-read))
      (cond ((zerop num-read)
             nil)
            ((< num-read n)
             (error 'bgzf-io-error
                    :text (format nil #.(txt "expected to read ~a bytes"
                                             "but only ~a were available")
                                  n num-read)))
            (t
             (let ((buffer (make-array n :element-type '(unsigned-byte 8))))
               (loop
                  for i from 0 below n
                  do (setf (aref buffer i)
                           (mem-aref array-ptr :unsigned-char i)))
               buffer))))))

(defun read-string (bgzf n &key null-terminated)
  "Reads N characters from handle BGZF and returns them as a Lisp
string. The NULL-TERMINATED keyword is used to indicate whether the C
string is null-terminated so that the terminator may be consumed."
  (let ((len (if null-terminated
                 (1- n)
               n)))
    (let ((bytes (read-bytes bgzf n)))
      (make-sb-string bytes 0 (1- len)))))

(defun bgzf-open (filespec &key (direction :input))
  (let ((ptr (bgzf-ffi:bgzf-open (namestring filespec) (ecase direction
                                                         (:input "r")
                                                         (:output "w")))))
    (when (null-pointer-p ptr)
      (error 'bgzf-io-error :errno unix-ffi:*error-number*
             :text (format nil "failed to open ~a" filespec)))
    (make-bgzf :file filespec :ptr ptr :open-p t)))

(defun bgzf-close (bgzf)
  (when (bgzf-open-p bgzf)
    (if (zerop (bgzf-ffi:bgzf-close (bgzf-ptr bgzf)))
        t
      (error 'bgzf-io-error :errno unix-ffi:*error-number*
             :text (format nil "failed to close ~a cleanly" bgzf)))))

(defun read-bam-magic (bgzf)
  "Reads the BAM magic number from the handle BGZF and returns T if it
is valid or raises a {define-condition malformed-file-error} if not."
  (let ((magic (read-bytes bgzf 4)))
    (unless (equalp +bam-magic+ magic)
      (error 'malformed-file-error :text "invalid BAM magic number"))
    t))

(defun read-bam-header (bgzf)
  "Returns the unparsed BAM header from the handle BGZF as a Lisp
string."
  (let ((len (decode-int32le (read-bytes bgzf 4))))
    (when (plusp len)
      (read-string bgzf len))))

(defun read-num-references (bgzf)
  "Returns the number of reference sequences described in the BAM
header of handle BGZF."
  (decode-int32le (read-bytes bgzf 4)))

(defun read-reference-meta (bgzf)
  "Returns the reference sequence metadata for a single reference
sequence described in the BAM header of handle BGZF. Two values are
returned, the reference name as a string and the reference length as
an integer,"
  (let ((len (decode-int32le (read-bytes bgzf 4))))
    (values (read-string bgzf len :null-terminated t)
            (decode-int32le (read-bytes bgzf 4)))))

(defun read-alignment (bgzf)
  "Reads one alignment block from handle BGZF, returns it as a Lisp
array of unsigned-byte 8. The handle is advanced to the next
alignment. If no more alignments are available, returns NIL."
  (let ((size-bytes (read-bytes bgzf 4)))
    (if (null size-bytes)
        nil
      (let ((block-size (decode-int32le size-bytes)))
        (read-bytes bgzf block-size)))))

(defun reference-id (alignment-record)
  (decode-int32le alignment-record 0))

(defun alignment-position (alignment-record)
  (decode-int32le alignment-record 4))

(defun read-name-length (alignment-record)
  (decode-uint8le alignment-record 8))

(defun mapping-quality (alignment-record)
  (decode-uint8le alignment-record 9))

(defun alignment-bin (alignment-record)
  (decode-uint16le alignment-record 10))

(defun cigar-length (alignment-record)
  (decode-uint16le alignment-record 12))

(defun alignment-flag (alignment-record)
  (decode-uint16le alignment-record 14))

(defun read-length (alignment-record)
  (decode-int32le alignment-record 16))

(defun mate-reference-id (alignment-record)
  (decode-int32le alignment-record 20))

(defun mate-alignment-position (alignment-record)
  (decode-int32le alignment-record 24))

(defun insert-length (alignment-record)
  (decode-int32le alignment-record 28))

(defun read-name (alignment-record)
  (decode-read-name alignment-record (read-name-length alignment-record)))

(defun alignment-cigar (alignment-record)
  (let* ((name-len (read-name-length alignment-record))
         (cigar-index (+ 32 name-len))
         (cigar-bytes (* 4 (cigar-length alignment-record))))
    (decode-cigar alignment-record cigar-index cigar-bytes)))

(defun seq-string (alignment-record)
  (multiple-value-bind (read-len cigar-index cigar-bytes seq-index
                        seq-bytes qual-index tag-index)
      (alignment-indices alignment-record)
    (declare (ignore cigar-index cigar-bytes seq-bytes qual-index tag-index))
    (decode-seq-string alignment-record seq-index read-len)))

(defun qual-string (alignment-record)
  (multiple-value-bind (read-len cigar-index cigar-bytes seq-index
                        seq-bytes qual-index tag-index)
      (alignment-indices alignment-record)
    (declare (ignore cigar-index cigar-bytes seq-bytes seq-index tag-index))
    (decode-qual-string alignment-record qual-index read-len)))

(defun alignment-tag-values (alignment-record)
  (multiple-value-bind (read-len cigar-index cigar-bytes seq-index
                        seq-bytes qual-index tag-index)
      (alignment-indices alignment-record)
    (declare (ignore read-len cigar-index cigar-bytes seq-bytes qual-index
                     seq-index))
    (decode-tag-values alignment-record tag-index)))

(defun alignment-core (alignment-record)
  (list (reference-id alignment-record)
        (alignment-position alignment-record)
        (read-name-length alignment-record)
        (mapping-quality alignment-record)
        (alignment-bin alignment-record)
        (cigar-length alignment-record)
        (alignment-flag alignment-record)
        (read-length alignment-record)
        (mate-reference-id alignment-record)
        (mate-alignment-position alignment-record)
        (insert-length alignment-record)))

(defun alignment-core-alist (alignment-record)
  (pairlis '(:reference-id :alignment-pos :read-name-length
             :mapping-quality :alignment-bin :cigar-length
             :alignment-flag :read-length :mate-reference-id
             :mate-alignment-position :insert-length)
           (alignment-core alignment-record)))

(defun alignment-indices (alignment-record)
  (let* ((read-len (read-length alignment-record))
         (name-len (read-name-length alignment-record))
         (cigar-index (+ 32 name-len))
         (cigar-bytes (* 4 (cigar-length alignment-record)))
         (seq-index (+ cigar-index cigar-bytes))
         (seq-bytes (ceiling read-len 2))
         (qual-index (+ seq-index seq-bytes))
         (tag-index (+ qual-index read-len)))
    (values read-len cigar-index cigar-bytes seq-index seq-bytes
            qual-index tag-index)))

(defun decode-read-name (alignment-record name-len)
  "Returns a string containing the template/read name of length
NAME-LEN, encoded at byte 32 in ALIGNMENT-RECORD."
  ;; The read name is null terminated and the terminator is included
  ;; in the name length
  (make-sb-string alignment-record 32 (- (+ 32 name-len) 2)))

(defun decode-seq-string (alignment-record index read-len)
  "Returns a string containing the alignment query sequence of length
READ-LEN. The sequence must be present in ALIGNMENT-RECORD at INDEX."
  (declare (optimize (speed 3)))
  (declare (type (simple-array (unsigned-byte 8) (*)) alignment-record)
           (type fixnum index read-len))
  (flet ((decode-base (nibble)
           (ecase nibble
             (0 #\=)
             (1 #\A)
             (2 #\C)
             (4 #\G)
             (8 #\T)
             (15 #\N))))
    (loop
       with seq = (make-array read-len :element-type 'base-char)
       for i from 0 below read-len
       for j of-type fixnum = (+ index (floor i 2))
       do (setf (char seq i)
                (decode-base (if (evenp i)
                                 (ldb (byte 4 4) (aref alignment-record j))
                               (ldb (byte 4 0) (aref alignment-record j)))))
       finally (return seq))))

(defun decode-qual-string (alignment-record index read-len)
  "Returns a string containing the alignment query sequence of length
READ-LEN. The sequence must be present in ALIGNMENT-RECORD at
INDEX. The SAM spec states that quality data are optional, with
absence indicated by 0xff. If the first byte of quality data is 0xff,
NIL is returned."
  (flet ((encode-phred (x)
           (code-char (+ 33 (min 93 x)))))
    (if (= #xff (aref alignment-record index))
        nil
      (loop
         with str = (make-array read-len :element-type 'base-char)
         for i from 0 below read-len
         for j from index below (+ index read-len)
         do (setf (char str i)
                  (encode-phred (decode-uint8le alignment-record j)))
         finally (return str)))))

(defun decode-cigar (alignment-record index cigar-len)
  "Returns a list of CIGAR operations from ALIGNMENT-RECORD at INDEX."
  (flet ((decode-len (uint32)
           (ash uint32 -4))
         (decode-op (uint32)
           (ecase (ldb (byte 4 0) uint32)
             (0 :m)
             (1 :i)
             (2 :d)
             (3 :n)
             (4 :s)
             (5 :h)
             (6 :p))))
    (loop
       for i from index below (1- (+ index cigar-len)) by 4
       collect (decode-uint32le alignment-record i) into cigar
       finally (return (mapcar (lambda (x)
                                 (list (decode-op x) (decode-len x))) cigar)))))

(defun decode-tag-values (alignment-record index)
  "Returns an alist of auxilliary data from ALIGNMENT-RECORD at
INDEX. The BAM two-letter data keys are transformed to Lisp keywords."
  (declare (optimize (speed 3)))
  (declare (type (simple-array (unsigned-byte 8)) alignment-record)
           (type fixnum index))
  (loop
     with tag-index of-type fixnum = index
     while (< tag-index (length alignment-record))
     collect (let* ((type-index (+ tag-index +tag-size+))
                    (type-code (code-char (aref alignment-record type-index)))
                    (val-index (1+ type-index))
                    (tag (make-sb-string alignment-record tag-index
                                         (1+ tag-index)))
                    (val (ecase type-code
                           (#\A         ; A printable character
                            (setf tag-index (+ val-index 1))
                            (code-char (aref alignment-record val-index)))
                           (#\C         ; C unsigned 8-bit integer
                            (setf tag-index (+ val-index 1))
                            (decode-uint8le alignment-record val-index))
                           ((#\H #\Z) ; H hex string, Z printable string
                            (let ((end (position +null-byte+ alignment-record
                                                 :start val-index)))
                              (setf tag-index (1+ end))
                              (make-sb-string alignment-record val-index
                                              (1- end))))
                           (#\I         ; I unsigned 32-bit integer
                            (setf tag-index (+ val-index 4))
                            (decode-uint32le alignment-record val-index))
                           (#\S         ; S unsigned short
                            (setf tag-index (+ val-index 2))
                            (decode-int16le alignment-record val-index))
                           (#\c         ; c signed 8-bit integer
                            (setf tag-index (+ val-index 1))
                            (decode-int8le alignment-record val-index))
                           (#\f         ; f single-precision float
                            (decode-float32le alignment-record val-index))
                           (#\i         ; i signed 32-bit integer
                            (setf tag-index (+ val-index 4))
                            (decode-int32le alignment-record val-index))
                           (#\s         ; s signed short
                            (setf tag-index (+ val-index 2))
                            (decode-int16le alignment-record val-index)))))
               (list tag val))))

(defun sequenced-pair-p (flags)
  (logbitp 0 flags))

(defun mapped-proper-pair-p (flags)
  (logbitp 1 flags))

(defun query-unmapped-p (flags)
  (logbitp 2 flags))

(defun mate-unmapped-p (flags)
  (logbitp 3 flags))

(defun query-forward-p (flags)
  (not (logbitp 4 flags)))

(defun mate-forward-p (flags)
  (not (logbitp 5 flags)))

(defun first-in-pair-p (flags)
  (logbitp 6 flags))

(defun second-in-pair-p (flags)
  (logbitp 7 flags))

(defun alignment-not-primary-p (flags)
  (logbitp 8 flags))

(defun pcr/optical-duplicate-p (flags)
  (logbitp 9 flags))


