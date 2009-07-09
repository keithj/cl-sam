;;;
;;; Copyright (C) 2009 Keith James. All rights reserved.
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

(defconstant +null-byte+ #x00
  "The termination byte for BAM strings.")
(defconstant +tag-size+ 2
  "The size of a BAM auxilliary tag in bytes.")

(defun reference-id (alignment-record)
  "Returns the reference sequence identifier of ALIGNMENT-RECORD. This
is an integer locally assigned to a reference sequence within the
context of a BAM file."
  (decode-int32le alignment-record 0))

(defun alignment-position (alignment-record)
  "Returns the 1-based sequence coordinate of ALIGNMENT-RECORD in the
reference sequence of the first base of the clipped read."
  (1+ (decode-int32le alignment-record 4)))

(defun read-name-length (alignment-record)
  "Returns the length in ASCII characters of the read name of
ALIGNMENT-RECORD."
  (decode-uint8le alignment-record 8))

(defun mapping-quality (alignment-record)
  "Returns the integer mapping quality of ALIGNMENT-RECORD."
  (decode-uint8le alignment-record 9))

(defun alignment-bin (alignment-record)
  "Returns an integer that indicates the alignment bin to which
ALIGNMENT-RECORD has been assigned."
  (decode-uint16le alignment-record 10))

(defun cigar-length (alignment-record)
  "Returns the number of CIGAR operations in ALIGNMENT-RECORD."
  (decode-uint16le alignment-record 12))

(defun alignment-flag (alignment-record &key (validate t))
  "Returns an integer whose bits are flags that describe properties of
the ALIGNMENT-RECORD. If the VALIDATE key is T (the default) the
flag's bits are checked for internal consistency."
  (let ((flag (decode-uint16le alignment-record 14)))
    (when (and validate (not (valid-flag-p flag)))
      (error 'dxi:malformed-field-error
             :record (alignment-core alignment-record :validate nil)
             :field flag
             :text (format nil (txt "pair-specific flag ~b"
                                    "set for an unpaired read ~s at ~a"
                                    "in reference ~d")
                           flag (read-name alignment-record)
                           (alignment-position alignment-record)
                           (reference-id alignment-record))))
    flag))

(defun sequenced-pair-p (flag)
  "Returns T if FLAG indicates that the read was sequenced as a member
of a pair, or NIL otherwise."
  (logbitp 0 flag))

(defun mapped-proper-pair-p (flag)
  "Returns T if FLAG indicates that the read was mapped as a member of
a properly oriented read-pair, or NIL otherwise."
  (logbitp 1 flag))

(defun query-unmapped-p (flag)
  "Returns T if FLAG indicates that the read was not mapped to a
reference, or NIL otherwise."
  (logbitp 2 flag))

(defun mate-unmapped-p (flag)
  "Returns T if FLAG indicates that the read's mate was not mapped to
a reference, or NIL otherwise."
  (logbitp 3 flag))

(defun query-forward-p (flag)
  "Returns T if FLAG indicates that the read was mapped to the forward
strand of a reference, or NIL if it was mapped to the reverse strand."
  (not (logbitp 4 flag)))

(defun mate-forward-p (flag)
  "Returns T if FLAG indicates that the read's mate was mapped to the
forward, or NIL if it was mapped to the reverse strand."
  (not (logbitp 5 flag)))

(defun first-in-pair-p (flag)
  "Returns T if FLAG indicates that the read was the first in a pair
of reads from one template, or NIL otherwise."
  (logbitp 6 flag))

(defun second-in-pair-p (flag)
  "Returns T if FLAG indicates that the read was the second in a pair
of reads from one template, or NIL otherwise."
  (logbitp 7 flag))

(defun alignment-not-primary-p (flag)
  "Returns T if FLAG indicates that the read mapping was not the
primary mapping to a reference, or NIL otherwise."
  (logbitp 8 flag))

(defun fails-platform-qc-p (flag)
  "Returns T if FLAG indicates that the read failed plaform quality
control, or NIL otherwise."
  (logbitp 9 flag))

(defun pcr/optical-duplicate-p (flag)
  "Returns T if FLAG indicates that the read is a PCR or optical
duplicate, or NIL otherwise."
  (logbitp 10 flag))

(defun valid-flag-p (flag)
  "Returns T if the paired-read-specific bits of FLAG are internally
consistent."
  (if (not (sequenced-pair-p flag))
      ;; If unpaired, the pair data bits should not be set
      (not (or (mapped-proper-pair-p flag)
               (mate-unmapped-p flag)
               (first-in-pair-p flag)
               (second-in-pair-p flag)))
    ;; If paired, must be a proper pair
    (or (and (first-in-pair-p flag) (not (second-in-pair-p flag)))
        (and (second-in-pair-p flag) (not (first-in-pair-p flag))))))

(defun read-length (alignment-record)
  "Returns the length of the read described by ALIGNMENT-RECORD."
  (decode-int32le alignment-record 16))

(defun mate-reference-id (alignment-record)
  "Returns the integer reference ID of ALIGNMENT-RECORD."
  (decode-int32le alignment-record 20))

(defun mate-alignment-position (alignment-record)
  "Returns the 1-based sequence position of the read mate's alignment
described by ALIGNMENT-RECORD."
  (1+ (decode-int32le alignment-record 24)))

(defun insert-length (alignment-record)
  "Returns the insert length described by ALIGNMENT-RECORD."
  (decode-int32le alignment-record 28))

(defun read-name (alignment-record)
  "Returns the read name string described by ALIGNMENT-RECORD."
  (decode-read-name alignment-record (read-name-length alignment-record)))

(defun alignment-cigar (alignment-record)
  "Returns the CIGAR record list of the alignment described by
ALIGNMENT-RECORD. CIGAR operations are given as a list, each member
being a list of a CIGAR operation keyword and an integer operation
length."
  (let* ((name-len (read-name-length alignment-record))
         (cigar-index (+ 32 name-len))
         (cigar-bytes (* 4 (cigar-length alignment-record))))
    (decode-cigar alignment-record cigar-index cigar-bytes)))

(defun seq-string (alignment-record)
  "Returns the sequence string described by ALIGNMENT-RECORD."
  (multiple-value-bind (read-len cigar-index cigar-bytes seq-index
                        seq-bytes qual-index tag-index)
      (alignment-indices alignment-record)
    (declare (ignore cigar-index cigar-bytes seq-bytes qual-index tag-index))
    (decode-seq-string alignment-record seq-index read-len)))

(defun quality-string (alignment-record)
  "Returns the sequence quality string described by ALIGNMENT-RECORD."
  (multiple-value-bind (read-len cigar-index cigar-bytes seq-index
                        seq-bytes qual-index tag-index)
      (alignment-indices alignment-record)
    (declare (ignore cigar-index cigar-bytes seq-bytes seq-index tag-index))
    (decode-quality-string alignment-record qual-index read-len)))

(defun alignment-tag-values (alignment-record)
  "Returns an alist of tag and values described by ALIGNMENT-RECORD."
  (multiple-value-bind (read-len cigar-index cigar-bytes seq-index
                        seq-bytes qual-index tag-index)
      (alignment-indices alignment-record)
    (declare (ignore read-len cigar-index cigar-bytes seq-bytes qual-index
                     seq-index))
    (decode-tag-values alignment-record tag-index)))

(defun alignment-core (alignment-record &key (validate t))
  "Returns a list of the core data described by ALIGNMENT-RECORD. The
list elements are comprised of reference-id, alignment-position,
read-name length, mapping-quality alignment-bin, cigar length,
alignment flag, read length, mate reference-id, mate
alignment-position and insert length."
  (list (reference-id alignment-record)
        (alignment-position alignment-record)
        (read-name-length alignment-record)
        (mapping-quality alignment-record)
        (alignment-bin alignment-record)
        (cigar-length alignment-record)
        (alignment-flag alignment-record :validate validate)
        (read-length alignment-record)
        (mate-reference-id alignment-record)
        (mate-alignment-position alignment-record)
        (insert-length alignment-record)))

(defun alignment-core-alist (alignment-record &key (validate t))
  "Returns the same data as {defun alignment-core} in the form of an
alist."
  (pairlis '(:reference-id :alignment-pos :read-name-length
             :mapping-quality :alignment-bin :cigar-length
             :alignment-flag :read-length :mate-reference-id
             :mate-alignment-position :insert-length)
           (alignment-core alignment-record :validate validate)))

(defun alignment-flag-alist (alignment-record &key (validate t))
  "Returns the bitwise flags of ALIGNMENT-RECORD in the form of an
alist. The primary purpose of this function is debugging."
  (let ((flag (alignment-flag alignment-record :validate validate)))
    (pairlis '(:sequenced-pair :mapped-proper-pair :query-unmapped
               :mate-unmapped :query-forward :mate-forward :first-in-pair
               :second-in-pair :alignment-not-primary :fails-platform-qc
               :pcr/optical-duplicate)
             (list (sequenced-pair-p flag)
                   (mapped-proper-pair-p flag)
                   (query-unmapped-p flag)
                   (mate-unmapped-p flag)
                   (query-forward-p flag)
                   (mate-forward-p flag)
                   (first-in-pair-p flag)
                   (second-in-pair-p flag)
                   (alignment-not-primary-p flag)
                   (fails-platform-qc-p flag)
                   (pcr/optical-duplicate-p flag)))))

(defun alignment-indices (alignment-record)
  "Returns 7 integer values which are byte-offsets within
ALIGNMENT-RECORD at which the various core data lie. See the SAM
spec."
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
           (type (unsigned-byte 32) index read-len))
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
       for j of-type (unsigned-byte 32) = (+ index (floor i 2))
       do (setf (char seq i)
                (decode-base (if (evenp i)
                                 (ldb (byte 4 4) (aref alignment-record j))
                               (ldb (byte 4 0) (aref alignment-record j)))))
       finally (return seq))))

(defun decode-quality-string (alignment-record index read-len)
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
  (declare (optimize (speed 3) (safety 0)))
  (declare (type (simple-array (unsigned-byte 8)) alignment-record)
           (type fixnum index))
  (loop
     with tag-index of-type fixnum = index
     while (< tag-index (length alignment-record))
     collect (let* ((type-index (+ tag-index +tag-size+))
                    (type-code (code-char (aref alignment-record type-index)))
                    (tag (make-sb-string alignment-record tag-index
                                         (1+ tag-index)))
                    (val-index (1+ type-index)))
               (declare (type fixnum val-index))
               (let  ((val (ecase type-code
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
                 (list tag val)))))
