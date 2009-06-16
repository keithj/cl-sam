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

(defconstant +null-byte+ #x00
  "The termination byte for BAM strings.")
(defconstant +tag-size+ 2
  "The size of a BAM auxilliary tag in bytes.")

(defun reference-id (alignment-record)
  (decode-int32le alignment-record 0))

(defun alignment-position (alignment-record)
  (1+ (decode-int32le alignment-record 4)))

(defun read-name-length (alignment-record)
  (decode-uint8le alignment-record 8))

(defun mapping-quality (alignment-record)
  (decode-uint8le alignment-record 9))

(defun alignment-bin (alignment-record)
  (decode-uint16le alignment-record 10))

(defun cigar-length (alignment-record)
  (decode-uint16le alignment-record 12))

(defun alignment-flag (alignment-record &key (validate t))
  (let ((flag (decode-uint16le alignment-record 14)))
    (when (and validate (not (valid-flag-p flag)))
      (error 'dxi:malformed-field-error
             :record (alignment-core alignment-record :validate nil)
             :field flag
             :text "pair-specific flags set for an unpaired read"))
    flag))

(defun sequenced-pair-p (flag)
  (logbitp 0 flag))

(defun mapped-proper-pair-p (flag)
  (logbitp 1 flag))

(defun query-unmapped-p (flag)
  (logbitp 2 flag))

(defun mate-unmapped-p (flag)
  (logbitp 3 flag))

(defun query-forward-p (flag)
  (not (logbitp 4 flag)))

(defun mate-forward-p (flag)
  (not (logbitp 5 flag)))

(defun first-in-pair-p (flag)
  (logbitp 6 flag))

(defun second-in-pair-p (flag)
  (logbitp 7 flag))

(defun alignment-not-primary-p (flag)
  (logbitp 8 flag))

(defun fails-platform-qc-p (flag)
  (logbitp 9 flag))

(defun pcr/optical-duplicate-p (flag)
  (logbitp 10 flag))

(defun valid-flag-p (flag)
  (if (not (sequenced-pair-p flag))
      ;; If unpaired, the pair data bits should not be set
      (and (not (mapped-proper-pair-p flag))
           (not (mate-unmapped-p flag))
           (not (mate-forward-p flag))
           (not (first-in-pair-p flag))
           (not (second-in-pair-p flag)))
    ;; If paired, must be a proper pair
    (or (and (first-in-pair-p flag) (not (second-in-pair-p flag)))
        (and (second-in-pair-p flag) (not (first-in-pair-p flag))))))

(defun read-length (alignment-record)
  (decode-int32le alignment-record 16))

(defun mate-reference-id (alignment-record)
  (decode-int32le alignment-record 20))

(defun mate-alignment-position (alignment-record)
  (1+ (decode-int32le alignment-record 24)))

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

(defun alignment-core (alignment-record &key (validate t))
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
  (pairlis '(:reference-id :alignment-pos :read-name-length
             :mapping-quality :alignment-bin :cigar-length
             :alignment-flag :read-length :mate-reference-id
             :mate-alignment-position :insert-length)
           (alignment-core alignment-record :validate validate)))

(defun alignment-flag-alist (alignment-record &key (validate t))
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
