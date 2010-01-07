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

(defun make-header-string (alist)
  "Returns a new SAM header string representing the header data in
ALIST."
 (with-output-to-string (s)
   (write-sam-header alist s)))

(defun write-sam-header (alist &optional (stream t))
  "Writes SAM header ALIST as a string to STREAM."
  (loop
     for rec in alist
     do (write-header-record rec stream)))

(defun write-header-record (header-record &optional (stream t))
  "Writes alist HEADER-RECORD to STREAM."
  (write-char #\@ stream)
  (write-string (symbol-name (first header-record)) stream)
  (write-char #\Tab stream)
  (dolist (x (intersperse (loop
                             for (tag . value) in (rest header-record)
                             collect (format nil "~a:~a" tag
                                             (typecase value
                                               (symbol (string-downcase
                                                        (symbol-name value)))
                                               (t value)))) #\Tab))
    (princ x stream))
  (terpri stream))

(declaim (inline write-sam-alignment))
(defun write-sam-alignment (alignment-record ref-table &optional (stream t))
  "Writes ALIGNMENT-RECORD to STREAM as SAM. REF-TABLE is a
hash-table created by {defun make-reference-table} ."
  (declare (optimize (speed 3) (safety 0)))
  (multiple-value-bind (read-len cigar-index cigar-bytes seq-index
                        qual-index tag-index)
      (alignment-indices alignment-record)
    (let ((position (1+ (the fixnum (alignment-position alignment-record))))
          (mate-position (1+ (the fixnum (mate-alignment-position
                                          alignment-record)))))
      (declare (type fixnum cigar-bytes position mate-position))
      (write-string (read-name alignment-record) stream)
      (write-char #\Tab stream)
      (princ (alignment-flag alignment-record) stream)
      (write-char #\Tab stream)
      (if (zerop position)
          (write-char #\* stream)
        (write-string
         (gethash (reference-id alignment-record) ref-table) stream))
      (write-char #\Tab stream)
      (princ position stream)
      (write-char #\Tab stream)
      (princ (mapping-quality alignment-record) stream)
      (write-char #\Tab stream)
      (cond ((zerop position)
             (write-char #\* stream))
            ((zerop cigar-bytes)
             (write-char #\* stream))
            (t
             (write-cigar alignment-record cigar-index cigar-bytes stream)))
      (write-char #\Tab stream)
      (cond ((zerop mate-position)
             (write-char #\* stream))
            ((= (the fixnum (reference-id alignment-record))
                (the fixnum (mate-reference-id alignment-record)))
             (write-char #\= stream))
            (t
             (write-string
              (gethash (mate-reference-id alignment-record) ref-table) stream)))
      (write-char #\Tab stream)
      (princ mate-position stream)
      (write-char #\Tab stream)
      (princ (insert-length alignment-record) stream)
      (write-char #\Tab stream)
      (write-seq-string alignment-record seq-index read-len stream)
      (write-char #\Tab stream)
      (write-quality-string alignment-record qual-index read-len stream)
      (write-tag-values alignment-record tag-index stream)
      (terpri stream))))

(declaim (inline write-cigar))
(defun write-cigar (alignment-record index num-bytes
                    &optional (stream t))
  "Writes the CIGAR string of ALIGNMENT-RECORD at INDEX as NUM-BYTES
bytes, to STREAM."
  (declare (optimize (speed 3)))
  (declare (type fixnum index num-bytes))
  (flet ((decode-len (uint32)
           (ash uint32 -4))
         (decode-op (uint32)
           (ecase (ldb (byte 4 0) uint32)
             (0 #\M)
             (1 #\I)
             (2 #\D)
             (3 #\N)
             (4 #\S)
             (5 #\H)
             (6 #\P))))
    (loop
       for i from index below (1- (+ index num-bytes)) by 4
       do (let ((x (decode-uint32le alignment-record i)))
            (princ (decode-len x) stream)
            (write-char (decode-op x) stream)))))

(declaim (inline write-seq-string))
(defun write-seq-string (alignment-record index num-bytes
                         &optional (stream t))
  "Writes the sequence string of ALIGNMENT-RECORD at INDEX as
NUM-BYTES bytes, to STREAM."
  (declare (optimize (speed 3)))
  (declare (type (simple-array (unsigned-byte 8)) alignment-record)
           (type fixnum index))
  (flet ((decode-base (nibble)
           (ecase nibble
             (0 #\=)
             (1 #\A)
             (2 #\C)
             (4 #\G)
             (8 #\T)
             (15 #\N))))
    (loop
       for i of-type fixnum from 0 below num-bytes
       for j of-type uint32 = (+ index (floor i 2))
       do (write-char (decode-base
                       (if (evenp i)
                           (ldb (byte 4 4) (aref alignment-record j))
                         (ldb (byte 4 0) (aref alignment-record j))))
                      stream))))

(declaim (inline write-quality-string))
(defun write-quality-string (alignment-record index num-bytes
                             &optional (stream t))
  "Writes the quality string of ALIGNMENT-RECORD at INDEX as
NUM-BYTES, to STREAM."
  (declare (optimize (speed 3)))
  (declare (type simple-octet-vector alignment-record)
           (type vector-index index num-bytes))
  (flet ((encode-phred (x)
           (code-char (+ 33 (min 93 x)))))
    (if (= #xff (aref alignment-record index))
        (write-char #\* stream)
      (loop
         for i from index below (+ index num-bytes)
         do (write-char (encode-phred (decode-uint8le alignment-record i))
                        stream)))))

(declaim (inline write-tag-values))
(defun write-tag-values (alignment-record index
                         &optional (stream t))
  "Writes the auxilliary data of ALIGNMENT-RECORD at INDEX to STREAM."
  (declare (optimize (speed 3) (safety 0)))
  (declare (type simple-octet-vector alignment-record)
           (type vector-index index))
  (loop
     with tag-index of-type fixnum = index
     while (< tag-index (length alignment-record))
     do (let* ((type-index (+ tag-index +tag-size+))
               (type-code (code-char (aref alignment-record type-index)))
               (tag  (make-sb-string alignment-record tag-index
                                     (1+ tag-index)))
               (val-index (1+ type-index)))
          (declare (type fixnum val-index))
          (write-char #\Tab stream)
          (write-string tag stream)
          (ecase type-code
            (#\A                        ; A printable character
             (setf tag-index (+ val-index 1))
             (write-string ":A:" stream)
             (write-char (code-char (aref alignment-record val-index)) stream))
            (#\C                        ; C unsigned 8-bit integer
             (setf tag-index (+ val-index 1))
             (write-string ":i:" stream)
             (princ (decode-uint8le alignment-record val-index) stream))
            ((#\H #\Z)              ; H hex string, Z printable string
             (let ((end (position +null-byte+ alignment-record
                                  :start val-index)))
               (setf tag-index (1+ end))
               (write-char #\: stream)
               (write-char type-code stream)
               (write-char #\: stream)
               (write-string (make-sb-string alignment-record val-index
                                             (1- end)) stream)))
            (#\I                        ; I unsigned 32-bit integer
             (setf tag-index (+ val-index 4))
             (write-string ":i:" stream)
             (princ (decode-uint32le alignment-record val-index) stream))
            (#\S                        ; S unsigned short
             (setf tag-index (+ val-index 2))
             (write-string ":i:" stream)
             (princ (decode-uint16le alignment-record val-index) stream))
            (#\c                        ; c signed 8-bit integer
             (setf tag-index (+ val-index 1))
             (write-string ":i:" stream)
             (princ (decode-int8le alignment-record val-index) stream))
            (#\f                        ; f single-precision float
             (setf tag-index (+ val-index 4))
             (write-string ":f:" stream)
             (princ (decode-float32le alignment-record val-index) stream))
            (#\i                        ; i signed 32-bit integer
             (setf tag-index (+ val-index 4))
             (write-string ":i:" stream)
             (princ (decode-int32le alignment-record val-index) stream))
            (#\s                        ; s signed short
             (setf tag-index (+ val-index 2))
             (write-string ":i:" stream)
             (princ (decode-int16le alignment-record val-index) stream))))))
