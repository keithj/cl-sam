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

(defconstant +max-num-bins+ (/ (1- (ash 1 18)) 7)
  "The maximum number of bins possible in the BAM indexing scheme.")
(defconstant +linear-bin-size+ (expt 2 14))

(defvar *bam-index-magic* (make-array 4 :element-type 'octet
                                      :initial-contents '(66 65 73 1))
  "The BAI index file magic header bytes.")

(defun read-index-magic (stream)
  "Reads the BAI magic number from STREAM and returns T if it is valid
or raises a {define-condition malformed-file-error} if not."
  (let ((magic (make-array 4 :element-type 'octet :initial-element 0)))
    (read-sequence magic stream)
    (or (equalp *bam-index-magic* magic)
        (error 'malformed-file-error
               :format-control "invalid BAI magic number ~a"
               :format-arguments (list magic)))))

(defstruct bam-index
  "A BAM file index covering all reference sequences."
  (refs (make-array 0) :type simple-vector))

(defstruct ref-index
  "An index for a single reference sequence."
  (bins (make-array 0) :type simple-vector)
  (intervals (make-array 0) :type simple-vector))

(defstruct (bin (:print-object print-bin))
  "A bin within a BAM binning index."
  (num 0 :type fixnum :read-only t)
  (chunks (make-array 0) :type simple-vector))

(defstruct (chunk (:print-object print-chunk))
  "A chunk within a BAM index bin."
  (start 0 :type (unsigned-byte 64))
  (end 0 :type (unsigned-byte 64)))

(defstruct interval
  "A compressed interval within a BAM linear index."
  (num 1 :type fixnum)
  (offset 0 :type fixnum))

(defun read-bam-index (stream)
  "Reads a BAM (.bai) index from STREAM."
  (when (read-index-magic stream)
    (let ((bytes (make-array 4 :element-type 'octet :initial-element 0)))
      (read-sequence bytes stream)
      (let ((num-refs (decode-int32le bytes)))
        (loop
           with refs = (make-array num-refs)
           for i from 0 below num-refs
           do (setf (svref refs i) (read-ref-index stream))
           finally (return (make-bam-index :refs refs)))))))

(defun read-ref-index (stream)
  "Reads an index for a single reference sequence from STREAM."
  (let ((bytes (make-array 4 :element-type 'octet)))
    (flet ((read-index-size ()
             (read-sequence bytes stream)
             (decode-int32le bytes)))
      (let* ((num-bins (read-index-size))
             (bins (read-binning-index num-bins stream))
             (num-intervals (read-index-size))
             (intervals (read-linear-index num-intervals stream)))
        (make-ref-index :bins bins :intervals intervals)))))

(defun read-binning-index (num-bins stream)
  "Reads NUM-BINS bins from STREAM, returning a vector."
  (loop
     with bins = (make-array num-bins)
     for i from 0 below num-bins
     do (setf (svref bins i) (read-bin stream))
     finally (return bins)))

(defun read-bin (stream)
  "Reads an index bin from STREAM."
  (let ((bytes (make-array 8 :element-type 'octet)))
    (labels ((read-bin-num ()
               (read-sequence bytes stream :end 4)
               (decode-uint32le bytes))
             (read-num-chunks ()
               (read-sequence bytes stream :end 4)
               (decode-int32le bytes))
             (read-voffset ()
               (read-sequence bytes stream)
               (decode-uint64le bytes))
             (read-chunk ()
               (let* ((start (read-voffset))
                      (end (read-voffset)))
                 (make-chunk :start start :end end))))
      (let* ((bin-num (read-bin-num))
             (num-chunks (read-num-chunks)))
        (make-bin :num bin-num
                  :chunks (loop
                             with chunks = (make-array num-chunks)
                             for i from 0 below num-chunks
                             do (setf (svref chunks i) (read-chunk))
                             finally (return chunks)))))))

(defun read-linear-index (num-intervals stream)
  "Reads NUM-INTERVALS linear bin intervals from STREAM, returning
them in a vector. Identical intervals are run-length compressed."
  (let ((bytes (make-array 8 :element-type 'octet)))
    (flet ((read-ioffset ()
             (read-sequence bytes stream)
             (decode-uint64le bytes))
           (compress (ioffsets)
             (when ioffsets
               (let ((current (first ioffsets))
                     (run-length 1)
                     (result ()))
                 (dolist (ioffset (rest ioffsets)
                          (progn
                            (push (list run-length current) result)
                            (nreverse result)))
                     (cond ((= current ioffset)
                            (incf run-length))
                           (t
                            (push (list run-length current) result)
                            (setf current ioffset
                                  run-length 1))))))))
      (let ((intervals (compress (loop
                                    repeat num-intervals
                                    collect (read-ioffset)))))
        (make-array (length intervals)
                    :initial-contents
                    (mapcar (lambda (i)
                              (make-interval :num (first i) :offset (second i)))
                            intervals))))))

(defun print-bin (bin stream)
  "Prints a string representation of BIN to STREAM."
  (print-unreadable-object (bin stream :type t)
    (princ (bin-num bin) stream)
    (princ #\Space stream)
    (princ (bin-chunks bin) stream)))

(defun print-chunk (chunk stream)
  "Prints a string representation of CHUNK to STREAM."
  (print-unreadable-object (chunk stream :type t)
    (let ((start (chunk-start chunk))
          (end (chunk-end chunk)))
      (format stream "~d:~d..~d:~d"
              (bgzf-virtual-position start) (bgzf-virtual-offset start)
              (bgzf-virtual-position end) (bgzf-virtual-offset end)))))

(defun region-to-bins (start end)
  "Returns a bit vector with the bits for relevant bins set, given
a 0-based range.

Each bin may span 2^29, 2^26, 2^23, 2^20, 2^17 or 2^14 bp. Bin 0
spans a 512Mbp region, bins 1-8 span 64Mbp, 9-72 8Mbp, 73-584 1Mbp,
585-4680 128Kbp and bins 4681-37449 span 16Kbp regions."
  (check-arguments (<= 0 start end) (start end) "must satisfy (<= 0 start end)")
  ;; Offsets calculated thus:
  ;; (mapcar (lambda (x)
  ;;           (/ (1- (ash 1 x)) 7)) '(15 12 9 6 3))
  (let ((bins (make-array +max-num-bins+ :element-type 'bit
                          :initial-element 0)))
    (when (< start end)
      (let ((end (1- end)))
        (loop
           for  shift in '( -14 -17 -20 -23 -26)
           for offset in '(4681 585  73   9   1)
           do (loop
                 for k from (+ offset (ash start -14))
                 to (+ offset (ash end shift))
                 do (setf (sbit bins k) 1)))))
    bins))

(defun region-to-bin (start end)
  "Returns a bin number given a 0-based range."
  (check-arguments (<= 0 start end) (start end) "must satisfy (<= 0 start end)")
  (macrolet ((define-bins (offsets shifts)
               `(cond ,@(mapcar (lambda (offset shift)
                                  `((= (ash start ,shift)
                                       (ash end ,shift))
                                    (+ ,(/ (1- (ash 1 offset)) 7)
                                       (ash end ,shift)))) offsets shifts)
                      (t
                       0))))
    (let ((end (1- end)))
      (define-bins
          ( 15  12   9   6   3)
          (-14 -17 -20 -23 -26)))))
