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
(defconstant +linear-bin-size+ (expt 2 14)
  "The size in bases of intervals in the linear bin index.")

(defvar *bam-index-magic* (make-array 4 :element-type 'octet
                                      :initial-contents '(66 65 73 1))
  "The BAI index file magic header bytes.")

(defstruct bam-index
  "A BAM file index covering all reference sequences."
  (refs (make-array 0) :type simple-vector))

(defstruct ref-index
  "An index for a single reference sequence. An index is composed of a
vector of bins (the binning index) and a vector of intervals (the
linear index). The bins are sorted by increasing bin number."
  (bins (make-array 0) :type simple-vector)
  (intervals (make-array 0) :type simple-vector))

(defstruct (bin (:print-object print-bin))
  "A bin within a BAM binning index. A bin is a region on a reference
sequence that is identified by a number. The bin location and
numbering system used is the UCSC binning scheme by Richard Durbin and
Lincoln Stein. A bin contains a vector of chunks."
  (num 0 :type fixnum :read-only t)
  (chunks (make-array 0) :type simple-vector))

(defstruct (chunk (:print-object print-chunk))
  "A chunk within a BAM index bin. Chunks are groups of records that
lie within a bin and are close together within the BAM file."
  (start 0 :type (unsigned-byte 64))
  (end 0 :type (unsigned-byte 64)))

(defstruct interval
  "A compressed interval within a BAM linear index. An interval is a
range covering 2^14 bases that starts a particular virtual offset
within a BAM file. Identical, adjacent intervals are run-length
compressed by incrementing the interval number for each interval in
the run."
  (num 1 :type fixnum)
  (offset 0 :type fixnum))

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

(defun ref-index (bam-index ref-num)
  "Returns the reference index for reference number REF-NUM in
BAM-INDEX."
  (svref (bam-index-refs bam-index) ref-num))

(defun region-to-bins (start end)
  "Returns a bit vector with the bits for relevant bins set, given
a 0-based range.

Each bin may span 2^29, 2^26, 2^23, 2^20, 2^17 or 2^14 bp. Bin 0 spans
a 512 Mbp region, bins 1-8 span 64 Mbp, 9-72 8 Mbp, 73-584 1 Mbp,
585-4680 128 Kbp and bins 4681-37449 span 16 Kbp regions."
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

(defun region-to-interval (coord)
  "Given a coordinate COORD, returns the array index of its interval
in the BAM linear index."
  (floor start +linear-bin-size+))

(defun find-bins (bam-index ref-num start end)
  (let ((bin-nums (region-to-bins start end))
        (bins (ref-index-bins (ref-index bam-index ref-num))))
    (loop
       for i from 0 below +max-num-bins+
       for bin = (when (plusp (sbit bin-nums i))
                   (binary-search bins i :key #'bin-num))
       when bin
       collect bin)))
