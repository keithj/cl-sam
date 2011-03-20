;;;
;;; Copyright (c) 2010-2011 Keith James. All rights reserved.
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

(defun write-index-magic (stream)
  (write-sequence *bam-index-magic* stream))

(defgeneric write-bam-index (index stream)
  (:method ((index bam-index) stream)
    (let ((bytes (make-array 8 :element-type 'octet :initial-element 0)))
      (write-index-magic stream)
      (%write-int32 (length (bam-index-refs index)) bytes stream)
      (+ 4 4 (reduce #'+ (map 'vector (lambda (ref)
                                        (write-ref-index ref stream))
                              (bam-index-refs index))))))
  (:method ((index samtools-bam-index) stream)
    (let ((bytes (make-array 8 :element-type 'octet :initial-element 0)))
      (write-index-magic stream)
      (%write-int32 (length (bam-index-refs index)) bytes stream)
      (+ 4 4 (reduce #'+ (map 'vector (lambda (ref)
                                        (write-ref-index ref stream))
                              (bam-index-refs index)))
         (length (%write-int64
                  (samtools-bam-index-unassigned index) bytes stream))))))

(defgeneric write-ref-index (ref-index stream)
  (:method ((ref-index ref-index) stream)
    (let ((bytes (make-array 4 :element-type 'octet :initial-element 0)))
      (%write-int32 (length (ref-index-bins ref-index)) bytes stream)
      (+ 4 (write-binning-index (ref-index-bins ref-index) stream)
         (write-linear-index (ref-index-intervals ref-index) stream))))
  (:method ((ref-index samtools-ref-index) stream)
    (let ((bytes (make-array 8 :element-type 'octet :initial-element 0)))
      (%write-int32 (1+ (length (ref-index-bins ref-index))) bytes stream)
      (+ 4 (write-binning-index (ref-index-bins ref-index) stream)
         (write-bin
          (let ((chunk1 (make-chunk
                         :start (samtools-ref-index-start ref-index)
                         :end (samtools-ref-index-end ref-index)))
                (chunk2 (make-chunk
                         :start (samtools-ref-index-mapped ref-index)
                         :end (samtools-ref-index-unmapped ref-index))))
            (make-bin :num +samtools-kludge-bin+
                      :chunks (make-array 2 :initial-contents
                                          (list chunk1 chunk2)))) stream)
         (write-linear-index (ref-index-intervals ref-index) stream)))))

(defun write-binning-index (bins stream)
  (reduce #'+ (map 'vector (lambda (bin)
                             (write-bin bin stream)) bins)))

(defun write-bin (bin stream)
  (let ((bytes (make-array 8 :element-type 'octet))
        (num-chunks (length (bin-chunks bin))))
    (flet ((write-chunk (chunk)
             (%write-int64 (chunk-start chunk) bytes stream)
             (%write-int64 (chunk-end chunk) bytes stream )))
      (%write-int32 (bin-num bin) bytes stream)
      (%write-int32 num-chunks bytes stream)
      (map nil #'write-chunk (bin-chunks bin))
      (+ 4 4 (* 16 num-chunks)))))

(defun write-linear-index (intervals stream)
  (let ((bytes (make-array 8 :element-type 'octet))
        (num-intervals (length intervals)))
    (flet ((write-interval (interval)
             (%write-int64 interval bytes stream)))
      (%write-int32 num-intervals bytes stream)
      (map nil #'write-interval intervals))
    (+ 4 (* 8 num-intervals))))

(defun %write-int32 (n buffer stream)
  (write-sequence (encode-int32le n buffer) stream :end 4))

(defun %write-int64 (n buffer stream)
  (write-sequence (encode-int64le n buffer) stream))
