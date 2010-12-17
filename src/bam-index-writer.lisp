;;;
;;; Copyright (C) 2010 Keith James. All rights reserved.
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

(defun index-bam-file (filespec)
  "Returns a new BAM-INDEX object, given pathname designator FILESPEC."
  (with-bgzf (bam filespec)
    (multiple-value-bind (header num-refs ref-meta)
        (read-bam-meta bam)
      (declare (ignore header))
      (let ((chunks (make-hash-table :test #'equal))
            (intervals (make-hash-table)))
        (dotimes (n num-refs)
          (let ((ref-len (second (assocdr n ref-meta))))
            (setf (gethash n intervals) (make-array
                                         (ceiling ref-len +linear-bin-size+)
                                         :element-type 'fixnum
                                         :initial-element 0))))
        (loop
           for cstart = (bgzf-tell bam)
           for aln = (read-alignment bam)
           while aln
           do (let ((pos (alignment-position aln)))
                ;; FIXME -- do something with -1 pos reads?
                ;; FIXME -- do we want to override the stored
                ;; bin number sometimes?
                (unless (minusp pos)
                  (let* ((ref-num (reference-id aln))
                         (len (alignment-reference-length aln))
                         (stored-bin-num (alignment-bin aln))
                         (bin-num (if (zerop stored-bin-num)
                                      (region-to-bin pos (+ pos len))
                                      stored-bin-num))
                         (ckey (list ref-num bin-num))
                         (cend (+ cstart 4 (length aln))))
                    (let* ((cseq (or (gethash ckey chunks)
                                     (setf (gethash ckey chunks)
                                           (make-array 0 :fill-pointer 0
                                                       :adjustable t))))
                           (clen (length cseq))
                           (clast (when (plusp clen)
                                    (aref cseq (1- clen)))))
                      (if (and clast (voffset-merge-p
                                      (chunk-start clast) cstart))
                          (setf (chunk-end clast) cend)
                          (vector-push-extend (make-chunk :start cstart
                                                          :end cend) cseq)))
                    ;; FIXME -- unaligned reads with a pos but
                    ;; alignment length on reference == 0, are
                    ;; nominally length 1
                    (let* ((iseq (gethash ref-num intervals))
                           (istart (floor pos +linear-bin-size+))
                           (iend (if (zerop len)
                                     (1+ istart)
                                     (+ istart (floor len +linear-bin-size+)))))
                      (loop
                         for i from istart to iend
                         when (or (zerop (aref iseq i))
                                  (< cstart (aref iseq i)))
                         do (setf (aref iseq i) cstart)))))))
        (build-bam-index chunks intervals)))))

(defun build-bam-index (chunks intervals)
  "Returns a new bam index given a hash-table CHUNKS of chunk vectors
keyed on ref-num and bin-num and a hash-table INTERVALS of interval
vectors, keyed on ref-num."
  (let* ((num-refs (hash-table-count intervals))
         (bins (make-hash-table :size num-refs))
         (ref-indices (make-array num-refs)))
    (loop
       for (ref-num bin-num) being the hash-keys of chunks
       using (hash-value chunks)
       do (let ((bvec (gethash ref-num bins))
                (bin (make-bin :num bin-num
                               :chunks (make-array (length chunks)
                                                   :initial-contents chunks))))
            (if bvec
                (vector-push-extend bin bvec)
                (setf (gethash ref-num bins)
                      (make-array 1 :fill-pointer 0 :adjustable t
                                  :initial-element bin)))))
    (loop
       for ref-num being the hash-keys of bins
       using (hash-value bvec)
       do (setf (svref ref-indices ref-num)
                (make-ref-index
                 :num ref-num
                 :bins (sort (make-array (length bvec) :initial-contents bvec)
                             #'< :key #'bin-num)
                 :intervals (gethash ref-num intervals))))
    (make-bam-index :refs ref-indices)))

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
