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

;;
;; Indexing algorithm:
;;
;; get file position of the alignment start
;;
;; foreach alignment
;;   get alignment position
;;
;;   if alignment is unassigned (ref-id is -1)
;;     count unassigned
;;   else if alignment is unmapped (flag #x04 is set)
;;     count unmapped
;;   else
;;     count mapped
;;
;;   file chunk start = file position before reading
;;   file chunk end = file position after reading
;;
;;   if alignment position is not negative
;;     get the alignment bin number
;;     get the previous chunks for this bin
;;     get file position of the alignment end
;;     calculate file chunk of alignment
;;     if current chunk is close to previous chunk
;;         extend previous chunk to current chunk end
;;       else
;;         add current chunk to list for this bin
;;
;;     get alignment length on reference
;;     calculate intervals of alignment
;;     foreach interval
;;       if chunk start for interval < current chunk start
;;         chunk start for interval = current chunk start
;;
;; Notes:
;;
;; Each time a new reference is encountered during indexing, a
;; ref-index is created from the accumulated index data of the
;; previous reference.

;; Unassigned reads are identified by a reference id of -1.
;; Unassigned reads are located at the end of a sorted BAM file.
;; Unmapped reads may be assigned to a reference, but may have a pos
;; of -1. Unmapped reads with a pos, have an alignment length of zero.
;;

(defun index-bam-file (filespec)
  "Returns a new BAM-INDEX object, given pathname designator FILESPEC."
  (with-bgzf (bam filespec)
    (multiple-value-bind (header num-refs ref-meta)
        (read-bam-meta bam)
      (declare (ignore header))
      (let ((mapped 0)
            (unmapped 0)
            (unassigned 0)
            (indices ()))
        (labels ((ref-len (ref-id)
                   (check-arguments (assoc ref-id ref-meta) (ref-id)
                                    "no such reference in this BAM file")
                   (second (assocdr ref-id ref-meta)))
                 (make-intervals (ref-id)
                   (make-array (ceiling (ref-len ref-id) +linear-bin-size+)
                               :element-type 'fixnum :initial-element 0))
                 (add-ref-index (ref-id ref-start ref-end chunks intervals)
                   (push (%make-ref-index
                          ref-id ref-start ref-end chunks intervals
                          mapped unmapped) indices))
                 (interval-start (pos)
                   (floor pos +linear-bin-size+))
                 (interval-end (pos len)
                   (if (zerop len)
                       (interval-start pos)
                       (floor (1- (+ pos len)) +linear-bin-size+)))
                 (bin-num (stored-bin-num pos len)
                   (if (zerop stored-bin-num)
                       (region-to-bin pos (+ pos len))
                       stored-bin-num)))
          (do* ((file-pos (bgzf-tell bam) (bgzf-tell bam))
                (aln (read-alignment bam) (read-alignment bam))
                (ref-start file-pos)
                (previous-ref 0)
                (previous-pos 0)
                (chunks (make-hash-table))
                (intervals (make-intervals 0))
                (last-interval 0)
                (eof nil))
               (eof (%make-bam-index num-refs indices unassigned))
            (cond ((and (null aln) (zerop unassigned))
                   (add-ref-index previous-ref ref-start (bgzf-tell bam)
                                  chunks
                                  (subseq intervals 0 (1+ last-interval)))
                   (setf eof t))
                  ((null aln)
                   (setf eof t))
                  (t
                   (let ((ref-id (reference-id aln))
                         (pos (alignment-position aln))
                         (flag (alignment-flag aln)))
                     (if (and (= previous-ref ref-id) (> previous-pos pos))
                         (bam-sort-error previous-ref previous-pos ref-id pos
                                         "previous position > current")
                         (setf previous-pos pos))
                     (cond ((and (/= previous-ref ref-id) (minusp ref-id))
                            (add-ref-index previous-ref ref-start file-pos
                                           chunks
                                           (subseq intervals 0
                                                   (1+ last-interval)))
                            (setf previous-ref ref-id))
                           ((and (/= previous-ref ref-id) (minusp previous-ref))
                            (bam-sort-error previous-ref previous-pos ref-id pos
                                            (txt "alignment to reference ~a"
                                                 "after unassigned reads")
                                            ref-id))
                           ((/= previous-ref ref-id)
                            (add-ref-index previous-ref ref-start file-pos
                                           chunks
                                           (subseq intervals 0
                                                   (1+ last-interval)))
                            (clrhash chunks)
                            (setf previous-ref ref-id
                                  ref-start file-pos
                                  intervals (make-intervals ref-id)
                                  last-interval 0
                                  mapped 0
                                  unmapped 0)))
                     (cond ((minusp ref-id)
                            (incf unassigned))
                           ((query-unmapped-p flag)
                            (incf unmapped))
                           (t
                            (incf mapped)))
                     (when (not (minusp pos))
                       (let* ((len (alignment-reference-length aln))
                              (bin-num (bin-num (alignment-bin aln) pos len))
                              (bin-chunks (gethash bin-num chunks))
                              (chunk (first bin-chunks))
                              (cstart file-pos)
                              (cend (bgzf-tell bam))
                              (istart (interval-start pos))
                              (iend (interval-end pos len)))
                         (if (and chunk
                                  (voffset-merge-p (chunk-end chunk) cstart))
                             (setf (chunk-end chunk) cend)
                             (setf (gethash bin-num chunks)
                                   (cons (make-chunk :start cstart
                                                     :end cend) bin-chunks)))
                         (setf last-interval (max last-interval iend))
                         (loop
                            for i from istart to iend
                            when (or (zerop (aref intervals i))
                                     (< cstart (aref intervals i)))
                            do (setf (aref intervals i) cstart)))))))))))))

(defun %make-ref-index (reference-id ref-start ref-end chunks intervals
                        &optional (mapped 0) (unmapped 0))
  (flet ((make-chunks (elts)
           (make-array (length elts) :initial-contents (nreverse elts)))
         (fill-intervals (elts)
           (loop
              with previous = 0
              for i from 0 below (length elts)
              do (if (zerop (aref elts i))
                     (setf (aref elts i) previous)
                     (setf previous (aref elts i)))
              finally (return elts))))
    (let ((bins (loop
                   for key being the hash-keys of chunks using (hash-value val)
                   collect (make-bin :num key :chunks (make-chunks val)))))
      (make-samtools-ref-index :id reference-id
                               :bins (sort (make-array (length bins)
                                                       :initial-contents bins)
                                           #'< :key #'bin-num)
                               :intervals (fill-intervals intervals)
                               :start ref-start :end ref-end
                               :mapped mapped :unmapped unmapped))))

(defun %make-bam-index (num-refs ref-indices &optional (unassigned 0))
  (let* ((empty-indices (loop
                           for i from (length ref-indices) below num-refs
                           collect (make-ref-index :id i)))
         (ivec (replace (replace (make-array num-refs) ref-indices)
                        empty-indices :start1 (length ref-indices))))
    (make-samtools-bam-index :refs (sort ivec #'< :key #'ref-index-id)
                             :unassigned unassigned)))

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
