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

(defparameter *voffset-merge-distance* (expt 2 15)
  "If two index chunks are this number of bytes or closer to each
other, they should be merged.")

(defun index-bam-file (filespec)
  (with-bgzf (bam filespec)
    (multiple-value-bind (header num-refs ref-meta)
        (read-bam-meta bam)
      (declare (ignore header))
      (let* ((ref-indices (make-array num-refs))
             (bam-index (make-bam-index :refs ref-indices)))
        (dotimes (n num-refs)
          (let ((ref-len (second (assocdr n ref-meta))))
            ;; FIXME -- limit the bin vector length to the minimum
            ;; required, instead of using +max-num-bins+
            (setf (aref ref-indices n)
                  (make-ref-index
                   :num n
                   :bins (make-array +max-num-bins+ :initial-element nil)
                   :intervals (make-array (ceiling ref-len +linear-bin-size+)
                                          :initial-element 0)))))
        (loop
           for voffset = (bgzf-tell bam)
           for aln = (read-alignment bam)
           while aln
           do (let ((pos (alignment-position aln)))
                ;; FIXME -- Do something with -1 pos reads?
                (unless (minusp pos)
                  ;; Update binning index
                  (let* ((stored-bin-num (alignment-bin aln))
                         (ref-num (reference-id aln))
                         (len (alignment-reference-length aln))
                         ;; FIXME -- do we want to override the stored
                         ;; bin number sometimes?
                         (bin-num (if (zerop stored-bin-num)
                                      (region-to-bin pos (+ pos len))
                                      stored-bin-num))
                         (ref-index (ref-index bam-index ref-num))
                         (bin (or (ref-index-bin ref-index bin-num)
                                  (setf (ref-index-bin ref-index bin-num)
                                        (make-bin :num bin-num))))
                         (voffset2 (+ 4 (length aln))))
                    (extend-chunks bin voffset voffset2)
                    ;; Update linear index
                    ;; FIXME -- unaligned reads with a pos but
                    ;; alignment length on reference == 0, are
                    ;; nominally length 1
                    (let* ((intervals (ref-index-intervals ref-index))
                           (start (floor pos +linear-bin-size+))
                           (end (if (zerop len)
                                    (1+ start)
                                    (+ start (floor len +linear-bin-size+)))))
                      (loop
                         for i from start to end
                         when (or (zerop (aref intervals i))
                                  (< voffset (aref intervals i)))
                         do (setf (aref intervals i) voffset)))))))
        bam-index))))

(defun voffset-merge-p (voffset1 voffset2)
  "Returns T if BGZF virtual offsets should be merged into a single
range to save disk seeks."
  (let ((coffset1 (bgzf-coffset voffset1))
        (coffset2 (bgzf-coffset voffset2)))
    (or (= coffset1 coffset2)
        (< coffset1 (+ coffset2 *voffset-merge-distance*)))))

(defun extend-chunks (bin voffset ulen)
  "Extends the chunk range in BIN, either by stretching the end chunk,
or by adding a new one, to cover the additional range from BGZF
virtual offset VOFFSET over length ULEN of uncompressed bytes,
returning BIN. ULEN is normally the uncompressed length of an
alignment record being added to BIN."
  (let* ((chunks (bin-chunks bin))
         (clen (length chunks)))
    (if (plusp clen)
        (let ((last (aref chunks (1- clen))))
          (if (voffset-merge-p (chunk-start last) voffset)
              (setf (chunk-end last) (+ voffset ulen))
              (append-chunk bin (make-chunk :start voffset
                                            :end (+ voffset ulen)))))
        (append-chunk bin (make-chunk :start voffset
                                      :end (+ voffset ulen))))))
