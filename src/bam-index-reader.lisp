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

;;;
;;; Implementation notes:
;;;
;;; The binning index
;;;
;;; The largest bin number for a standard index should be
;;; 37449. However, samtools creates an undocumented extra bin 37450
;;; for each reference containing two chunks which are not really
;;; chunks, but a data kludge where the start and end fields have the
;;; meaning described in the SAM spec:
;;;
;;; First chunk:
;;; "start" field contains the file offset of the start of the reference.
;;;   "end" field contains the file offset of the end of the reference.
;;;
;;; Second chunk:
;;;  "start" field contains the number of mapped reads on the reference
;;;   "end" field contains the number of unmapped reads on the reference
;;;
;;; Sometimes unmapped reads are assigned a reference and mapping
;;; position for sorting purposes. Also note that the extra bin IS NOT
;;; ALWAYS PRESENT.
;;;
;;; In addition, there are 8 bytes appended to the end of the index
;;; containing the number of unmapped reads that have no reference
;;; and (or?) no coordinates. Your guess is as good as mine.

(defun read-index-magic (stream)
  "Reads the BAI magic number from STREAM and returns T if it is valid
or raises a {define-condition malformed-file-error} if not."
  (let ((magic (make-array 4 :element-type 'octet :initial-element 0)))
    (read-sequence magic stream)
    (or (equalp *bam-index-magic* magic)
        (error 'malformed-file-error
               :format-control "invalid BAI magic number ~a"
               :format-arguments (list magic)))))

(defun read-bam-index (stream)
  "Reads a BAM (.bai) index from STREAM."
  (when (read-index-magic stream)
    (let ((bytes (make-array 4 :element-type 'octet :initial-element 0)))
      (read-sequence bytes stream)
      (let ((num-refs (decode-int32le bytes)))
        (loop
           with refs = (make-array num-refs)
           for i from 0 below num-refs
           do (setf (svref refs i) (read-ref-index i stream))
           finally (return (make-bam-index :refs refs)))))))

(defun read-ref-index (ref-num stream)
  "Reads an index for a single reference sequence from STREAM."
  (let ((bytes (make-array 4 :element-type 'octet)))
    (flet ((read-index-size ()
             (read-sequence bytes stream)
             (decode-int32le bytes)))
      (let* ((num-bins (read-index-size))
             (bins (read-binning-index num-bins stream))
             (num-intervals (read-index-size))
             (intervals (read-linear-index num-intervals stream)))
        (if (zerop num-bins)
            (make-ref-index :num ref-num :bins bins :intervals intervals)
            (let ((last-bin (svref bins (1- num-bins)))) ; kludge
              (cond ((= +samtools-kludge-bin+ (bin-num last-bin))
                     (let ((x (svref (bin-chunks last-bin) 0))
                           (y (svref (bin-chunks last-bin) 1)))
                       (make-samtools-ref-index :bins (subseq bins 0
                                                              (1- num-bins))
                                                :intervals intervals
                                                :start (chunk-start x)
                                                :end (chunk-end x)
                                                :mapped (chunk-start y)
                                                :unmapped (chunk-end y))))
                    ((> (bin-num last-bin) +max-num-bins+)
                     (error 'malformed-record-error
                            :record last-bin
                            :format-control "bin number out of range"))
                    (t
                     (make-ref-index :num ref-num :bins bins
                                     :intervals intervals)))))))))

(defun read-binning-index (num-bins stream)
  "Reads NUM-BINS bins from STREAM, returning a vector of bins, sorted
by increasing bin number."
  (loop
     with bins = (make-array num-bins)
     for i from 0 below num-bins
     do (setf (svref bins i) (read-bin stream))
     finally (return (sort bins #'< :key #'bin-num))))

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
