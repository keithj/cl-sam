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
        (let ((bin
               (make-bin :num bin-num
                         :chunks (loop
                                    with chunks = (make-array num-chunks)
                                    for i from 0 below num-chunks
                                    do (setf (svref chunks i) (read-chunk))
                                    finally (return chunks)))))
                  (when (= 37450 bin-num)
                    (warn "bin out of range ~a" bin))
                  bin)))))

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
