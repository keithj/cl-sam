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

(defun view-sam (bam-filespec sam-filespec)
  (with-bgzf-file (bgzf bam-filespec)
    (with-open-file (out sam-filespec
                         :direction :output :element-type 'base-char
                         :external-format :ascii
                         :if-exists :supersede)
      (stream-view-sam bgzf out))))

(defun stream-view-sam (bgzf &optional (stream t))
  (unless (zerop (bgzf-tell bgzf))    ; Rewind unless already at start
    (bgzf-seek bgzf 0))
  (multiple-value-bind (header num-refs ref-meta)
      (read-bam-meta bgzf)
    (declare (ignore num-refs))
    (write-sam-header (ensure-order
                       (ensure-sam-version
                        (make-sam-header header)) :sort) stream)
    (loop
       with refs = (make-reference-table ref-meta)
       for alignment = (read-alignment bgzf)
       while alignment
       count alignment into total
       do (write-sam-alignment alignment refs stream)
       finally (return total))))

(defun flagstat (bam-filespec)
  (with-bgzf-file (bgzf bam-filespec :direction :input)
    (read-bam-meta bgzf)
    (loop
       for alignment = (read-alignment bgzf)
       while alignment
       for flag = (alignment-flag alignment)
       count flag into total
       count (fails-platform-qc-p flag) into qc-fail
       count (pcr/optical-duplicate-p flag) into duplicates
       count (query-mapped-p flag) into mapped
       count (sequenced-pair-p flag) into seq-paired
       count (first-in-pair-p flag) into read1
       count (second-in-pair-p flag) into read2
       count (mapped-proper-pair-p flag) into proper-paired
       count (and (query-mapped-p flag)
                  (mate-mapped-p flag)) into both-mapped
       count (mate-unmapped-p flag) into singletons
       finally (format t #.(str "~d in total~%"
                                "~d QC failure~%"
                                "~d duplicates~%" 
                                "~d mapped (~$%)~%"
                                "~d paired in sequencing~%" 
                                "~d read1~%"
                                "~d read2~%" 
                                "~d properly paired (~$%)~%"
                                "~d both mapped~%" 
                                "~d singletons (~$%)~%")
                       total qc-fail duplicates mapped
                       (* 100 (/ mapped total))
                       seq-paired read1 read2 proper-paired
                       (* 100 (/ proper-paired total))
                       both-mapped singletons
                       (* 100 (/ singletons total))))))

