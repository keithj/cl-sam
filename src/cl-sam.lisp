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
  (let ((*print-pretty* nil))
    (with-bgzf (bam bam-filespec)
      (with-open-file (sam sam-filespec
                           :direction :output :element-type 'base-char
                           :external-format :ascii :if-exists :supersede)
        (stream-view-sam bam sam)))))

(defun stream-view-sam (bam &optional (stream t))
  (unless (zerop (bgzf-tell bam))     ; Rewind unless already at start
    (bgzf-seek bam 0))
  (multiple-value-bind (header num-refs ref-meta)
      (read-bam-meta bam)
    (declare (ignore num-refs))
    (declare (optimize (speed 3) (safety 0)))
    (write-sam-header (ensure-order
                       (ensure-sam-version
                        (make-sam-header header)) :sort) stream)
    (loop
       with refs = (make-reference-table ref-meta)
       for alignment = (read-alignment bam)
       while alignment
       count alignment into total
       do (write-sam-alignment alignment refs stream)
       finally (return total))))

(defun flagstat (bam-filespec &optional (stream t))
  "Writes SAM flag counts from BAM-FILESPEC to STREAM."
  (with-bgzf (bam bam-filespec :direction :input)
    (read-bam-meta bam)
    (let ((repaired 0))
      (handler-bind
          ((malformed-record-error
            (lambda (c)
              (incf repaired)
              (use-value (repair-mapping-flags (record-of c))))))
        (loop
           for alignment = (read-alignment bam :validate t)
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
           finally (format stream #.(str "~d in total~%"
                                         "~d QC failure~%"
                                         "~d duplicates~%" 
                                         "~d mapped (~$%)~%"
                                         "~d paired in sequencing~%" 
                                         "~d read1~%"
                                         "~d read2~%" 
                                         "~d properly paired (~$%)~%"
                                         "~d both mapped~%" 
                                         "~d singletons (~$%)~%"
                                         "~d with repaired mapping flags~%")
                           total qc-fail duplicates mapped
                           (* 100 (/ mapped total))
                           seq-paired read1 read2 proper-paired
                           (* 100 (/ proper-paired total))
                           both-mapped singletons
                           (* 100 (/ singletons total))
                           repaired))))))
