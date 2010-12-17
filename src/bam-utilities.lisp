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

(defun repair-mapping-flags (aln)
  "Returns alignment ALN if it has correct mapped-proper-pair,
query-mapped and mate-mapped flags, or a modified copy with fixed
flags. Also sets the user tag ZF:I:<original flag> if the flags are
changed.

Alignments produced by the BWA aligner have been observed to contain
invalid flags. This is caused by BWA creating mappings that overhang
the end of the reference. The underlying cause is reference
concatenation in the Burrows-Wheeler index.

This function attempts to fix invalid mapped-proper-pair flags in
cases where query-unmapped and/or mate-unmapped flags are set. Such
unmapped reads may also have reference-ids set and good mapping
scores. This function does not correct these other fields.

It is necessary to use this function in order to obtain correct flag
counts when using samtools flagstat or {defun flagstat} ."
  (let ((flag (alignment-flag aln)))
    (cond ((not (sequenced-pair-p flag))
           aln)
          ((and (mapped-proper-pair-p flag)
                (or (query-unmapped-p flag) (mate-unmapped-p flag)))
           (copy-alignment-record
            aln :alignment-flag (flag-bits flag :not-mapped-proper-pair)
            :tag-values (acons :zf flag (alignment-tag-values aln))))
          (t
           aln))))

(defun generate-bam-file (filespec &key (num-refs 10) (ref-length 1000)
                          (start 0) end (step-size 10)
                          (read-length 50) (insert-length 250)
                          (mapping-quality 50))
  "Writes a very uniform BAM file to the file denoted by pathname
designator FILESPEC, or testing purposes."
  (let* ((tmp (tmp-pathname))
         (end (or end ref-length))
         (read-group "generated_group")
         (ref-meta (loop
                      for i from 0 below num-refs
                      collect (list i (format nil "ref_~d" i) ref-length)))
         (header (with-output-to-string (s)
                   (write-sam-header
                    (apply #'list
                           (hd-record :version "1.3" :sort-order :coordinate)
                           (rg-record read-group "generated")
                           (loop
                              for m in ref-meta
                              collect (sq-record (second m) (third m)))) s))))
    (with-bam (bam (header num-refs ref-meta) tmp :direction :output
                   :if-exists :supersede)
      (do ((i start (+ step-size i)))
          ((>= (+ i insert-length (* 2 read-length)) end) filespec)
        (let ((seq1 (make-string read-length :initial-element #\a))
              (seq2 (make-string read-length :initial-element #\t))
              (quality (make-string read-length :initial-element #\e))
              (pos1 i)
              (pos2 (+ i read-length insert-length))
              (cigar (list (cons :m read-length))))
          (consume bam (make-alignment-record 
                        (format nil "r~9,'0d" i) seq1
                        (flag-bits 0 :sequenced-pair :first-in-pair
                                   :mapped-proper-pair
                                   :query-forward :mate-reverse
                                   :query-mapped :mate-mapped)
                        :quality-str quality
                        :reference-id 0 :mate-reference-id 0
                        :alignment-pos pos1 :mate-alignment-pos pos2
                        :insert-length insert-length
                        :mapping-quality mapping-quality :cigar cigar
                        :tag-values `((:rg . ,read-group) (:nm . 0))))
          (consume bam (make-alignment-record 
                        (format nil "r~9,'0d" i) seq2
                        (flag-bits 0 :sequenced-pair :second-in-pair
                                   :mapped-proper-pair
                                   :query-reverse :mate-forward
                                   :query-mapped :mate-mapped)
                        :quality-str quality
                        :reference-id 0 :mate-reference-id 0
                        :alignment-pos pos2 :mate-alignment-pos pos1
                        :insert-length insert-length
                        :mapping-quality mapping-quality :cigar cigar
                        :tag-values `((:rg . ,read-group) (:nm . 0)))))))
    (sort-bam-file tmp filespec)
    filespec))
