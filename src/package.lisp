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

(in-package :cl-user)

(defpackage :bgzf-ffi
  (:use #:common-lisp #:cffi)
  (:export #:bgzf-open
           #:bgzf-close
           #:bgzf-read
           #:bgzf-write
           #:bgzf-seek
           #:bgzf-tell
           #:bgzf-check-eof)
  (:documentation "A foreign function interface to the block gzip
  functions provided by the Broad Institute. See bgzf.h in the
  SAMtools package for function documentation."))

(defpackage :sam
  (:use #:common-lisp #:cffi
        #:deoxybyte-utilities #:deoxybyte-io)
  (:export
   ;; Conditions
   #:bgzf-io-error
   #:bam-error
   #:bam-sort-error

   ;; BGZF handle API
   #:bgzf
   #:bzgf-p
   #:bgzf-file
   #:bgzf-open-p
   #:bgzf-open
   #:bgzf-close
   #:bgzf-seek
   #:bgzf-tell
   #:bgzf-eof-p
   #:with-bgzf

   ;; Gray streams API
   #:bgzip-stream
   #:bgzip-input-stream
   #:bgzip-output-stream
   #:bgzip-open
   #:with-open-bgzip
   
   ;; Low-level BAM reading API
   #:read-bam-magic
   #:read-bam-header
   #:read-num-references
   #:read-reference-meta
   #:read-bam-meta
   #:read-alignment
   #:read-bam-terminator

   ;; Low-level BAM writing API
   #:write-bam-magic
   #:write-bam-header
   #:write-num-references
   #:write-reference-meta
   #:write-bam-meta
   #:write-alignment

   ;; BAM IO
   #:with-bam
   #:make-bam-input
   #:make-bam-output

   ;; Low-level BAM API
   #:reference-id
   #:alignment-position
   #:read-name-length
   #:mapping-quality
   #:alignment-bin
   #:cigar-length
   #:alignment-flag
   #:read-length
   #:mate-reference-id
   #:mate-alignment-position
   #:insert-length

   #:read-name
   #:alignment-cigar
   #:seq-string
   #:quality-string
   #:alignment-tag-values

   #:sequenced-pair-p
   #:mapped-proper-pair-p
   #:query-mapped-p
   #:query-unmapped-p
   #:mate-mapped-p
   #:mate-unmapped-p
   #:query-forward-p
   #:query-reverse-p
   #:mate-forward-p
   #:mate-reverse-p
   #:first-in-pair-p
   #:second-in-pair-p
   #:alignment-not-primary-p
   #:fails-platform-qc-p
   #:pcr/optical-duplicate-p

   #:valid-flag-p
   #:valid-mapped-pair-p
   #:valid-mapped-proper-pair-p
   #:valid-pair-num-p

   #:alignment-core
   #:alignment-core-alist
   #:alignment-flag-alist

   #:make-alignment-record
   #:copy-alignment-record
   #:make-reference-table
   #:flag-bits
   #:define-alignment-tag
   #:alignment-tag-documentation

   ;; BAM sorting
   #:sort-bam-file
   #:alignment-record<
   #:alignment-name<
   #:alignment-strand<

   ;; BAM indexing
   #:with-bam-index
   #:read-bam-index
   #:write-bam-index
   #:index-bam-file
   #:find-bins

   ;; BAM index
   #:bam-index
   #:ref-index
   #:bin
   #:chunk
   #:empty-bin-p
   #:bam-index-refs
   #:ref-index-num
   #:ref-index-bins
   #:ref-index-intervals
   #:bin-num
   #:bin-chunks
   #:chunk-start
   #:chunk-end
   #:ref-index-bin
   #:bin-chunk

   ;; Low-level SAM API
   #:*sam-version*
   #:valid-sam-version-p
   #:valid-reference-name-p
   #:make-sam-header
   #:make-header-record
   #:header-records
   #:header-type
   #:header-tags
   #:header-value
   #:ensure-mandatory-header-tags
   #:ensure-valid-header-tags
   #:ensure-valid-programs
   #:previous-programs
   #:last-programs

   #:hd-record
   #:sq-record
   #:rg-record
   #:pg-record
   #:add-pg-record

   #:write-sam-header
   #:write-header-record
   #:write-sam-alignment

   ;; SAM header convenience
   #:name-sorted-p
   #:coordinate-sorted-p

   ;; SAM header editing
   #:merge-sam-headers
   #:merge-header-records
   #:subst-sam-version
   #:subst-sort-order
   #:subst-group-order

   ;; Utilities
   #:view-sam
   #:flagstat
   #:repair-mapping-flags
   #:generate-bam-file
   #:generate-reference-file
   #:alignment-generator)
  (:documentation "cl-sam is a Common Lisp toolkit for manipulation of
DNA sequence alignment data stored in the Sequence Alignment/Map (SAM)
format <http://samtools.sourceforge.net>. The SAM specficiation
describes text (SAM) and binary (BAM) formats.

cl-sam uses the BGZF block compression code from the SAMtools C
toolkit, but is otherwise an independent Lisp implementation capable
of editing every aspect of BAM files or creating BAM data de novo.

The following example implements something similar to samtools
flagstat:

;;; (with-bgzf (bgzf \"example.bam\" :direction :input)
;;;   (multiple-value-bind (header num-refs ref-meta)
;;;       (read-bam-meta bgzf)
;;;     (format t \"BAM header: ~s~%\" header)
;;;     (format t \"Number of references: ~d~%\" num-refs)
;;;     (loop
;;;        for (id name length) in ref-meta
;;;        do (format t \"Reference name: ~s Length: ~d~%\" id name length)))
;;;   (loop
;;;      for alignment = (read-alignment bgzf)
;;;      while alignment
;;;      for flag = (alignment-flag alignment)
;;;      count flag into total
;;;      count (fails-platform-qc-p flag) into qc-fail
;;;      count (pcr/optical-duplicate-p flag) into duplicates
;;;      count (not (query-unmapped-p flag)) into mapped
;;;      count (sequenced-pair-p flag) into seq-paired
;;;      count (first-in-pair-p flag) into read1
;;;      count (second-in-pair-p flag) into read2
;;;      count (mapped-proper-pair-p flag) into proper-paired
;;;      count (and (not (query-unmapped-p flag))
;;;                 (not (mate-unmapped-p flag))) into both-mapped
;;;      count (mate-unmapped-p flag) into singletons
;;;      finally (format t (str \"~d in total~%\"
;;;                             \"~d QC failure~%\"
;;;                             \"~d duplicates~%\"
;;;                             \"~d mapped (~$%)~%\"
;;;                             \"~d paired in sequencing~%\"
;;;                             \"~d read1~%\"
;;;                             \"~d read2~%\"
;;;                             \"~d properly paired (~$%)~%\"
;;;                             \"~d both mapped~%\"
;;;                             \"~d singletons (~$%)~%\")
;;;                      total qc-fail duplicates mapped
;;;                      (* 100 (/ mapped total))
;;;                      seq-paired read1 read2 proper-paired
;;;                      (* 100 (/ proper-paired total))
;;;                      both-mapped singletons
;;;                      (* 100 (/ singletons total)))))"))
