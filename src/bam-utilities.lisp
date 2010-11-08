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

(defmacro with-bam ((var (&optional header num-refs ref-meta) filespec
                         &rest args &key (compress t) (null-padding 0)
                         &allow-other-keys)
                    &body body)
  "Evaluates BODY with VAR bound to a newly opened BAM iterator on
pathname designator FILESPEC. The direction (:input versus :output) is
determined by the stream-opening arguments in ARGS.

On reading, HEADER, NUM-REFS and REF-META will be automatically bound
to the BAM file metadata i.e. the metadata are read automatically and
the iterator positioned before the first alignment record.

On writing, HEADER, NUM-REFS and REF-META should be bound to
appropriate values for the BAM file metadata, which will be
automatically written to the underlying stream.

The COMPRESS and NULL-PADDING keyword arguments are only applicable on
writing where they control whether the header block should be
compressed and whether the header string should be padded with nulls
to allow space for expansion.

For example:

To count all records in a BAM file:

;;; (with-bam (in (header) \"in.bam\")
;;;   (declare (ignore header))
;;;   (loop
;;;      while (has-more-p in)
;;;      count (next in)))

To copy only records with a mapping quality of >= 30 to another BAM
file:

;;; (with-bam (in (h n r) \"in.bam\")
;;;   (with-bam (out (h n r) \"out.bam\" :direction :output)
;;;     (let ((q30 (discarding-if (lambda (x)
;;;                                 (< (mapping-quality x) 30)) in)))
;;;       (loop
;;;          while (has-more-p q30)
;;;          do (consume out (next q30))))))"
  (with-gensyms (bgzf)
    (let ((default-header (with-output-to-string (s)
                            (write-sam-header
                             `((:HD (:VN . ,*sam-version*))) s))))
      `(with-bgzf (,bgzf ,filespec ,@(remove-key-values
                                      '(:compress :null-padding) args))
         ,@(if (search '(:direction :output) args)
               `((write-bam-meta ,bgzf
                                 ,(or header default-header)
                                 ,(or num-refs 0) ,ref-meta
                                 :compress ,compress
                                 :null-padding ,null-padding)
                 (let ((,var (make-bam-output ,bgzf)))
                   ,@body))
               `((multiple-value-bind (,@(when header `(,header))
                                       ,@(when num-refs `(,num-refs))
                                       ,@(when ref-meta `(,ref-meta)))
                     (read-bam-meta ,bgzf)
                   ,@(cond (ref-meta
                            `((declare (ignorable ,header ,num-refs))))
                           (num-refs
                            `((declare (ignorable ,header)))))
                   (let ((,var (make-bam-input ,bgzf)))
                     ,@body))))))))

(defun make-bam-input (bam)
  "Returns a generator function that returns BAM alignment records for
BAM stream BAM. The standard generator interface functions NEXT and
HAS-MORE-P may be used in operations on the returned generator.

For example:

To count all records in a BAM file:

;;; (with-bgzf (bam \"in.bam\")
;;;   (read-bam-meta bam)
;;;   (let ((in (make-bam-input bam)))
;;;     (loop
;;;        while (has-more-p in)
;;;        count (next in))))

To copy only records with a mapping quality of >= 30 to another BAM
file:

;;; (with-bgzf (bam-in \"in.bam\")
;;;   (with-bgzf (bam-out \"out.bam\" :direction :output)
;;;     (multiple-value-bind (header num-refs ref-meta)
;;;         (read-bam-meta bam-in)
;;;       (write-bam-meta bam-out header num-refs ref-meta))
;;;     (let ((fn (discarding-if (lambda (x)
;;;                                (< (mapping-quality x) 30))
;;;                              (make-bam-input bam-in))))
;;;       (loop
;;;          while (has-more-p fn)
;;;          do (write-alignment bam-out (next fn))))))"
  (let ((current (read-alignment bam)))
    (defgenerator
        (more (not (null current)))
        (next (prog1
                  current
                (setf current (read-alignment bam)))))))

(defun make-bam-output (bam)
  "Returns a consumer function that accepts an argument of a BAM
record and writes it to BAM output stream BAM. The standard consumer
interface function CONSUME may be used in operations on the returned
consumer."
  (lambda (aln)
    (write-alignment bam aln)))

(defun repair-mapping-flags (aln)
  "Returns alignment ALN if it has correct mapped-proper-pair,
query-mapped and mate-mapped flags, or a modified copy with fixed
flags. Also sets the user tag ZF:I:<original flag> if the flags are
changed.

Alignments produced by the BWA have been observed to contain invalid
flags. This is caused by BWA creating mappings that overhang the end
of the reference. The underlying cause is reference concatenation in
the Burrows-Wheeler index.

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
