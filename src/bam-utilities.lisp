;;;
;;; Copyright (c) 2010-2011 Keith James. All rights reserved.
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

(defun generate-reference-file (filespec name length
                                &optional (fn (lambda ()
                                                (char "acgt" (random 4)))))
  "Creates an artificial reference genome Fasta file to accompany a
generated BAM file.

Arguments:

- filespec (pathname designator): The file to create.
- name (string): The name of the reference sequence.
- length (integer): The length of the reference sequence.

Optional:

- fn (function): A function that returns a single base character. The
  default randomly returns a, c, g or t in equal proportions.

Returns:

- The filespec."
  (with-open-file (ref filespec :direction :output :if-exists :supersede)
    (write-char #\> ref)
    (write-line name ref)
    (loop
       for start from 0 below length by 60
       for end = (min (+ start 60) length)
       do (loop
             repeat (- end start)
             do (write-char (funcall fn) ref)
             finally (terpri ref))))
  filespec)

(defun encode-phred-quality (q)
  "Returns the character encoding Phred quality Q."
  (code-char (+ q 33)))

(defun default-seq-string (length)
  (make-string length :initial-element #\A))

(defun default-qual-string (length)
  (make-string length :initial-element #\e))

(defun alignment-generator (reference-id read-group
                            &key (read-length 50) name-suffix
                            (insert-length 250) (start 0)
                            (end (+ insert-length (* 2 read-length) 1))
                            (step-size 10) (mapping-quality 50)
                            (seq-fn #'default-seq-string)
                            (quality-fn #'default-qual-string))
  "Returns a new standard generator function which returns pairs of
BAM alignment records.

Arguments:

- reference-id (integer): The BAM reference identifier.
- read-group (string): The read group name.

Key:

- read-length (integer): The read length.
- insert-length (integer): The insert length, defined as the distance
  between the inner boundaries of the paired alignments.
- start (integer): The position at which to start generating alignments.
- end (integer): The position at which to stop generating alignments,
  expressed in Lisp vector indices (zero-based, half-open), rather
  than BAM's (zero-based, closed).
- step-size (integer): The distance between each successive pair of
  alignments.
- mapping-quality (integer): The mapping quality of the alignments.

Returns

- A function."
  (let ((i start)
        read1 read2)
    (defgenerator
        (more (< (+ i insert-length (* 2 read-length)) end))
        (next (let ((seq (funcall seq-fn read-length))
                    (quality (funcall quality-fn read-length))
                    (read-name (format nil "r~9,'0d~@[.~a~]" i name-suffix))
                    (pos1 i)
                    (pos2 (1- (+ i read-length insert-length)))
                    (cigar (list (cons :m read-length))))
                (setf read1 (make-alignment-record
                             read-name seq
                             (flag-bits 0 :sequenced-pair :first-in-pair
                                        :mapped-proper-pair
                                        :query-forward :mate-reverse
                                        :query-mapped :mate-mapped)
                             :quality-str quality
                             :reference-id reference-id
                             :mate-reference-id reference-id
                             :alignment-pos pos1 :mate-alignment-pos pos2
                             :insert-length insert-length
                             :mapping-quality mapping-quality :cigar cigar
                             :tag-values `((:rg . ,read-group) (:nm . 0)))
                      read2 (make-alignment-record
                             read-name seq
                             (flag-bits 0 :sequenced-pair :second-in-pair
                                        :mapped-proper-pair
                                        :query-reverse :mate-forward
                                        :query-mapped :mate-mapped)
                             :quality-str quality
                             :reference-id reference-id
                             :mate-reference-id reference-id
                             :alignment-pos pos2 :mate-alignment-pos pos1
                             :insert-length insert-length
                             :mapping-quality mapping-quality :cigar cigar
                             :tag-values `((:rg . ,read-group) (:nm . 0)))
                      i (+ step-size i))
                (values read1 read2))))))

(defun default-orphanizer-test ()
  "Returns a new orphanizer test that causes alternate last fragments
to be omitted."
  (let (orphanize)
    (lambda (aln)
      (let ((flag (alignment-flag aln)))
        (cond ((and orphanize (last-frag-p flag))
               (setf orphanize nil)
               t)
              ((last-frag-p flag)
               (setf orphanize t)
               nil)
              (t
               nil))))))

(defun alignment-orphanizer (gen &key (test (default-orphanizer-test)))
  "Given an alignment generator function GEN, returns a new function
that omits some alignments, yielding orphans. Any alignments for which
TEST returns T will be omitted. The default is to omit alternate last
fragment alignments."
 (defgenerator
     (more (has-more-p gen))
     (next (multiple-value-bind (fwd rev)
               (next gen)
             (values
              (unless (funcall test fwd)
                fwd)
              (unless (funcall test rev)
                rev))))))

(defun generate-bam-file (filespec num-refs ref-length read-groups
                          &rest aln-generators)
  "Writes a very uniform BAM file to the file denoted by pathname
designator FILESPEC, for testing purposes.

Arguments:

- filespec (pathname designator): The BAM file to be written.
- num-refs (integer): The number of reference sequences in the BAM header.
- read-groups (list of string): The read group names in the BAM header. Each
  reference will have a set of reads generated in each read group.

Rest:

- aln-generators (list of function): A list of alignment generator functions
  created with {defun alignment-generator} function. All of these functions
  will be called until exhausted, so they should represent finite sequences.

Returns

- filespec."
  (check-arguments (and (integerp num-refs) (plusp num-refs)) (num-refs)
                   "must be a positive integer")
  (check-arguments (integerp ref-length) (ref-length)
                   "must be an integer")
  (let* ((tmp (tmp-pathname))
         (ref-meta (loop
                      for i from 0 below num-refs
                      collect (list i (format nil "ref_~d" i) ref-length)))
         (header (with-output-to-string (s)
                   (write-sam-header
                    (concatenate 'list
                                 (list (hd-record :version *sam-version*
                                                  :sort-order :coordinate))
                                 (mapcar (lambda (rg)
                                           (rg-record rg "generated"))
                                         read-groups)
                                 (loop
                                    for m in ref-meta
                                    collect (sq-record (second m)
                                                       (third m)))) s))))
    (with-bam (bam (header num-refs ref-meta) tmp :direction :output
                   :if-exists :supersede)
      (dolist (fn aln-generators)
        (loop
           while (has-more-p fn)
           do (multiple-value-bind (fwd rev) ; One may be nil to give orphans
                  (next fn)
                (when fwd
                  (consume bam fwd))
                (when rev
                  (consume bam rev))))))
    (sort-bam-file tmp filespec)
    filespec))
