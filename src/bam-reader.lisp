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

(defmacro with-bam ((var (&optional header num-refs ref-meta) filespec
                         &rest args &key (compress t) (null-padding 0)
                         index ref-num (start 0) (end (expt 2 29))
                         &allow-other-keys)
                    &body body)
  "Evaluates BODY with VAR bound to a new BAM generator function on
pathname designator FILESPEC. The direction (:input versus :output) is
determined by the stream-opening arguments in ARGS. The standard
generator interface functions NEXT and HAS-MORE-P may be used in
operations on the returned generator.

On reading, HEADER, NUM-REFS and REF-META will be automatically bound
to the BAM file metadata i.e. the metadata are read automatically and
the iterator positioned before the first alignment record.

Optionally, iteration may be restricted to a specific reference,
designated by REF-NUM and to alignments that start between reference
positions START and END. Furthermore, a BAM-INDEX object INDEX may be
provided to allow the stream to seek directly to the desired region of
the file.

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
                                      '(:compress :null-padding :index
                                        :ref-num :start :end) args))
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
                   (let ((,var (make-bam-input ,bgzf
                                               :index ,index
                                               :ref-num ,ref-num
                                               :start ,start
                                               :end ,end)))
                     ,@body))))))))

(defun make-bam-input (bam &key index ref-num (start 0) (end (expt 2 29)))
  "Returns a generator function that returns BAM alignment records for
BAM stream BAM. Optionally, iteration may be restricted to a specific
reference, designated by REF-NUM and to alignments that start between
reference positions START and END. Furthermore, a BAM-INDEX object
INDEX may be provided to allow the stream to seek directly to the
desired region of the file. The standard generator interface functions
NEXT and HAS-MORE-P may be used in operations on the returned
generator.

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
  (flet ((out-of-bounds-p (aln)
           (or (not aln)
               (let ((pos (alignment-position aln)))
                 (or (< pos start) (> pos end))))))
    (let ((chunks (and index (make-bam-chunk-queue index ref-num start end))))
      (cond ((and chunks (queue-empty-p chunks))
             (defgenerator
                 (more nil)
                 (next nil)))
            (chunks
             (bgzf-seek bam (max (aref (ref-index-intervals
                                        (ref-index index ref-num))
                                       (region-to-interval start))
                                 (chunk-start (queue-first chunks))))
             (let ((fn (make-bam-chunk-input bam chunks end)))
               (discarding-if #'out-of-bounds-p fn)))
            (t
             (let ((fn (let ((current (read-alignment bam)))
                         (defgenerator
                             (more (not (null current)))
                             (next (prog1
                                       current
                                     (setf current (read-alignment bam))))))))
               (discarding-if #'out-of-bounds-p fn)))))))

(defun read-bam-magic (bgzf)
  "Reads the BAM magic number from the handle BGZF and returns T if it
is valid or raises a {define-condition malformed-file-error} if not."
  (let ((magic (read-bytes bgzf 4)))
    (or (equalp *bam-magic* magic)
        (error 'malformed-file-error :file (bgzf-pathname bgzf)
               :format-control "invalid BAM magic number ~a"
               :format-arguments (list magic)))))

(defun read-bam-header (bgzf)
  "Returns the unparsed BAM header from the handle BGZF as a Lisp
string."
  (let ((len (decode-int32le (read-bytes bgzf 4))))
    (when (plusp len)
      (read-string bgzf len))))

(defun read-num-references (bgzf)
  "Returns the number of reference sequences described in the BAM
header of handle BGZF."
  (decode-int32le (read-bytes bgzf 4)))

(defun read-reference-meta (bgzf)
  "Returns the reference sequence metadata for a single reference
sequence described in the BAM header of handle BGZF. Two values are
returned, the reference name as a string and the reference length as
an integer,"
  (let ((len (decode-int32le (read-bytes bgzf 4))))
    (values (ensure-valid-reference-name (read-string bgzf len
                                                      :null-terminated t))
            (decode-int32le (read-bytes bgzf 4)))))

(defun read-alignment (bgzf &key validate)
  "Reads one alignment block from handle BGZF, returns it as a Lisp
array of unsigned-byte 8. The handle is advanced to the next
alignment. If no more alignments are available, returns NIL.

This is the preferred function to use for validation. i.e. validate
early. A USE-VALUE restart is provided so that an invalid alignment
may be modified and re-read without unwinding the stack."
  (declare (optimize (speed 3)))
  (let ((size-bytes (read-bytes bgzf 4)))
    (if (null size-bytes)
        nil
        (let* ((block-size (decode-int32le size-bytes))
               (aln (read-bytes bgzf block-size)))
          (restart-case
              (if validate
                  (and (alignment-flag aln :validate t) aln)
                  aln)
            (use-value (value) value))))))

(defun read-bam-meta (bgzf)
  "Reads all BAM metadata from handle BGZF, leaving the handle
pointing at the first alignment record. Returns the header string, the
number of references and a list of reference sequence metadata. The
list contains one element per reference sequence, each element being a
list of reference identifier, reference name and reference length."
  (read-bam-magic bgzf)
  (let ((header (read-bam-header bgzf))
        (num-refs (read-num-references bgzf)))
    (values header num-refs
            (loop
               for ref-id from 0 below num-refs
               collect (cons ref-id (multiple-value-list
                                     (read-reference-meta bgzf)))))))

(defun read-bam-terminator (bgzf)
  "Reads the EOF from handle BGZF, returning T if it is present, or
raising a {define-condition malformed-file-error} otherwise."
  (or (bgzf-eof-p bgzf)
      (error 'malformed-file-error :file (bgzf-pathname bgzf)
             :format-control "BGZF EOF was missing")))

