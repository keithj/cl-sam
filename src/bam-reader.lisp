;;;
;;; Copyright (c) 2009-2011 Keith James. All rights reserved.
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

(defun make-bam-full-scan-input (bam)
  "Returns a generator function that returns BAM alignment records for
BAM stream BAM by scanning the entire stream. This to be used in cases
where a full scan is necessary, or as a fallback when neither index,
nor regions are available. The standard generator interface functions
NEXT and HAS-MORE-P may be used in operations on the returned
generator."
  (let ((aln (read-alignment bam)))
    (defgenerator
        (more (not (null aln)))
        (next (prog1
                  aln
                (setf aln (read-alignment bam)))))))

(defun make-bam-scan-input (bam regions)
  "Returns a generator function that returns BAM alignment records for
BAM stream BAM by scanning for alignments in the list of REGIONS. Note
that REGIONS use zero-based, closed base coordinates. This to be used
in cases where a full scan is necessary, or as a fallback when an
index is not available. The standard generator interface functions
NEXT and HAS-MORE-P may be used in operations on the returned
generator."
  (let ((regions (queue-append (make-queue) regions))
        ref-id rstart rend)
    (labels ((next-region ()
             (unless (queue-empty-p regions)
               (let ((region (queue-dequeue regions)))
                 (setf ref-id (region-ref region)
                       rstart (region-start region)
                       rend (region-end region)))))
             (before-region-p (aln)
               (let ((id (reference-id aln)))
               (or (< id ref-id)
                   (and (= id ref-id)
                        (< (alignment-reference-end aln) rstart)))))
             (after-region-p (aln)
               (or (> (reference-id aln) ref-id)
                   (> (alignment-position aln) rend)))
             (find-alignment ()
               (do (aln endp)
                   (endp aln)
                 (let ((next (read-alignment bam)))
                   (cond ((null next)
                          (setf endp t))
                         ((after-region-p next)
                          (next-region))
                         ((before-region-p next)
                          nil)
                         (t
                          (setf aln next
                                endp t)))))))
      (cond ((queue-empty-p regions)
             (defgenerator
                 (more nil)
                 (next nil)))
            (t
             (next-region)
             (let ((aln (find-alignment)))
               (defgenerator
                   (more (not (null aln)))
                   (next (prog1
                             aln
                           (setf aln (find-alignment)))))))))))

(defun make-bam-index-input (bam index regions)
  "Returns a generator function that returns BAM alignment records for
BAM stream BAM by scanning for alignments in the list of REGIONS.
Note that REGIONS are zero-based, closed base coordinates. The
standard generator interface functions NEXT and HAS-MORE-P may be used
in operations on the returned generator."
  (let ((regions (queue-append (make-queue) regions))
        (chunks (make-queue))
        chunk ref-id rstart rend)
    (labels ((next-region ()
               (unless (queue-empty-p regions)
                 (let ((region (queue-dequeue regions)))
                   (setf ref-id (region-ref region)
                         rstart (region-start region)
                         rend (region-end region))
                   (queue-clear chunks)
                   (queue-append chunks (make-bam-chunks index region))
                   (unless (queue-empty-p chunks)
                     ;; FIXME -- When switching region, check whether
                     ;; the next chunk is close enough to the current
                     ;; bgzf position such that a seek is not necessary
                     (setf chunk (queue-dequeue chunks))
                     (bgzf-seek bam (max (find-interval
                                          (ref-index index ref-id) rstart)
                                         (chunk-start chunk)))))))
             (next-chunk ()
               (unless (queue-empty-p chunks)
                 (setf chunk (queue-dequeue chunks))
                 (bgzf-seek bam (chunk-start chunk))))
             (before-region-p (aln)
               (let ((id (reference-id aln)))
                 (or (< id ref-id)
                     (and (= id ref-id)
                          (< (alignment-reference-end aln) rstart)))))
             (after-region-p (aln)
               (or (> (reference-id aln) ref-id)
                   (> (alignment-position aln) rend)))
             (find-alignment ()
               (do (aln endp)
                   (endp aln)
                 (let ((next (read-alignment bam)))
                   (cond ((null next)
                          (setf endp t))
                         ((after-region-p next)
                          (next-region))
                         ((before-region-p next)
                          nil)
                         (t
                          (setf aln next
                                endp t)))
                   (when (and (bgzf-empty-p bam) ; ensure buffer is flushed
                              (>= (bgzf-tell bam) (chunk-end chunk)))
                     (if (queue-empty-p chunks)
                         (next-region)
                         (next-chunk)))))))
      (cond ((queue-empty-p regions)
             (defgenerator
                 (more nil)
                 (next nil)))
            (t
             (next-region)
             (let ((aln (find-alignment)))
               (defgenerator
                   (more (not (null aln)))
                   (next (prog1
                             aln
                           (setf aln (find-alignment)))))))))))

(defun make-region< (ref-meta)
  "Returns a new comparator function which sorts regions according to
the BAM reference metadata REF-META. Ranges are sorted first by the
reference order in the BAM file, then by region start and finally by
region end."
  (let ((lookup (loop
                   with x = (make-array (length ref-meta) :element-type 'fixnum)
                   for elt in ref-meta
                   for i = 0 then (1+ i)
                   do (setf (aref x (first elt)) i)
                   finally (return x))))
    (lambda (r1 r2)
      (or (< (aref lookup (region-ref r1)) (aref lookup (region-ref r2)))
          (< (region-start r1) (region-start r2))
          (and (= (region-start r1) (region-start r2))
               (< (region-end r1) (region-end r2)))))))

(defun normalise-regions (regions ref-meta)
  "Returns a list of REGIONS, normalised according to the BAM
reference metadata REF-META. Reference names in REGIONS are converted
to their corresponding reference ID and REGIONS are sorted first by
the reference order in the BAM file, then by region start and finally
by region end. Overlapping REGIONS are merged."
  (flet ((normalise (regions)
           (let ((name-to-id (let ((x (make-hash-table :test #'equal)))
                               (dolist (elt ref-meta x)
                                 (setf (gethash (second elt) x) (first elt))))))
             (let ((normalised (make-hash-table)))
               (dolist (region (mapcar #'ensure-valid-region regions))
                 (let* ((ref (region-ref region))
                        (ref-num (gethash ref name-to-id ref))
                        (norm (region ref-num (region-start region)
                                      (region-end region))))
                   (check-arguments (integerp ref-num) (regions)
                                    "unknown reference name in ~a" region)
                   (if (gethash ref-num normalised)
                       (push norm (gethash ref-num normalised))
                       (setf (gethash ref-num normalised) (list norm)))))
               (loop
                  for v being the hash-values of normalised nconc v))))
         (merge-overlaps (regions)
           (if (endp (rest regions))
               regions
               (let ((r1 (copy-region (first regions)))
                     merged)
                 (dolist (r2 (rest regions) (nreverse (cons r1 merged)))
                   (cond ((equalp r1 r2)
                          nil)
                         ((> (region-start r2) (region-end r1))
                          (push r1 merged)
                          (setf r1 (copy-region r2)))
                         ((> (region-end r2) (region-end r1))
                          (setf (region-end r1) (region-end r2)))))))))
    (merge-overlaps (normalise regions))))

(defun partition-regions (regions)
  "Partitions a list of normalised REGIONS into lists, each pertaining
to a specific reference sequence. The overall order of the REGIONS is
maintained."
  (let ((current-ref (region-ref (first regions)))
        partitioned current)
    (dolist (region regions (nreverse (cons current partitioned)))
      (cond ((= current-ref (region-ref region))
             (push region current))
            (t
             (push current partitioned)
             (setf current-ref (region-ref region)
                   current (list region)))))))

(defun ensure-valid-region (region)
  "Returns REGION if it is a valid region designator, or raises an
{define-condition invalid-argument-error} .

A region designator may be any of the following:

 - A list of a string and two integers, being a reference name,
   start position and end position, respectively.

 - A list of three integers, being a reference identifier, start
   position and end position, respectively.

 - A region struct.

The start position must be less than, or equal to the end position."
  (check-arguments (or (listp region)
                       (region-p region)) (region)
                       "region designator must be a list or region")
  (flet ((check-region (ref start end)
           (check-arguments (or (stringp ref) (integerp ref)) (ref)
                            "reference designator must be a string or integer")
           (check-arguments (and (integerp start)
                                 (integerp end)
                                 (<= start end)) (start end)
                                 "must satisfy (<= start end)")
           t))
    (etypecase region
      (list (and (apply #'check-region region) (apply #'region region)))
      (region (and (check-region (region-ref region) (region-start region)
                                 (region-end region))
                  region)))))

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
  (let ((size-bytes (read-bytes bgzf +alignment-size-tag+)))
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
               for ref-num from 0 below num-refs
               collect (cons ref-num (multiple-value-list
                                      (read-reference-meta bgzf)))))))

(defun read-bam-terminator (bgzf)
  "Reads the EOF from handle BGZF, returning T if it is present, or
raising a {define-condition malformed-file-error} otherwise."
  (or (bgzf-eof-p bgzf)
      (error 'malformed-file-error :file (bgzf-pathname bgzf)
             :format-control "BGZF EOF was missing")))
