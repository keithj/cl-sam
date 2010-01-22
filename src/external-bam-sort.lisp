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

(defclass bgzf-handle-mixin ()
  ((bgzf :initform nil
         :initarg :bgzf
         :accessor bgzf-of)))

(defclass bam-sort-input-stream (sort-input-stream bgzf-handle-mixin)
  ())

(defclass bam-sort-output-stream (sort-output-stream bgzf-handle-mixin)
  ())

(defclass bam-merge-stream (merge-stream wrapped-stream-mixin)
  ())

(defmethod stream-read-element ((stream bam-sort-input-stream))
  (read-alignment (bgzf-of stream)))

(defmethod stream-write-element ((alignment-record vector)
                                 (stream bam-sort-output-stream))
  (write-alignment (bgzf-of stream) alignment-record))

(defmethod make-merge-stream ((stream bam-sort-input-stream) predicate
                              &key key (buffer-size 100000))
  (declare (optimize (speed 3) (safety 1)))
  (declare (type vector-index buffer-size))
  (let ((alignments (loop
                       with buf = (make-array buffer-size :element-type t
                                              :initial-element nil)
                       for i from 0 below buffer-size
                       for alignment = (stream-read-element stream)
                       while alignment
                       do (setf (svref buf i) alignment)
                       finally (return
                                 (if (= i buffer-size)
                                     buf
                                   (subseq buf 0 i))))))
    (cond ((plusp (length alignments))
           (let ((alignments (stable-sort alignments predicate :key key)))
             (loop
                with out = (open (make-tmp-pathname :basename "bam-merge-sort")
                                 :direction :io :element-type 'octet)
                with alen-bytes = (make-array 4 :element-type 'octet
                                              :initial-element 0)
                for i from 0 below (length alignments)
                do (let ((alignment (svref alignments i)))
                     (declare (type simple-octet-vector alignment))
                     (encode-int32le (length alignment) alen-bytes)
                     (write-sequence alen-bytes out)
                     (write-sequence alignment out))
                finally (cond ((file-position out 0)
                               (return (make-instance 'bam-merge-stream
                                                      :stream out)))
                              (t
                               (error 'file-error :pathname out))))))
          (t
           nil))))

(defmethod initialize-instance :after ((stream bam-merge-stream) &key)
  (with-accessors ((s stream-of) (h stream-head-of))
      stream
    (setf h (%read-bam-alignment s))))

(defmethod stream-merge ((stream bam-merge-stream))
  (with-accessors ((s stream-of) (h stream-head-of))
      stream
    (setf h (%read-bam-alignment s))))

(declaim (inline alignment-record<))
(defun alignment-record< (alignment-record1 alignment-record2)
  "Returns T if ALIGNMENT-RECORD1 sorts before
ALIGNMENT-RECORD2. Sorting semantics are not fully defined in the SAM
spec, however, an informal consensus on sequence order sorting seems
to be:

- mapped reads should first be sorted by their reference in the order
  in which reference sequences appear in the header
- unmapped reads should sort after mapped reads
- reads mapped to the same reference must appear in ascending order of
  their alignment position

This function compares first by reference sequence, then alignment
position, by alignment strand and finally by read name. The output is
identical to a coordinate sort performed by Picard 1.07."
  (declare (optimize (speed 3) (safety 1)))
  (let ((ref1 (reference-id alignment-record1))
        (ref2 (reference-id alignment-record2)))
    (declare (type int32 ref1 ref2))
    (cond ((= -1 ref1 ref2)             ; unmapped reads
           nil)
          ((= -1 ref1)                  ; unmapped read
           nil)
          ((= -1 ref2)                  ; unmapped read
           t)
          ((= ref1 ref2)
           (let ((pos1 (alignment-position alignment-record1))
                 (pos2 (alignment-position alignment-record2)))
             (declare (type int32 pos1 pos2))
             (or (< pos1 pos2)
                 (and (= pos1 pos2)
                      (let ((q1 (query-forward-p
                                 (alignment-flag alignment-record1)))
                            (q2 (query-forward-p
                                 (alignment-flag alignment-record2))))
                        (if (eq q1 q2)
                            (let ((name1 (read-name alignment-record1))
                                  (name2 (read-name alignment-record2)))
                              (declare (type simple-base-string name1 name2))
                              (string< name1 name2))
                          (and q1 (not q2))))))))
          (t
           (< ref1 ref2)))))

(defun alignment-name< (alignment-record1 alignment-record2)
  "Returns T if ALIGNMENT-RECORD1 sorts before ALIGNMENT-RECORD2 by
read name. Sorting semantics of read names are not fully defined in
the SAM spec; whether you should expect natural order or numeric order
of read name strings is not defined.

This function compares numerically by read name, then by alignment
template region for reads paired in sequencing and finally by
alignment strand."
  (declare (optimize (speed 3)))
  (let ((name1 (read-name alignment-record1))
        (name2 (read-name alignment-record2)))
    (declare (type simple-base-string name1 name2))
    (or (string< name1 name2)
        (and (string= name1 name2)
             (let ((flag1 (alignment-flag alignment-record1))
                   (flag2 (alignment-flag alignment-record2)))
               (or (and (sequenced-pair-p flag1)
                        (first-in-pair-p flag1) (second-in-pair-p flag2))
                   (and (query-forward-p flag1) (query-reverse-p flag2))))))))

(defun alignment-name-natural< (alignment-record1 alignment-record2)
  (declare (optimize (speed 3) (safety 0)))
  (let* ((len1 (read-name-length alignment-record1))
         (len2 (read-name-length alignment-record2))
         (start 32)
         (end1 (1- (+ start len1)))
         (end2 (1- (+ start len2))))
    (do* ((i start (1+ i))
          (c1 (code-char (aref alignment-record1 i)))
          (c2 (code-char (aref alignment-record2 i))))
         ((< i (min end1 end2)))
      (declare (type fixnum i))
      (cond ((and (digit-char-p c1) (digit-char-p c2))
             (multiple-value-bind (n1 j1)
                 (parse-digits alignment-record1 i end1)
               (multiple-value-bind (n2 j2)
                   (parse-digits alignment-record2 i end2)
                 (declare (type fixnum n1 n2 j1 j2))
                 (if (< n1 n2)
                     (return t)
                   (setf i (min j1 j2))))))
            ((char< c1 c2)
             (return t))
            ((char> c1 c2)
             (return nil))
            (t
             nil)))))

(defun sort-bam-file (in-filespec out-filespec
                      &key (sort-order :coordinate) (buffer-size 1000000))
  "Sorts a BAM file by coordinate or by read name.

Arguments:

- in-filespec (pathname designator): The input BAM file.
- out-filespec (pathname designator): The output BAM file.

Key:

- sort-order (symbol): The sort order, either :coordinate
  or :queryname .

- buffer-size (fixnum): The maximum number of reads to sort in memory
  at any one time, defaulting to 1000000.

Returns:

- The number of alignments sorted.
- The number of files used in the external merge sort."
  (with-bgzf (bgzf-in (pathstring in-filespec) :direction :input)
    (with-bgzf (bgzf-out (pathstring out-filespec) :direction :output
                         :if-exists :supersede)
      (multiple-value-bind (header num-refs ref-meta)
          (read-bam-meta bgzf-in)
        (let ((predicate (ecase sort-order
                           (:coordinate #'alignment-record<)
                           (:queryname #'alignment-name<)))
              (header-str (make-header-string
                           (subst-sam-version
                            (subst-sort-order
                             (ensure-order (make-sam-header header) :sort)
                             sort-order)))))
          (write-bam-meta bgzf-out header-str num-refs ref-meta)
          (sort-bam-alignments bgzf-in bgzf-out predicate
                               :buffer-size buffer-size))))))

(defun sort-bam-alignments (bgzf-in bgzf-out predicate
                            &key key (buffer-size 1000000))
  "Sorts alignments from block gzip input stream BGZF-IN and writes
them to block gzip output stream BGZF-OUT, sorted by PREDICATE. A
function KEY may be supplied to transform alignments into arguments
for PREDICATE. The BUFFER-SIZE argument declares the maximum number of
alignments that will be sorted in memory at any time, defaulting to
1000000."
  (let ((sort-in (make-instance 'bam-sort-input-stream :bgzf bgzf-in))
        (sort-out (make-instance 'bam-sort-output-stream :bgzf bgzf-out)))
    (external-merge-sort sort-in sort-out predicate
                         :key key :buffer-size buffer-size)))

(declaim (inline %read-bam-alignment))
(defun %read-bam-alignment (stream)
  (declare (optimize (speed 3)))
  (let ((alen-bytes (make-array 4 :element-type 'octet :initial-element 0)))
    (if (zerop (read-sequence alen-bytes stream))
        nil
      (let ((record-length (decode-int32le alen-bytes)))
        (if (minusp record-length)
            (error 'malformed-record-error
                   :text "BAM record reported a negative record length")
          (let ((record (make-array record-length :element-type 'octet
                                    :initial-element 0)))
            (read-sequence record stream)
            record))))))

(let ((buffer (make-array 100 :element-type 'base-char :initial-element #\Nul)))
  (defun parse-digits (bytes start end)
    (declare (optimize (speed 3)))
    (declare (type simple-octet-vector bytes)
             (type vector-index start end))
    (let ((len (- end start)))
      (when (> len (length buffer))
        (setf buffer (make-array len :element-type 'base-char
                                 :initial-element #\Nul)))
      (loop
         for i from start below end
         for j = 0 then (1+ j)
         do (setf (char buffer j) (code-char (aref bytes i)))
         finally (return (parse-integer buffer :end j))))))
