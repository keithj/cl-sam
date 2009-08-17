;;;
;;; Copyright (C) 2009 Keith James. All rights reserved.
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
  (let ((out (open (make-tmp-pathname :basename "bam-merge-sort")
                   :direction :io
                   :element-type '(unsigned-byte 8)))
        (alignments (make-array buffer-size :adjustable t :fill-pointer 0)))
    (declare (type (vector (simple-array (unsigned-byte 8))) alignments))
    (loop
       for i from 0 below buffer-size
       for alignment = (stream-read-element stream)
       while alignment
       do (vector-push alignment alignments))
    (cond ((plusp (length alignments))
           (loop
              with alen-bytes = (make-array 4 :element-type '(unsigned-byte 8))
              for alignment across (sort alignments predicate :key key)
              do (progn
                   (encode-int32le (length alignment) alen-bytes)
                   (write-sequence alen-bytes out)
                   (write-sequence alignment out))
              finally (if (file-position out 0)
                          (return
                            (make-instance 'bam-merge-stream :stream out))
                        (error 'file-error :pathname out))))
          (t
           (close out :abort t)
           nil))))

(defmethod initialize-instance :after ((stream bam-merge-stream) &key)
  (with-accessors ((s stream-of) (h stream-head-of))
      stream
    (setf h (%read-bam-alignment s))))

(defmethod stream-merge ((stream bam-merge-stream))
  (with-accessors ((s stream-of) (h stream-head-of))
      stream
    (setf h (%read-bam-alignment s))))

(defun alignment-position< (alignment-record1 alignment-record2)
  (< (alignment-position alignment-record1)
     (alignment-position alignment-record2)))

(defun bam-alignment-position< (alignment-record1 alignment-record2)
  (declare (optimize (speed 3)))
  (let ((ref1 (reference-id alignment-record1))
        (ref2 (reference-id alignment-record2)))
    (declare (type int32 ref1 ref2))
    (cond ((= ref1 ref2)
           (alignment-position< alignment-record1 alignment-record2))
          ((minusp ref1)
           nil)
          ((minusp ref2)
           t)
          (t
           t))))

;; FIXME -- edit the header to include proper sort metadata
(defun sort-bam-file (in-filespec out-filespec predicate
                      &key key (buffer-size 5000000))
  (with-bgzf-file (bgzf-in (namestring in-filespec) :direction :input)
    (with-bgzf-file (bgzf-out (namestring out-filespec) :direction :output)
      (multiple-value-bind (header num-refs ref-meta)
          (read-bam-meta bgzf-in)
        (write-bam-meta bgzf-out header num-refs ref-meta)
        (let ((sort-in (make-instance 'bam-sort-input-stream :bgzf bgzf-in))
              (sort-out (make-instance 'bam-sort-output-stream :bgzf bgzf-out)))
          (external-merge-sort sort-in sort-out predicate
                               :key key :buffer-size buffer-size))))))

(defun %read-bam-alignment (stream)
  (let ((alen-bytes (make-array 4 :element-type '(unsigned-byte 8))))
    (if (zerop (read-sequence alen-bytes stream))
        nil
      (let ((record-length (decode-int32le alen-bytes)))
        (if (minusp record-length)
            (error 'malformed-record-error
                   :text "BAM record reported a negative record length")
          (let ((record (make-array record-length
                                    :element-type '(unsigned-byte 8)
                                    :initial-element 0)))
            (copy-array alen-bytes 0 3
                        record 0)
            (read-sequence record stream)
            record))))))
