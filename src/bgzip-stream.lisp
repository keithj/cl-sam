;;;
;;; Copyright (C) 2009-2010 Genome Research Ltd. All rights reserved.
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

(defconstant +bgzip-buffer-size+ 4096
  "Buffer size for {defclass bgzip-input-stream} internal buffer.")

(deftype bgzip-buffer ()
  "Buffer type for {defclass bgzip-input-stream} internal buffer."
  '(simple-array octet (4096)))

(deftype bgzip-buffer-index ()
  "Index type for {defclass bgzf-input-stream} internal buffer."
  '(integer 0 4096))

(defclass bgzf-handle-mixin ()
  ((bgzf :initform nil
         :initarg :bgzf
         :documentation "The BGZF file handle.")))

(defclass bgzip-stream (fundamental-binary-stream bgzf-handle-mixin)
  ()
  (:documentation "A BGZF stream capable of reading or writing block
compressed data."))

(defclass bgzip-input-stream (bgzip-stream fundamental-binary-input-stream)
  ((buffer :initarg :buffer
           :initform nil
           :documentation "The Lisp buffer from which data are read.")
   (num-bytes :initform 0
              :documentation "The number of bytes that were read into
the buffer from the stream.")
   (offset :initform 0
           :documentation "The offset in the byte buffer from which
the next byte is to be read."))
  (:documentation "A stream that reads from a BGZF file."))

(defmacro with-open-bgzip ((var filespec &rest args) &body body)
  `(let ((,var (apply #'bgzip-open ,filespec ,args)))
     (unwind-protect
          (progn
            ,@body)
       (when ,var
         (close ,var)))))

(defun bgzip-open (filespec &key (direction :input))
  "Opens a block gzip stream for FILESPEC.

Key:

- direction (symbol): One of :input (the default) or :output

Returns:

- A {defclass bgzip-stream}"
  (ecase direction
    (:input
     (make-instance 'bgzip-input-stream
                    :bgzf (bgzf-open filespec :direction direction)
                    :buffer (make-array +bgzip-buffer-size+ :element-type 'octet
                                        :initial-element 0)))
    (:output (error "BGZF output streams are not implemented yet."))))

(defmethod stream-element-type ((stream bgzip-stream))
  'octet)

(defmethod stream-close ((stream bgzip-stream) &key abort)
  (declare (ignore abort))
  (when (open-stream-p stream)
    (unwind-protect
         (with-slots (bgzf)
             stream
           (if (bgzf-close bgzf)
               t
               (error 'bgzf-io-error :bgzf bgzf
                      :text "failed to close file cleanly")))
      (call-next-method))))

(defmethod stream-file-position ((stream bgzip-input-stream) &optional position)
  (with-slots (bgzf offset num-bytes)
      stream
    (flet ((num-bytes-buffered ()
             (- num-bytes offset)))
      (cond (position
             (unless (bgzf-seek bgzf position)
               (error 'bgzf-io-error :bgzf bgzf :text "failed to seek in file"))
             (setf offset 0
                   num-bytes 0)
             t)
            (t
             (let ((position (bgzf-tell bgzf)))
               (- position (num-bytes-buffered))))))))

(defmethod stream-read-byte ((stream bgzip-input-stream))
  (if (and (buffer-empty-p stream) (zerop (fill-buffer stream)))
      :eof
    (with-slots (buffer offset)
        stream
      (prog1
          (aref buffer offset)
        (incf offset)))))

#+(or :sbcl :ccl)
(defmethod stream-read-sequence ((stream bgzip-input-stream) sequence
                                 &optional (start 0) end)
  (%stream-read-sequence stream sequence start end))

#+lispworks
(defmethod stream-read-sequence ((stream bgzip-input-stream) sequence
                                 start end)
  (%stream-read-sequence stream sequence start end))

(defun buffer-empty-p (stream)
  (declare (optimize (speed 3) (safety 1)))
  (with-slots (offset num-bytes)
      stream
    (= (the fixnum offset) (the fixnum num-bytes))))

(defun fill-buffer (stream)
  (with-slots (bgzf buffer offset num-bytes)
      stream
    (declare (optimize (speed 3) (safety 1)))
    (declare (type bgzip-buffer buffer))
    (let ((n (length buffer)))
      (multiple-value-bind (buffer num-read)
          (read-bytes bgzf n :buffer buffer)
        (declare (ignore buffer))
        (declare (type bgzip-buffer-index num-read))
        (setf offset 0
              num-bytes num-read)))))

;; (declaim (inline %stream-read-sequence))
(defun %stream-read-sequence (stream sequence &optional (start 0) end)
   ;; (declare (optimize (speed 3) (safety 1)))
  (macrolet ((define-copy-op (seq-type seq-accessor
                              &key (speed 1) (safety 2))
               `(let ((seq-index start))
                  (declare (optimize (speed ,speed) (safety ,safety)))
                  (declare (type ,seq-type sequence)
                           (type fixnum seq-index))
                    (let ((end (or end (length sequence))))
                      (declare (type fixnum end))
                      (loop
                         while (and (not (buffer-empty-p stream))
                                    (< seq-index end))
                         do (loop
                               for i of-type fixnum from seq-index below end
                               for j of-type fixnum from offset below num-bytes
                               do (progn
                                    (setf (,seq-accessor sequence i)
                                          (aref buffer j))
                                    (incf seq-index)
                                    (incf offset))
                               finally (when (buffer-empty-p stream)
                                         (fill-buffer stream)))
                         finally (return seq-index))))))
    (if (and (buffer-empty-p stream) (zerop (the fixnum (fill-buffer stream))))
        0
        (with-slots (buffer offset num-bytes)
            stream
          (declare (type bgzip-buffer buffer)
                   (type bgzip-buffer-index offset num-bytes)
                   (type fixnum start))
          (typecase sequence
            (simple-octet-vector
             (define-copy-op simple-octet-vector aref
               :speed 3 :safety 0))
            (simple-vector
             (define-copy-op simple-vector svref
               :speed 3 :safety 0))
            ((simple-array * (*))
             (define-copy-op (simple-array * (*)) aref))
            (list
             (define-copy-op list elt
               :speed 3 :safety 0))
            (t
             (define-copy-op sequence elt)))))))
