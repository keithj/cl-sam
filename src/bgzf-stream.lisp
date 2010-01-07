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

(defconstant +bgzf-buffer-size+ 8192
  "Buffer size for {defclass bgzf-input-stream} internal buffer.")

(deftype bgzf-buffer ()
  "Buffer type for {defclass bgzf-input-stream} internal buffer."
  `(simple-array (unsigned-byte 8) (,+bgzf-buffer-size+)))

(deftype bgzf-buffer-index ()
  "Index type for {defclass bgzf-input-stream} internal buffer."
  `(integer 0 ,+bgzf-buffer-size+))

(defclass bgzf-stream (fundamental-binary-stream)
  ((bgzf :initarg :bgzf
         :reader bgzf-of
         :documentation "The BGZF file handle."))
  (:documentation "A BGZF stream capable of reading or writing block
compressed data."))

(defclass bgzf-input-stream (bgzf-stream fundamental-binary-input-stream)
  ((buffer :initarg :buffer
           :initform nil
           :reader buffer-of
           :documentation "The Lisp buffer from which data are read.")
   (num-bytes :initform 0
              :accessor num-bytes-of
              :documentation "The number of bytes that were read into
the buffer from the stream.")
   (offset :initform 0
           :accessor offset-of
           :documentation "The offset in the byte buffer from which
the next byte is to be read."))
  (:documentation "A stream that reads from a BGZF file."))

(defun bgzf-stream-open (filespec &key (direction :input))
  (ecase direction
    (:input
     (make-instance 'bgzf-input-stream
                    :bgzf (bgzf-open filespec :direction direction)
                    :buffer (make-array +bgzf-buffer-size+ :element-type 'octet
                                        :initial-element 0)))
    (:output (error "BGZF output streams are not implemented yet."))))

(defmethod stream-element-type ((stream bgzf-stream))
  'octet)

(defmethod close ((stream bgzf-stream) &key abort)
  (declare (ignore abort))
  (when (open-stream-p stream)
    (unwind-protect 
         (if (bgzf-close (bgzf-of stream))
             t
           (error 'bgzf-io-error :text "failed to close file cleanly"))
      (call-next-method))))

(defmethod stream-file-position ((stream bgzf-input-stream) &optional position)
  (cond (position
         (unless (bgzf-seek (bgzf-of stream) position)
           (error 'bgzf-io-error :text "failed to seek in file"))
         (setf (offset-of stream) 0
               (num-bytes-of stream) 0)
         t)
        (t
         (let ((position (bgzf-tell (bgzf-of stream))))
           (- position (num-bytes-buffered stream))))))

(defmethod stream-read-byte ((stream bgzf-input-stream))
  (if (and (buffer-empty-p stream) (zerop (fill-buffer stream)))
      :eof
    (with-accessors ((buffer buffer-of) (offset offset-of))
        stream
      (prog1
          (aref buffer offset)
        (incf offset)))))

(defmethod stream-read-sequence ((stream bgzf-input-stream) sequence
                                 &optional (start 0) end)
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
      (with-accessors ((buffer buffer-of) (offset offset-of)
                       (num-bytes num-bytes-of))
          stream
        (declare (type bgzf-buffer buffer)
                 (type bgzf-buffer-index offset num-bytes)
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

(defun buffer-empty-p (stream)
  (declare (optimize (speed 3) (safety 1)))
  (= (the fixnum (offset-of stream)) (the fixnum (num-bytes-of stream))))

(defun num-bytes-buffered (stream)
  (- (num-bytes-of stream) (offset-of stream)))

(defun fill-buffer (stream)
  (with-accessors ((bgzf bgzf-of) (buffer buffer-of) (offset offset-of)
                   (num-bytes num-bytes-of))
      stream
    (declare (optimize (speed 3) (safety 1)))
    (declare (type bgzf-buffer buffer))
    (let ((n (length buffer)))
      (multiple-value-bind (buffer num-read)
          (read-bytes bgzf n :buffer buffer)
        (declare (ignore buffer))
        (declare (type bgzf-buffer-index num-read))
        (setf offset 0
              num-bytes num-read)))))
