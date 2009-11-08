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

(declaim (ftype (function (bgzf bam-alignment fixnum) fixnum) write-bytes))

(defun write-bytes (bgzf bytes n)
  (declare (optimize (speed 3)))
  (declare (type (simple-array (unsigned-byte 8) (*)) bytes)
           (type fixnum n))
  (with-foreign-object (array-ptr :unsigned-char n)
    (loop
       for i from 0 below n
       do (setf (mem-aref array-ptr :unsigned-char i) (aref bytes i)))
    (let ((num-written (bgzf-ffi:bgzf-write (bgzf-ptr bgzf) array-ptr n)))
      (declare (type fixnum num-written))
      (if (< num-written n)
          (error 'bgzf-io-error
                 :text (format nil "expected to write ~a bytes but only ~a were written"
                               n num-written))
        n))))
