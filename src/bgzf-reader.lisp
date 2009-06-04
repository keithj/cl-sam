;;;
;;; Copyright (C) 2009 Keith James. All rights reserved.
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

(defstruct bgzf
  (file nil :type t)
  (ptr nil :type t)
  (open-p nil :type t))

(defun bgzf-open (filespec &key (direction :input))
  (let ((ptr (bgzf-ffi:bgzf-open (namestring filespec) (ecase direction
                                                         (:input "r")
                                                         (:output "w")))))
    (when (null-pointer-p ptr)
      (error 'bgzf-io-error :errno unix-ffi:*error-number*
             :text (format nil "failed to open ~a" filespec)))
    (make-bgzf :file filespec :ptr ptr :open-p t)))

(defun bgzf-close (bgzf)
  (when (bgzf-open-p bgzf)
    (if (zerop (bgzf-ffi:bgzf-close (bgzf-ptr bgzf)))
        t
      (error 'bgzf-io-error :errno unix-ffi:*error-number*
             :text (format nil "failed to close ~a cleanly" bgzf)))))

(defun read-bytes (bgzf n)
  "Reads N bytes from the handle BGZF and returns them as a Lisp array
of unsigned-byte 8. If fewer than N bytes are available an a
{define-condition bgzf-read-error} is raised."
  (declare (optimize (speed 3)))
  (declare (type fixnum n))
  (with-foreign-object (array-ptr :unsigned-char n)
    (let ((num-read (bgzf-ffi:bgzf-read (bgzf-ptr bgzf) array-ptr n)))
      (declare (type fixnum num-read))
      (cond ((zerop num-read)
             nil)
            ((< num-read n)
             (error 'bgzf-io-error
                    :text (format nil #.(txt "expected to read ~a bytes"
                                             "but only ~a were available")
                                  n num-read)))
            (t
             (let ((buffer (make-array n :element-type '(unsigned-byte 8))))
               (loop
                  for i from 0 below n
                  do (setf (aref buffer i)
                           (mem-aref array-ptr :unsigned-char i)))
               buffer))))))

(defun read-string (bgzf n &key null-terminated)
  "Reads N characters from handle BGZF and returns them as a Lisp
string. The NULL-TERMINATED keyword is used to indicate whether the C
string is null-terminated so that the terminator may be consumed."
  (let ((len (if null-terminated
                 (1- n)
               n)))
    (let ((bytes (read-bytes bgzf n)))
      (make-sb-string bytes 0 (1- len)))))
