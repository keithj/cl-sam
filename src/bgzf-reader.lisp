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

(defmacro with-bgzf-file ((var filespec &key (direction :input)) &body body)
  "Executes BODY with VAR bound to a BGZF handle structure created by
opening the file denoted by FILESPEC.

Arguments:

- var (symbol): The symbol to be bound.
- filespec (pathname designator): The file to open.

Key:

- direction (keyword): The direction, one of either :read or :write ."
  `(let ((,var (bgzf-open ,filespec :direction ,direction)))
    (unwind-protect
         (progn
           ,@body)
      (when ,var
        (bgzf-close ,var)))))

(defun bgzf-open (filespec &key (direction :input))
  "Opens a block gzip file for reading or writing.

Arguments:

- filespec (pathname designator): The file to open.

Key:

- direction (keyword): The direction, one of either :read or :write .

Returns:

- A BGZF structure."
  (let ((ptr (bgzf-ffi:bgzf-open (namestring filespec) (ecase direction
                                                         (:input "r")
                                                         (:output "w")))))
    (when (null-pointer-p ptr)
      (error 'bgzf-io-error :errno unix-ffi:*error-number*
             :text (format nil "failed to open ~a" filespec)))
    (make-bgzf :file filespec :ptr ptr :open-p t)))

(defun bgzf-close (bgzf)
  "Closes an open block gzip handle.

Arguments:

- bgzf (bgzf structure): The handle to close.

- Returns: T on success."
  (when (bgzf-open-p bgzf)
    (if (zerop (bgzf-ffi:bgzf-close (bgzf-ptr bgzf)))
        t
      (error 'bgzf-io-error :errno unix-ffi:*error-number*
             :text (format nil "failed to close ~a cleanly" bgzf)))))

(defun bgzf-seek (bgzf position)
  "Seeks with the file encapsulated by a block gzip handle.

Arguments:

- bgzf (bgzf structure): The handle to seek.
- position (integer): The position to seek. Only values previously
  returned by {defun bgzf-tell} may be used.

- Returns: The new position."
  (zerop
   (bgzf-ffi:bgzf-seek (bgzf-ptr bgzf) position
                       (foreign-enum-value'unix-ffi:seek-directive :seek-set))))

(defun bgzf-tell (bgzf)
  "Returns the current position in the encapsulated file of a block gzip
handle.

Arguments:

- bgzf (bgzf structure): The handle.

- Returns: The file position."
  (let ((position (bgzf-ffi:bgzf-tell (bgzf-ptr bgzf))))
    (when (minusp position)
      (error 'bgzf-io-error :errno unix-ffi:*error-number*
             :text (format nil "failed to find position in ~a" bgzf)))
    position))

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
