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

(defun make-header-string (alist)
  "Returns a new SAM header string representing the header data in
ALIST."
 (with-output-to-string (s)
   (write-sam-header alist s)))

(defun write-sam-header (alist &optional (stream t))
  "Writes SAM header ALIST as a string to STREAM."
  (loop
     for rec in alist
     do (write-header-record rec stream)))

(defun write-header-record (header-record &optional (stream t))
  "Writes alist HEADER-RECORD to STREAM."
  (write-char #\@ stream)
  (write-string (symbol-name (first header-record)) stream)
  (write-char #\Tab stream)
  (dolist (x (intersperse (loop
                             for (tag . value) in (rest header-record)
                             collect (format nil "~a:~a" tag
                                             (typecase value
                                               (symbol (string-downcase
                                                        (symbol-name value)))
                                               (t value)))) #\Tab))
    (princ x stream))
  (terpri stream))
