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

(in-package :bgzf-ffi)

(define-foreign-library libbgzf
  (t (:default "libbgzf")))

(use-foreign-library libbgzf)

(defcfun ("bgzf_open" bgzf-open) :pointer
  (path :string)
  (mode :string))

(defcfun ("bgzf_close" bgzf-close) :int
  (fp :pointer))

;; (declaim (inline bgzf-read))
(defcfun ("bgzf_read" bgzf-read) :int
  (fp :pointer)
  (data :pointer)
  (length :int))

(defcfun ("bgzf_write" bgzf-write) :int
  (fp :pointer)
  (data :pointer)
  (length :int))

(defcfun ("bgzf_tell" bgzf-tell) :int64
  (fp :pointer))

(defcfun ("bgzf_seek" bgzf-seek) :int
  (fp :pointer)
  (pos unix-ffi:off-t)
  (where :int))
