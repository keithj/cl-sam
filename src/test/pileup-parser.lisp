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

(in-package :sam-test)

;; This is used to parse SAMtools pileup output as part of testing.
(define-line-parser parse-pileup-record #\Tab
  ((:reference-name :type :string)
   (:reference-position :type :integer)
   (:reference-base :type :string)
   (:read-depth :type :integer)
   (:read-bases :type :string)
   (:base-qualities :type :string)))

(defun read-pileup-record (stream)
  (let ((line (find-line stream #'content-string-p)))
    (if (stringp line)
        (parse-pileup-record line)
      nil)))
