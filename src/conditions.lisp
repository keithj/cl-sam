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

(define-condition bgzf-io-error (io-error simple-text-condition)
  ((errno :initform nil
          :initarg :errno
          :reader errno-of
          :documentation "The C error number."))
  (:report (lambda (condition stream)
             (format stream "BGZF error~@[ ~d~]~@[: ~a~]"
                     (errno-of condition) (message-of condition))))
  (:documentation "A condition raised when an error occurs reading
  from or writing to a BGZF file."))

(define-condition bam-error (error formatted-condition)
  ()
  (:documentation "The parent type of all BAM  error conditions."))

(define-condition bam-sort-error (bam-error)
  ((prev-position :initform nil
                  :initarg :prev-position
                  :reader prev-position-of)
   (position :initform nil
             :initarg :position
             :reader position-of)
   (prev-reference :initform nil
                   :initarg :prev-reference
                   :reader prev-reference-of)
   (reference :initform nil
              :initarg :reference
              :reader reference-of))
  (:report (lambda (condition stream)
             (format stream "BAM sort error ~a ~a~@[: ~a~]"
                     (list :previous (prev-reference-of condition)
                           (prev-position-of condition))
                     (list :current (reference-of condition)
                           (position-of condition))
                     (message-of condition)))))
