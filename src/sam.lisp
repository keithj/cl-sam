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

(defparameter *mandatory-tags* (pairlis '(:hd :sq :rg :pg)
                                        '((:vn) (:sn :ln) (:id :sm) (:id)))
  "A mapping that describes the mandatory tags for each SAM header
  record type.")

(defmacro define-tag-parser (name (record-type var) tag-specs)
  "Defines a tag parsing function NAME that parses tag values for sam
header RECORD-TYPE."
  `(progn
    (defun ,name (str)
      (unless (and (< 3 (length str))
                   (char= #\: (char str 2)))
        (error 'malformed-field-error :field str
               :text (format nil "missing colon in tag:value")))
      (let ((tag (subseq str 0 2))
            (,var (subseq str 3)))
        (cond ,@(loop
                   for (tag-val sym form) in tag-specs
                   collect `((string= ,tag-val tag)
                             (list ,sym ,(if form
                                            `,form
                                          `,var))))
              (t
               (error 'malformed-field-error :field str
                      :text (format nil "invalid ~a tag ~a"
                                    ,record-type tag))))))))

(define-tag-parser parse-hd-tag ("HD" value)
  (("VN" :vn)
   ("SO" :so (cond ((string= "unsorted" value)
                    :unsorted)
                   ((string= "queryname" value)
                    :queryname)
                   ((string= "coordinate" value)
                    :coordinate)
                   (t
                    (error 'malformed-field-error :field value
                           :text "invalid SO value"))))
   ("GO" :go (cond ((string= "none" value)
                    :none)
                   ((string= "query" value)
                    :query)
                   ((string= "reference" value)
                    :reference)
                   (t
                    (error 'malformed-field-error :field value
                           :text "invalid GO value"))))))

(define-tag-parser parse-sq-tag ("SQ" value)
  (("SN" :sn) ("LN" :ln (parse-integer value))
   ("AS" :as) ("M5" :m5) ("UR" :ur) ("SP" :sp)))

(define-tag-parser parse-rg-tag ("RG" value)
  (("ID" :id) ("SM" :sm) ("LB" :lb) ("DS" :ds) ("PU" :pu) ("PI" :pi)
   ("CN" :cn) ("DT" :dt) ("PL" :pl)))

(define-tag-parser parse-pg-tag ("PG" value)
  (("ID" :id) ("VN" :vn) ("CL" :cl)))

(defun mandatory-tags (record-type)
  "Returns a list of the mandatory tags for SAM header RECORD-TYPE."
  (rest (assoc record-type *mandatory-tags*)))

(defun ensure-mandatory-tags (header-record)
  "Returns HEADER-RECORD if all mandatory tags are present, or raises
a {define-condition malformed-record-error} ."
  (let ((tag-type (first header-record))
        (tags (rest header-record)))
    (if (subsetp (mandatory-tags tag-type) (mapcar #'first tags))
        header-record
      (error 'malformed-record-error
             :record header-record
             :text (format nil "missing mandatory ~a tags"
                           (set-difference 
                            (mandatory-tags tag-type)
                            (mapcar #'first tags)))))))

(defun parse-header-record (str)
  "Parses a single SAM header record STR and returns an alist. May
raise a {define-condition malformed-record-error} or
{define-condition malformed-field-error} ."
  (flet ((tags (fn)
           (mapcar fn (rest (string-split str #\Tab)))))
    (cond ((starts-with-string-p str "@HD")
           (cons :hd (tags #'parse-hd-tag)))
          ((starts-with-string-p str "@SQ")
           (cons :sq (tags #'parse-sq-tag)))
          ((starts-with-string-p str "@RG")
           (cons :rg (tags #'parse-rg-tag)))
          ((starts-with-string-p str "@PG")
           (cons :pg (tags #'parse-pg-tag)))
          (t
           (error 'malformed-record-error
                  :record str
                  :text "invalid SAM header record")))))
