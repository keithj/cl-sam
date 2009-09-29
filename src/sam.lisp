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

(defparameter *valid-tags* (pairlis '(:hd :sq :rg :pg)
                                    '((:vn :so :go)
                                      (:sn :ln :as :m5 :ur :sp)
                                      (:id :sm :lb :ds :pu :pi :cn :dt :pl)
                                      (:id :vn :cl))))

(defmacro define-tag-parser (name (record-type var) tag-specs)
  "Defines a tag parsing function NAME that parses tag values for SAM
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
                             (cons ,sym ,(if form
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

(defun make-header-record (str)
  "Parses a single SAM header record STR and returns a list. May
raise a {define-condition malformed-record-error} or {define-condition
malformed-field-error} . SAM tags are converted to keyword symbols.

Given a header record of

;;; \"@SQ	SN:AL096846	LN:6490	SP:Schizosaccharomyces pombe\"

the returned list will be

;;;  (:SQ (:SN . \"AL096846\") (:LN . 6490)
;;;       (:SP . \"Schizosaccharomyces pombe\"))

thus the list's first element is a keyword describing the record type
and the rest of the list is itself an alist of record keys and values."
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

(defun header-record-type (record)
  "Returns a symbol indicating the record-type of HEADER-RECORD, being
one of :HD , :SQ , :RG or :PG ."
  (first record))

(defun header-record-tags (record)
  "Returns an alist of the tag-values of HEADER-RECORD."
  (rest record))

(defun mandatory-tags (record-type)
  "Returns a list of the mandatory tags for SAM header
RECORD-TYPE. Both RECORD-TYPE. and the returned tags are represented
as symbols."
  (rest (assoc record-type *mandatory-tags*)))

(defun ensure-mandatory-tags (record)
  "Returns HEADER-RECORD if all mandatory tags are present, or raises
a {define-condition malformed-record-error} ."
  (let ((tag-type (header-record-type record))
        (tags (header-record-tags record)))
    (unless (subsetp (mandatory-tags tag-type) (mapcar #'first tags))
      (error 'malformed-record-error
             :record record
             :text (let ((diff (set-difference (mandatory-tags tag-type)
                                               (mapcar #'first tags))))
                     (format nil "~r missing mandatory tag~:p ~a"
                             (length diff) diff))))
    record))

(defun valid-tags (record-type)
  "Returns a list of the valid tags for SAM header RECORD-TYPE. Both
RECORD-TYPE and the returned tags are represented as symbols."
  (rest (assoc record-type *valid-tags*)))

(defun ensure-valid-tags (record)
  "Checks list HEADER-RECORD for tag validity and returns
HEADER-RECORD or raises a {define-condition malformed-record-error} if
invalid tags are present."
  (let ((tag-type (header-record-type record))
        (tags (header-record-tags record)))
    (unless (subsetp (valid-tags tag-type) (mapcar #'first tags))
      (error 'malformed-record-error
             :record record
             :text (let ((diff (set-difference (mapcar #'first tags)
                                               (valid-tags tag-type))))
                     (format nil "~r invalid tag~:p ~a"
                             (length diff) diff))))
    record))

(defun merge-header-records (record1 record2)
  "Returns a new header record created by merging the tags of
header-records RECORD1 and RECORD2. Records may be safely merged if
they have the same record-type and do not contain any conflicting tag
values."
  (unless (eql (header-record-type record1) (header-record-type record2))
    (error 'invalid-operation-error
           :text (format nil
                         "invalid merge caused by different record-types in ~a"
                         (list record1 record2))))
  (let* ((merged (cons (header-record-type record1)
                       (remove-duplicates
                        (concatenate 'list
                                     (header-record-tags record1)
                                     (header-record-tags record2))
                        :test #'equalp)))
         (clashes (find-duplicate-tags merged)))
    (when clashes
      (error 'invalid-operation-error
             :text (format nil
                           "invalid merge caused by clashing tags ~a in ~a"
                           clashes (list record1 record2))))
    merged))

(defun make-sam-header (str)
  "Returns a list containing the data in header STR as Lisp
objects.

Given a header of

;;; \"@HD	VN:1.0
;;; @SQ	SN:AL096846	LN:6490	SP:Schizosaccharomyces pombe
;;; @PG	ID:bwa	VN:0.4.6	CL:aln -n 0.04 -o 1 -e -1 -i 5 -d 10 -k 2 -M 3 -O 11 -E 4\"

the returned list will be

;;; ((:HD (:VN . \"1.0\"))
;;;  (:SQ (:SN . \"AL096846\") (:LN . 6490)
;;;       (:SP . \"Schizosaccharomyces pombe\"))
;;;  (:PG (:ID . \"bwa\") (:VN . \"0.4.6\")
;;;       (:CL . \"aln -n 0.04 -o 1 -e -1 -i 5 -d 10 -k 2 -M 3 -O 11 -E 4\"))))

thus each list element is a list whose first element is a keyword
describing the record type. The rest of each list is itself an alist
of record keys and values."
  (with-input-from-string (s str)
    (loop
       for rec = (read-line s nil nil)
       while rec
       collect (make-header-record rec))))

(defun merge-sam-headers (&rest headers)
  "Returns a new SAM header that is the result of merging
HEADERS. Headers may be safely merged if none of their constituent
records contain conflicting tag values once merged."
  (multiple-value-bind (hd sq rg pg)
      (partition-by-type headers)
    (delete-if #'null
               (nconc (list (reduce #'merge-header-records
                                    (remove-duplicates hd :test #'equalp)))
                      (simplify-records sq :sn)
                      (simplify-records rg :id)
                      (simplify-records pg :id)))))

(defun partition-by-type (headers)
  "Collects all the header-records in HEADERS by record-type, sorts
them and returns 4 values that are lists of the collected :hd , :sq
, :rg and :pg header-records, respectively. Does not modify HEADERS."
  (flet ((sort-records (records tag)
           (stable-sort records #'string<
                        :key (lambda (x)
                               (assocdr tag (header-record-tags x))))))
    (let (hd sq rg pg)
      (dolist (header headers (values (nreverse hd)
                                      (sort-records sq :sn)
                                      (sort-records rg :id)
                                      (sort-records pg :id)))
        (dolist (record header)
          (ecase (header-record-type record)
            (:hd (push record hd))
            (:sq (push record sq))
            (:rg (push record rg))
            (:pg (push record pg))))))))

(defun simplify-records (records id-tag)
  "Returns a simplified list of header-records copied from
RECORDS. The simplifications are removal of perfect duplicates,
grouping of header-records by ID-TAG value and subsequent merging of
header-records that share that ID-TAG value."
  (let* ((unique (remove-duplicates records :test #'equalp))
         (groups (group-by-id unique id-tag))
         (simplified ()))
    (dolist (record unique (nreverse simplified))
      (let* ((id (assocdr id-tag (header-record-tags record)))
             (group (gethash id groups)))
        (if (endp (rest group))
            (push (first group) simplified)
          (push (reduce #'merge-header-records (nreverse group)) simplified))
        (remhash id groups)))))

(defun group-by-id (records id-tag)
  "Returns a hash-table of header-records taken from list RECORDS. The
hash-table keys are ID-TAG values taken from the records and the
hash-table values are lists of header-records that share that ID-TAG
value."
  (let ((groups (make-hash-table :test #'equalp)))
    (dolist (record records groups)
      (let* ((id (assocdr id-tag (header-record-tags record)))
             (group (gethash id groups)))
        (setf (gethash id groups) (if group
                                      (cons record group)
                                    (list record)))))))

(defun find-duplicate-tags (record)
  ;; hash-table size reflects the maximum number of possible tags in a
  ;; record according to the SAM spec
  (let ((tags (header-record-tags record))
        (duplicates (make-hash-table :size 10)))
    (loop
       for (tag . val) in tags
       for prev = (gethash tag duplicates)
       do (setf (gethash tag duplicates) (if prev
                                             (cons val prev)
                                           (list val)))
       finally (return (loop
                          for tag being the hash-keys of duplicates
                          using (hash-value vals)
                          unless (endp (rest vals))
                          collect (list tag (nreverse vals)))))))
