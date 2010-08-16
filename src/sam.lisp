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

(defparameter *sam-version* "0.1.2-draft (20090416)"
  "The SAM version written by cl-sam. This is a revision of the draft
  that does not support @CO header lines, which were introduced
  later.")

(defparameter *valid-header-types* '(:hd :sq :rg :pg)
  "A list of valid SAM header types.")
(defparameter *mandatory-header-tags*
  (pairlis *valid-header-types* '((:vn) (:sn :ln) (:id :sm) (:id)))
  "A mapping that describes the mandatory tags for each SAM header
  record type.")
(defparameter *valid-header-tags*
  (pairlis *valid-header-types* '((:vn :so :go)
                                  (:sn :ln :as :m5 :ur :sp)
                                  (:id :sm :lb :ds :pu :pi :cn :dt :pl)
                                  (:id :vn :cl)))
  "A mapping that describes the valid tags for each SAM header record
  type.")
(defparameter *valid-sort-orders* '(:unsorted :coordinate :queryname)
  "Valid values for SAM sort order tags.")
(defparameter *valid-group-orders* '(:none :query :reference)
  "Valid values for SAM group order tags.")

(defmacro define-header-tag-parser (name (header-type var) tag-specs)
  "Defines a tag parsing function NAME that parses tag values for SAM
header HEADER-TYPE."
  `(progn
    (defun ,name (str)
      (unless (and (< 3 (length str))
                   (char= #\: (char str 2)))
        (error 'malformed-field-error :field str
               :format-control "missing colon in tag:value"))
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
                      :format-control "invalid ~a tag ~a"
                      :format-arguments (list ,header-type tag))))))))

(define-header-tag-parser parse-hd-tag ("HD" value)
  (("VN" :vn)
   ("SO" :so (cond ((string= "unsorted" value)
                    :unsorted)
                   ((string= "queryname" value)
                    :queryname)
                   ((string= "coordinate" value)
                    :coordinate)
                   (t
                    (error 'malformed-field-error :field value
                           :format-control "invalid SO value"))))
   ("GO" :go (cond ((string= "none" value)
                    :none)
                   ((string= "query" value)
                    :query)
                   ((string= "reference" value)
                    :reference)
                   (t
                    (error 'malformed-field-error :field value
                           :format-control "invalid GO value"))))))

(define-header-tag-parser parse-sq-tag ("SQ" value)
  (("SN" :sn) ("LN" :ln (parse-integer value))
   ("AS" :as) ("M5" :m5) ("UR" :ur) ("SP" :sp)))

(define-header-tag-parser parse-rg-tag ("RG" value)
  (("ID" :id) ("SM" :sm) ("LB" :lb) ("DS" :ds) ("PU" :pu) ("PI" :pi)
   ("CN" :cn) ("DT" :dt) ("PL" :pl)))

(define-header-tag-parser parse-pg-tag ("PG" value)
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
    (let ((record
           (cond ((or (< (length str) 4) ; 3 type chars, 1 tab char
                      (char/= #\@ (char str 0))
                      (char/= #\Tab (char str 3)))
                  (error 'malformed-record-error :record str
                         :format-control "invalid SAM header record"))
                 ((starts-with-string-p str "@HD")
                  (cons :hd (tags #'parse-hd-tag)))
                 ((starts-with-string-p str "@SQ")
                  (cons :sq (tags #'parse-sq-tag)))
                 ((starts-with-string-p str "@RG")
                  (cons :rg (tags #'parse-rg-tag)))
                 ((starts-with-string-p str "@PG")
                  (cons :pg (tags #'parse-pg-tag)))
                 (t
                  (error 'malformed-record-error :record str
                         :format-control "invalid SAM header record type")))))
      ;; This is belt-and-braces because the parser already rejects
      ;; invalid tags
      (ensure-valid-header-tags (ensure-mandatory-header-tags record)))))

(defun header-type (record)
  "Returns a symbol indicating the header-type of HEADER-RECORD, being
one of :HD , :SQ , :RG or :PG ."
  (first record))

(defun header-tags (record)
  "Returns an alist of the tag-values of HEADER-RECORD."
  (rest record))

(defun mandatory-header-tags (header-type)
  "Returns a list of the mandatory tags for SAM header
HEADER-TYPE. Both HEADER-TYPE. and the returned tags are represented
as symbols."
  (rest (assoc header-type *mandatory-header-tags*)))

(defun valid-header-tags (header-type)
  "Returns a list of the valid tags for SAM header HEADER-TYPE. Both
HEADER-TYPE and the returned tags are represented as symbols."
  (rest (assoc header-type *valid-header-tags*)))

(defun ensure-mandatory-header-tags (record)
  "Returns HEADER-RECORD if all mandatory tags are present, or raises
a {define-condition malformed-record-error} ."
  (let* ((header-type (header-type record))
         (tags (header-tags record))
         (mandatory (mandatory-header-tags header-type))
         (tag-keys (mapcar #'first tags)))
    (unless (and mandatory (subsetp mandatory tag-keys))
      (let ((diff (set-difference mandatory tag-keys)))
      (error 'malformed-record-error :record record
             :format-control "~r missing mandatory tag~:p ~a"
             :format-arguments (list (length diff) diff))))
    record))

(defun ensure-valid-header-tags (record)
  "Checks list HEADER-RECORD for tag validity and returns
HEADER-RECORD or raises a {define-condition malformed-record-error} if
invalid tags are present."
  (let* ((header-type (header-type record))
         (tags (header-tags record))
         (valid (valid-header-tags header-type))
         (tag-keys (mapcar #'first tags)))
    (unless (subsetp tag-keys valid)
      (let ((diff (set-difference tag-keys valid)))
        (error 'malformed-record-error :record record
               :format-control "~r invalid tag~:p ~a"
               :format-arguments (list (length diff) diff))))
    record))

(defun merge-header-records (record1 record2)
  "Returns a new header record created by merging the tags of
header-records RECORD1 and RECORD2. Records may be safely merged if
they have the same header-type and do not contain any conflicting tag
values."
  (unless (eql (header-type record1) (header-type record2))
    (error 'invalid-operation-error
           :format-control (txt "invalid merge caused by different"
                                "record-types in ~a")
           :format-arguments (list record1 record2)))
  (let* ((merged (cons (header-type record1)
                       (remove-duplicates (concatenate 'list
                                                       (header-tags record1)
                                                       (header-tags record2))
                                          :test #'equalp)))
         (clashes (find-duplicate-header-tags merged)))
    (when clashes
      (error 'invalid-operation-error
             :format-control "invalid merge caused by clashing tags ~a in ~a"
             :format-arguments (list clashes (list record1 record2))))
    merged))

(defun make-sam-header (str)
  "Returns a list containing the data in header STR as Lisp
objects. If the header is NIL (there was no header) this function
returns NIL.

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
  (when str
    (with-input-from-string (s str)
      (loop
         for line = (read-line s nil nil)
         while line
         for (header-type . tags) = (make-header-record line)
         collect (let* ((tags (remove-duplicates tags :test #'equal))
                        (clashes (find-duplicate-header-tags tags))
                        (record (cons header-type tags)))
                   (when clashes
                     (error 'malformed-record-error :record record
                            :format-control "clashing tags ~a"
                            :format-arguments (list clashes)))
                   record)))))

(defun hd-record (&key (version *sam-version*) (sort-order :unsorted)
                  (group-order :none))
  (assert (stringp version) (version)
          "VERSION should be a string, but was ~a" version)
  (cons :hd (reverse (pairlis (valid-header-tags :hd)
                              (list version sort-order group-order)))))

(defun sq-record (seq-name seq-length &key assembly-identity seq-md5 seq-uri
                  seq-species)
  (assert (stringp seq-name)
          (seq-name)
          "SEQ-NAME should be a string, but was ~a" seq-name)
  (assert (and (integerp seq-length) (plusp seq-length))
          (seq-length)
          "SEQ-LENGTH should be a positive integer, but was ~a" seq-length)
  (cons :sq (remove-if #'null
                       (mapcar (lambda (key value)
                                 (when value
                                   (cons key value)))
                               (valid-header-tags :sq)
                               (list seq-name seq-length assembly-identity
                                     seq-md5 seq-uri seq-species)))))

(defun rg-record (identity sample &key library description
                  (platform-unit :lane) insert-size
                  sequencing-centre sequencing-date platform-tech)
    (assert (stringp identity)
            (identity)
            "IDENTITY should be a string, but was ~a" identity)
    (assert (stringp sample)
            (sample)
            "SAMPLE should be a string, but was ~a" sample)
    (assert (and (integerp insert-size) (plusp insert-size))
            (insert-size)
            "INSERT-SIZE should be a positive integer, but was ~a" insert-size)
    (cons :rg
          (remove-if #'null
                       (mapcar (lambda (key value)
                                 (when value
                                   (cons key value)))
                               (valid-header-tags :rg)
                               (list identity sample library description
                                     platform-unit insert-size
                                     sequencing-centre sequencing-date
                                     platform-tech)))))

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

(defun subst-sam-version (header &optional (version *sam-version*))
  "Returns a copy of HEADER with any header version tag initially
present substituted by string VERSION, defaulting to *SAM-VERSION* ."
  (let ((current (assoc :vn (assocdr :hd header))))
    ;; TODO -- warn if we are writing an older SAM version than we are
    ;; reading or if otherwise not backwardly compatible.
    (subst (cons :vn version) current header)))

(defun subst-sort-order (header order)
  "Returns a copy of HEADER with any sort order tag initially present
substituted by symbol ORDER, which must be one of the valid SAM sort
orders."
  (assert (member order *valid-sort-orders*) (order)
          "Invalid sort order ~a: expected one of ~a" order
          *valid-sort-orders*)
  (let* ((sorted (ensure-order header :sort))
         (current (assoc :so (assocdr :hd sorted))))
    (subst (cons :so order) current sorted)))

(defun subst-group-order (header order)
  "Returns a copy of HEADER with any group order tag initially present
substituted by symbol ORDER, which must be one of the valid SAM sort
orders."
  (assert (member order *valid-group-orders*) (order)
          "Invalid group order ~a: expected one of ~a" order
          *valid-group-orders*)
  (let* ((sorted (ensure-order header :group))
         (current (assoc :go (assocdr :hd sorted))))
    (subst (cons :go order) current sorted)))

(defun ensure-sam-version (header &optional (version *sam-version*))
  "Returns a copy of HEADER that is guaranteed to contain a version
  tag."
  (if (assoc :hd header)              ; Version is mandatory in header
      header
      (cons (list :hd (cons :vn version)) header)))

;;; SAM spec is silent on whether order within a record is
;;; important. For now we use acons and change the order because the
;;; spec doesn't forbid it.
(defun ensure-order (header domain)
  "Returns a copy of HEADER that is guaranteed to contain a sort tag
for DOMAIN, which must be one of :sort or :group ."
  (let ((hd (or (assoc :hd header)
                (list :hd (cons :vn *sam-version*)))) ; Default if :HD is absent
        (tag (ecase domain
               (:sort :so)
               (:group :go)))
        (value (ecase domain
               (:sort :unsorted)
               (:group :none))))
    (cond ((null header)
           (list hd))
          ((find tag (rest hd) :key #'car)
           header)
          (t
           (subst (cons :hd (acons tag value (rest hd)))
                  hd header :test #'equal)))))

(defun partition-by-type (headers)
  "Collects all the header-records in HEADERS by header-type, sorts
them and returns 4 values that are lists of the collected :hd , :sq
, :rg and :pg header-records, respectively. Does not modify HEADERS."
  (flet ((sort-records (records tag)
           (stable-sort records #'string<
                        :key (lambda (x)
                               (assocdr tag (header-tags x))))))
    (let (hd sq rg pg)
      (dolist (header headers (values (nreverse hd)
                                      (sort-records sq :sn)
                                      (sort-records rg :id)
                                      (sort-records pg :id)))
        (dolist (record header)
          (ecase (header-type record)
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
      (let* ((id (assocdr id-tag (header-tags record)))
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
      (let* ((id (assocdr id-tag (header-tags record)))
             (group (gethash id groups)))
        (setf (gethash id groups) (if group
                                      (cons record group)
                                      (list record)))))))

(defun find-duplicate-header-tags (record)
  "Returns a list of duplicate SAM header tags found in RECORD. Each
list element contains the tag, followed by a list of the conflicting
values."
  ;; hash-table size reflects the maximum number of possible tags in a
  ;; record according to the SAM spec
  (let ((tags (header-tags record))
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
