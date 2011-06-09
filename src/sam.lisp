;;;
;;; Copyright (c) 2009-2011 Keith James. All rights reserved.
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

(defparameter *sam-version* "1.4"
  "The SAM version written by cl-sam.")

(defparameter *valid-header-types* '(:hd :sq :rg :pg :co)
  "A list of valid SAM header types.")
(defparameter *mandatory-header-tags*
  (pairlis *valid-header-types* '((:vn) (:sn :ln) (:id) (:id) nil))
  "A mapping that describes the mandatory tags for each SAM header
  record type.")
(defparameter *valid-header-tags*
  (pairlis *valid-header-types* '((:vn :so :go)
                                  (:sn :ln :as :m5 :ur :sp)
                                  (:id :sm :lb :ds :pu :pi :cn :dt :pl)
                                  (:id :pn :vn :pp :cl)
                                  nil))
  "A mapping that describes the valid tags for each SAM header record
  type.")
(defparameter *valid-sort-orders* '(:unsorted :coordinate :queryname)
  "Valid values for SAM sort order tags.")
(defparameter *valid-group-orders* '(:none :query :reference)
  "Valid values for SAM group order tags. Group order is no longer a
valid tag in SAM version 1.3.")

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
        (cond ((find-if #'lower-case-p tag)
               (cons (intern tag :keyword) ,var))
              ,@(loop
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
  (("VN" :vn (if (valid-sam-version-p value)
                 value
                 (error 'malformed-field-error :field value
                        :format-control "invalid SAM version number")))
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
  (("SN" :sn (if (valid-reference-name-p value)
                 value
                 (error 'malformed-field-error :field value
                        :format-control "invalid SN value")))
   ("LN" :ln (parse-integer value))
   ("AS" :as) ("M5" :m5) ("SP" :sp) ("UR" :ur)))

(define-header-tag-parser parse-rg-tag ("RG" value)
  (("ID" :id) ("CN" :cn) ("DS" :ds)  ("DT" :dt) ("FO" :fo) ("KS" :ks)
   ("LB" :lb) ("PI" :pi) ("PG" :pg) ("PL" :pl) ("PU" :pu) ("SM" :sm)))

(define-header-tag-parser parse-pg-tag ("PG" value)
  (("ID" :id) ("CL" :cl) ("PN" :pn) ("PP" :pp) ("VN" :vn)))

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
                  (cons :hd (remove :go (tags #'parse-hd-tag)
                                    :key #'first))) ; parse, but discard :go
                 ((starts-with-string-p str "@SQ")
                  (cons :sq (tags #'parse-sq-tag)))
                 ((starts-with-string-p str "@RG")
                  (cons :rg (tags #'parse-rg-tag)))
                 ((starts-with-string-p str "@PG")
                  (cons :pg (tags #'parse-pg-tag)))
                 ((starts-with-string-p str "@CO")
                  (cons :co (second (string-split str #\Tab))))
                 (t
                  (error 'malformed-record-error :record str
                         :format-control "invalid SAM header record type")))))
      ;; This is belt-and-braces because the parser already rejects
      ;; invalid tags
      (ensure-valid-header-tags (ensure-mandatory-header-tags record)))))

(defun header-records (header header-type)
  "Returns a list of all records of HEADER-TYPE from HEADER."
  (remove-if (lambda (x)
               (not (eql x header-type))) header :key #'first))

(defun header-type (record)
  "Returns a symbol indicating the header-type of RECORD, being one
of :HD , :SQ , :RG or :PG ."
  (first record))

(defun header-tags (record)
  "Returns an alist of the tag-values of RECORD."
  (if (atom (rest record))
      nil
      (rest record)))

(defun header-value (record tag)
  "Returns the value associated with TAG in RECORD."
  (assocdr tag (header-tags record)))

(defun user-header-tag-p (tag)
  "Returns T if TAG is a SAM 1.4 user-defined header tag. User-defined
tags are recognised by containing lower case letters."
  (find-if #'lower-case-p (string tag)))

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
    (unless (subsetp mandatory tag-keys)
      (let ((diff (set-difference mandatory tag-keys)))
        (error 'malformed-record-error :record record
               :format-control "~r missing mandatory tag~:p ~a"
               :format-arguments (list (length diff) diff))))
    record))

(defun ensure-valid-header-tags (record)
  "Checks list HEADER-RECORD for tag validity and returns
HEADER-RECORD or raises a {define-condition malformed-record-error} if
invalid tags are present. Ignores any user tags (lower case tags)."
  (let* ((header-type (header-type record))
         (tags (remove-if #'user-header-tag-p (header-tags record)
                          :key #'first))
         (valid (valid-header-tags header-type))
         (tag-keys (mapcar #'first tags)))
    (unless (subsetp tag-keys valid)
      (let ((diff (set-difference tag-keys valid)))
        (error 'malformed-record-error :record record
               :format-control "~r invalid tag~:p ~a"
               :format-arguments (list (length diff) diff))))
    record))

(defun ensure-valid-programs (header)
  "Returns HEADER if its PG records have unique ID tag values and
valid PP tag values. A valid PP tag must point to the ID of one of the
other PG records in the header."
  (let* ((pg-records (header-records header :pg))
         (programs (group-by-tag pg-records :id)))
    (ensure-unique-tag-values
     (dolist (record pg-records header)
       (let ((pp (header-value record :pp)))
         (unless (or (null pp) (gethash pp programs))
           (error 'malformed-field-error :field pp :record record
                  :format-control "previous program ~s does not exist"
                  :format-arguments (list pp))))) :pg :id)))

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
      (ensure-valid-programs
       (ensure-unique-tag-values
        (loop
           for line = (read-line s nil nil)
           while line
           for (header-type . content) = (make-header-record line)
           collect (etypecase content
                     (string (cons header-type content)) ; :CO header
                     (list (let* ((tags (remove-duplicates content
                                                           :test #'equal))
                                  (clashes (find-duplicate-header-tags tags))
                                  (record (cons header-type tags)))
                             (when clashes
                               (error 'malformed-record-error :record record
                                      :format-control "clashing tags ~a"
                                      :format-arguments (list clashes)))
                             record)))) :sq :sn)))))

(defun name-sorted-p (header)
  "Returns T if parsed HEADER indicates that the file is sorted by
name, or NIL otherwise."
  (check-arguments (listp header) (header) "expected a parsed header")
  (eql :queryname (header-value (assoc :hd header) :so)))

(defun coordinate-sorted-p (header)
  "Returns T if parsed HEADER indicates that the file is sorted by
coordinate, or NIL otherwise."
  (check-arguments (listp header) (header) "expected a parsed header")
  (eql :coordinate (header-value (assoc :hd header) :so)))

(defun valid-sam-version-p (str)
  "Returns T if SAM version string STR matches /^[0-9]+.[0-9]$/, or
NIL otherwise."
  (let ((parts (string-split str #\.)))
    (and (= 2 (length parts))
         (every #'digit-char-p (first parts))
         (every #'digit-char-p (second parts)))))

(defun valid-reference-name-p (str)
  "Returns T if STR is a valid reference sequence name matching the
regex [!-)+-<>-~][!-~]* , or NIL otherwise."
  (labels ((name-char-p (c)             ; regex [!-)+-<>-~][!-~]*
             (< 32 (char-code c) 127))
           (start-char-p (c)
             (and (name-char-p c)
                  (not (member c '(#\* #\=))))))
    (let ((len (length str)))
      (and (plusp len)
           (start-char-p (char str 0))
           (loop
              for i from 1 below len
              always (name-char-p (char str i)))))))

(defun previous-programs (header identity)
  "Returns a list of PG ID values from HEADER that are previous
programs with respect to PG ID IDENTITY. The list is ordered with the
most recently used program first i.e. reverse chronological order."
  (let ((by-id (group-by-tag (header-records header :pg) :id)))
    (flet ((pp (id)
             (header-value (first (gethash id by-id)) :pp)))
      (loop
         for id = (pp identity) then (pp id)
         while id
         collect id))))

(defun last-programs (header)
  "Returns a list of the PG ID values from HEADER that are the last
programs to act on the data. i.e. these are the leaf programs in the
previous program tree."
  (let* ((prog-ids (mapcar (lambda (rec)
                             (header-value rec :id))
                           (header-records header :pg)))
         (all-paths (mapcar (lambda (id)
                              (previous-programs header id)) prog-ids)))
    (stable-sort
     (set-difference prog-ids (remove-duplicates
                               (apply #'concatenate 'list all-paths)
                               :test #'equal)
                     :test #'equal) #'string<)))

(defun header-record (record-type &rest args)
  "Returns a new header record of HEADER-TYPE. ARGS are tag values in
the same order as the tag returned by {defun valid-header-tags} ,
which is the same order as they are presented in the SAM spec."
  (cons record-type (remove-if #'null
                               (mapcar (lambda (key value)
                                         (when value
                                           (cons key value)))
                                       (valid-header-tags record-type) args))))

(defun hd-record (&key (version *sam-version*) (sort-order :unsorted))
  "Returns a new HD record."
  (assert (stringp version) (version)
          "VERSION should be a string, but was ~a" version)
  ;; (:vn :so :go)
  (cons :hd (reverse (pairlis (valid-header-tags :hd)
                              (list version sort-order "none")))))

(defun sq-record (seq-name seq-length &key assembly-identity seq-md5 seq-uri
                  seq-species)
  "Returns a new SQ record."
  (assert (stringp seq-name) (seq-name)
          "SEQ-NAME should be a string, but was ~a" seq-name)
  (assert (and (integerp seq-length) (plusp seq-length)) (seq-length)
          "SEQ-LENGTH should be a positive integer, but was ~a" seq-length)
  ;; (:sn :ln :as :m5 :ur :sp)
  (header-record :sq seq-name seq-length assembly-identity seq-md5 seq-uri
                 seq-species))

(defun rg-record (identity sample &key library description
                  (platform-unit :lane) (insert-size 0)
                  sequencing-centre sequencing-date platform-tech)
  "Returns a new RG record."
  (assert (stringp identity) (identity)
          "IDENTITY should be a string, but was ~a" identity)
  (assert (stringp sample) (sample)
          "SAMPLE should be a string, but was ~a" sample)
  (assert (and (integerp insert-size) (or (zerop insert-size)
                                          (plusp insert-size))) (insert-size)
          "INSERT-SIZE should be zero or a positive integer, but was ~a"
          insert-size)
  ;; (:id :sm :lb :ds :pu :pi :cn :dt :pl)
  (header-record :rg identity sample library description platform-unit
                 insert-size sequencing-centre sequencing-date platform-tech))

(defun pg-record (identity &key program-name program-version previous-program
                  command-line)
  "Returns a new PG record."
  (assert (stringp identity) (identity)
          "IDENTITY should be a string, but was ~a" identity)
  (assert (stringp program-version) (program-version)
          "PROGRAM-VERSION should be a string, but was ~a" program-version)
  ;; (:id :pn :vn :pp :cl)
  (header-record :pg identity program-name program-version previous-program
                 command-line))

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
orders. If there is no HD record in HEADER, one is added."
  (assert (member order *valid-sort-orders*) (order)
          "Invalid sort order ~a: expected one of ~a" order
          *valid-sort-orders*)
  (let* ((sorted (ensure-order header :sort))
         (current (assoc :so (assocdr :hd sorted))))
    (subst (cons :so order) current sorted)))

(defun subst-group-order (header order)
  "Returns a copy of HEADER with any group order tag initially present
substituted by symbol ORDER, which must be one of the valid SAM sort
orders. If there is no HD record in HEADER, one is added."
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

(defun ensure-unique-tag-values (header record-type tag)
  "Returns HEADER if all records of RECORD-TYPE have unique TAG
values, with respect to each other, or raises a {define-condition
malformed-field-error} ."
  (loop
     with values = (make-hash-table :test #'equal)
     for record in (header-records header record-type)
     do (let ((value (header-value record tag)))
          (cond ((null value)
                 nil)
                ((gethash value values)
                 (error 'malformed-field-error :field value :record record
                       :format-control "duplicate ~a ~s"
                       :format-arguments (list tag value)))
                (t
                 (setf (gethash value values) t))))
     finally (return header)))

;;; SAM spec is silent on whether order within a record is
;;; important. For now we use acons and change the order because the
;;; spec doesn't forbid it.
(defun ensure-order (header domain)
  "Returns a copy of HEADER that is guaranteed to contain a sort tag
for DOMAIN, which must be one of :sort or :group ."
  (let* ((header (ensure-sam-version header))
         (hd (assoc :hd header))
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
                        :key (lambda (record)
                               (header-value record tag)))))
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

(defun simplify-records (records header-tag)
  "Returns a simplified list of header-records copied from
RECORDS. The RECORDS must all be of the same type. The simplifications
are: removal of perfect duplicates, grouping of header-records by
HEADER-TAG value and subsequent merging of header-records that share
that HEADER-TAG value. RECORDS may be an empty list."
  (when records
    (check-arguments (= 1 (length (remove-duplicates records :key #'first)))
                     (records) "header records must be the same type")
    (check-arguments (member header-tag
                             (valid-header-tags (header-type (first records))))
                     (header-tag)
                     "~a is not valid tag for a ~a header" header-tag
                     (header-type (first records)))
    (let* ((unique (remove-duplicates records :test #'equalp))
           (groups (group-by-tag unique header-tag))
           (simplified ()))
      (dolist (record unique (nreverse simplified))
        (let* ((field (header-value record header-tag))
               (group (gethash field groups)))
          (if (endp (rest group))
              (push (first group) simplified)
              (push (reduce #'merge-header-records
                            (nreverse group)) simplified))
          (remhash field groups))))))

(defun group-by-tag (records header-tag)
  "Returns a hash-table of header-records taken from list RECORDS. The
hash-table keys are HEADER-TAG values taken from the records and the
hash-table values are lists of header-records that share that
HEADER-TAG value."
  (check-arguments (or (null records)
                       (= 1 (length (remove-duplicates records :key #'first))))
                   (records) "header records must be the same type")
  (let ((groups (make-hash-table :test #'equalp)))
    (dolist (record records groups)
      (let* ((field (header-value record header-tag))
             (group (gethash field groups)))
        (setf (gethash field groups) (if group
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

(defun find-duplicate-records (header header-type)
  "Returns a list of any duplicate records of HEADER-TYPE in HEADER."
  (check-arguments (listp header) (header) "expected a parsed header")
  (let ((records (header-records header header-type)))
    (set-difference records (remove-duplicates records :test #'equal))))

(defun add-pg-record (header new-record)
  "Returns a copy of HEADER with PG record NEW-RECORD added. If the ID
of NEW-RECORD clashes with existing IDs, all IDs are remapped to new,
generated sequential integer IDs, starting at 0. PP links are also
updated. NEW-RECORD must have its PP field set by the caller."
  (multiple-value-bind (hd sq rg pg)
      (partition-by-type (list header))
    (declare (ignore sq))
    (nconc hd (header-records header :sq) ; use SQ records in original order
           rg (update-pg-records pg new-record))))

(defun update-pg-records (current-records new-record)
  "Returns a copy of CURRENT-RECORDS with NEW-RECORD added."
  (let* ((current-records current-records)
         (new-id (header-value new-record :id))
         (new-pp (header-value new-record :pp))
         (current-ids (mapcar (lambda (rec)
                                (header-value rec :id)) current-records))
         (num-ids (iota (length current-records)))
         (current-pps (mapcar (lambda (rec)
                                (header-value rec :pp)) current-records)))
    (check-arguments (or (null new-pp) (find new-pp current-ids
                                             :test #'string=))
                     (new-record)
                     "previous program ~s not found" new-pp)
    (labels ((map-id (id)
               "Map old ID to new, numeric identifier"
               (when id                 ; null maps to null
                 (elt num-ids (position id current-ids :test #'equal))))
             (map-kv (k v rec)
               "Substitute conses in REC bearing IDs"
               (if v
                   (subst (cons k (map-id v)) (cons k v) rec :test #'equal)
                   rec)))
      (nreverse
       (if (member new-id current-ids :test #'string=) ; ID clash
           (let ((mod-records (mapcar (lambda (id pp record)
                                        (map-kv :pp pp (map-kv :id id record)))
                                      current-ids current-pps current-records))
                 (num-id (format nil "~d" (length current-records)))
                 (current-pp (header-value new-record :pp)))
             (cons (subst (cons :pp (map-id current-pp)) ; substitute PP
                          (cons :pp current-pp)
                          (subst (cons :id num-id) ; substitute ID
                                 (cons :id new-id) new-record :test #'equal)
                          :test #'equal)
                   (reverse mod-records)))
           (cons new-record (reverse current-records)))))))
