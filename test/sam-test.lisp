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

(in-package :sam-test)

(defun canonical-header (header)
  (mapcar (lambda (record)
            (cons (header-type record)
                  (sort (header-tags record) #'string<
                        :key (lambda (x)
                               (symbol-name (first x))))))
          (copy-tree header)))

(defun header-equal (header1 header2)
  (equal (canonical-header header1) (canonical-header header2)))

(addtest (cl-sam-tests) valid-sam-version-p/1
  (ensure (valid-sam-version-p "1.0"))
  (ensure (valid-sam-version-p "1.3"))
  (ensure (not (valid-sam-version-p "1.0.0"))))

(addtest (cl-sam-tests) valid-reference-name-p/1
  (ensure (valid-reference-name-p "x"))
  (ensure (valid-reference-name-p "x*"))
  (ensure (valid-reference-name-p "x="))
  (ensure (not (valid-reference-name-p "")))
  (ensure (not (valid-reference-name-p " ")))
  (ensure (not (valid-reference-name-p "*")))
  (ensure (not (valid-reference-name-p "="))))

(addtest (cl-sam-tests) make-header-record/1
  (ensure (equalp '(:sq (:sn . "AL096846") (:ln . 6490)
                    (:sp . "Schizosaccharomyces pombe"))
                  (make-header-record
                   "@SQ	SN:AL096846	LN:6490	SP:Schizosaccharomyces pombe"))))

(addtest (cl-sam-tests) make-header-record/2
  (ensure-condition malformed-record-error
    (make-header-record "@INVALID-TYPE	SN:AL096846")))

(addtest (cl-sam-tests) make-header-record/3
  (ensure (equalp '(:co . "Test comment.")
                  (make-header-record "@CO	Test comment."))))

(addtest (cl-sam-tests) make-header-record/4
  (ensure-condition malformed-field-error
    (make-header-record "@HD	VN:x"))
  (ensure-condition malformed-field-error
    (make-header-record "@HD	VN:10"))
  (ensure-condition malformed-field-error
    (make-header-record "@HD	VN:0..0")))

(addtest (cl-sam-tests) make-header-record/5
  (ensure-condition malformed-field-error
    (make-header-record "@SQ	SN:"))
  (ensure-condition malformed-field-error
    (make-header-record "@SQ	SN:="))
  (ensure-condition malformed-field-error
    (make-header-record "@SQ	SN:*")))

(addtest (cl-sam-tests) ensure-mandatory-header-tags/1
  (ensure (ensure-mandatory-header-tags '(:hd (:vn . "1.3"))))
  (ensure (ensure-mandatory-header-tags '(:sq (:sn . "foo") (:ln . 100))))
  (ensure (ensure-mandatory-header-tags '(:rg (:id . "foo") (:sm . "bar"))))
  ;; :sm tag is no longer mandatory in SAM 1.3
  (ensure (ensure-mandatory-header-tags '(:rg (:id . "foo"))))
  (ensure (ensure-mandatory-header-tags '(:pg (:id . "foo"))))
  ;; CO has no tags
  (ensure (ensure-mandatory-header-tags '(:co . "foo"))))

(addtest (cl-sam-tests) ensure-mandatory-header-tags/2
  (ensure-condition malformed-record-error
    (ensure-mandatory-header-tags '(:hd nil)))
  (ensure-condition malformed-record-error
    (ensure-mandatory-header-tags '(:sq (:sn . "foo"))))
  (ensure-condition malformed-record-error
    (ensure-mandatory-header-tags '(:pg nil))))

(addtest (cl-sam-tests) ensure-valid-header-tags/1
  (ensure (ensure-mandatory-header-tags '(:hd (:vn . "1.3")))))

(addtest (cl-sam-tests) ensure-valid-header-tags/2
  (ensure-condition malformed-record-error
    (ensure-valid-header-tags '(:hd (:invalid . "invalid")))))

(addtest (cl-sam-tests) ensure-valid-header-tags/3
  ;; CO has no tags
  (ensure-valid-header-tags '(:co . "foo")))

(addtest (cl-sam-tests) previous-programs/1
  (let ((header (make-sam-header "@HD	VN:1.3
@SQ	SN:AL096846	LN:6490	SP:Schizosaccharomyces pombe
@PG	ID:1	PN:program a	VN:version 1.0	CL:run_program_a -a 1 -b 2
@PG	ID:2	PN:program a	VN:version 1.0	CL:run_program_a -x 1 -y 2
@PG	ID:3	PN:program b	VN:version 1.0	CL:run_program_b -a 1 -b 2	PP:1
@PG	ID:4	PN:program b	VN:version 1.0	CL:run_program_b -x 1 -y 2	PP:2
@PG	ID:5	PN:program c	VN:version 1.0	CL:run_program_c -a 1 -b 2	PP:2
@PG	ID:6	PN:program d	VN:version 1.0	CL:run_program_d -a 1 -b 2	PP:5
@PG	ID:7	PN:program e	VN:version 1.0	CL:run_program_e -a 1 -b 2	PP:6
@PG	ID:8	PN:program f	VN:version 1.0	CL:run_program_f -a 1 -b 2")))
    (ensure (null (previous-programs header "8")))
    (ensure (equalp '("6" "5" "2") (previous-programs header "7")))
    (ensure (equalp '("5" "2") (previous-programs header "6")))
    (ensure (equalp '("2") (previous-programs header "5")))
    (ensure (equalp '("2") (previous-programs header "4")))
    (ensure (equalp '("1") (previous-programs header "3")))
    (ensure (null (previous-programs header "2")))
    (ensure (null (previous-programs header "1")))))

(addtest (cl-sam-tests) last-programs/1
  (ensure (null (last-programs (make-sam-header "@HD	VN:1.3"))))
  (ensure (equal "1" (last-programs (make-sam-header "@HD	VN:1.3
@PG	ID:1	PN:program a	VN:version 1.0	CL:run_program_a -a 1 -b 2")))))

(addtest (cl-sam-tests) last-programs/1
  (let ((header (make-sam-header "@HD	VN:1.3
@SQ	SN:AL096846	LN:6490	SP:Schizosaccharomyces pombe
@PG	ID:1	PN:program a	VN:version 1.0	CL:run_program_a -a 1 -b 2
@PG	ID:2	PN:program a	VN:version 1.0	CL:run_program_a -x 1 -y 2
@PG	ID:3	PN:program b	VN:version 1.0	CL:run_program_b -a 1 -b 2	PP:1
@PG	ID:4	PN:program b	VN:version 1.0	CL:run_program_b -x 1 -y 2	PP:2
@PG	ID:5	PN:program c	VN:version 1.0	CL:run_program_c -a 1 -b 2	PP:2
@PG	ID:6	PN:program d	VN:version 1.0	CL:run_program_d -a 1 -b 2	PP:5
@PG	ID:7	PN:program e	VN:version 1.0	CL:run_program_e -a 1 -b 2	PP:6
@PG	ID:8	PN:program f	VN:version 1.0	CL:run_program_f -a 1 -b 2")))
    (ensure (equalp '("3" "4" "7" "8") (last-programs header)))))

(addtest (cl-sam-tests) add-pg-record/1
  (let ((header (make-sam-header "@HD	VN:1.3
@SQ	SN:AL096846	LN:6490	SP:Schizosaccharomyces pombe
@PG	ID:x	PN:program a	VN:version 1.0	CL:run_program_a -a 1 -b 2
@PG	ID:y	PN:program b	VN:version 1.0	CL:run_program_b -a 1 -b 2	PP:x")))
    ;; No PP, no ID clash
    (ensure (header-equal
             '((:HD (:VN . "1.3"))
               (:SQ (:SN . "AL096846") (:LN . 6490)
                (:SP . "Schizosaccharomyces pombe"))
               (:PG (:ID . "x") (:PN . "program a") (:VN . "version 1.0")
                (:CL . "run_program_a -a 1 -b 2"))
               (:PG (:ID . "y") (:PN . "program b") (:VN . "version 1.0")
                (:CL . "run_program_b -a 1 -b 2") (:PP . "x"))
               (:PG (:ID . "z") (:PN . "program c") (:VN . "2.0")
                (:CL . "run_program_c -a 1 -b 2")))
             (add-pg-record
              header
              (pg-record "z"
                         :program-name "program c"
                         :program-version "2.0"
                         :command-line "run_program_c -a 1 -b 2")))
            :report "No PP, no ID clash failed")
    ;; No PP, ID clash
    (ensure (header-equal
             '((:HD (:VN . "1.3"))
               (:SQ (:SN . "AL096846") (:LN . 6490)
                (:SP . "Schizosaccharomyces pombe"))
               (:PG (:ID . 0) (:PN . "program a") (:VN . "version 1.0")
                (:CL . "run_program_a -a 1 -b 2"))
               (:PG (:ID . 1) (:PN . "program b") (:VN . "version 1.0")
                (:CL . "run_program_b -a 1 -b 2") (:PP . 0))
               (:PG (:ID . "2")  (:PN . "program c") (:VN . "2.0")
                (:CL . "run_program_c -a 1 -b 2")))
             (add-pg-record
              header
              (pg-record "x"
                         :program-name "program c"
                         :program-version "2.0"
                         :command-line "run_program_c -a 1 -b 2")))
            :report "No PP, ID clash failed")
    ;; PP, no ID clash
    (ensure (header-equal
             '((:HD (:VN . "1.3"))
               (:SQ (:SN . "AL096846") (:LN . 6490)
                (:SP . "Schizosaccharomyces pombe"))
               (:PG (:ID . "x") (:PN . "program a") (:VN . "version 1.0")
                (:CL . "run_program_a -a 1 -b 2"))
               (:PG (:ID . "y") (:PN . "program b") (:VN . "version 1.0")
                (:CL . "run_program_b -a 1 -b 2") (:PP . "x"))
               (:PG (:ID . "z") (:PN . "program c") (:VN . "2.0") (:PP . "y")
                (:CL . "run_program_c -a 1 -b 2")))
             (add-pg-record
              header
              (pg-record "z"
                         :program-name "program c"
                         :program-version "2.0"
                         :command-line "run_program_c -a 1 -b 2"
                         :previous-program "y")))
            :report "PP, no ID clash failed")
    ;; PP, ID clash
    (ensure (header-equal
             '((:HD (:VN . "1.3"))
               (:SQ (:SN . "AL096846") (:LN . 6490)
                (:SP . "Schizosaccharomyces pombe"))
               (:PG (:ID . 0) (:PN . "program a") (:VN . "version 1.0")
                (:CL . "run_program_a -a 1 -b 2"))
               (:PG (:ID . 1) (:PN . "program b") (:VN . "version 1.0")
                (:CL . "run_program_b -a 1 -b 2") (:PP . 0))
               (:PG (:ID . "2") (:PN . "program c") (:VN . "2.0")
                (:PP . 1) (:CL . "run_program_c -a 1 -b 2")))
             (add-pg-record header
                            (pg-record "x"
                                       :program-name "program c"
                                       :program-version "2.0"
                                       :command-line "run_program_c -a 1 -b 2"
                                       :previous-program "y")))
            :report "PP, ID clash failed")))

(addtest (cl-sam-tests) add-pg-record/2
  (let ((header (make-sam-header "@HD	VN:1.3
@SQ	SN:AL096846	LN:6490	SP:Schizosaccharomyces pombe
@PG	ID:x	PN:program a	VN:version 1.0	CL:run_program_a -a 1 -b 2")))
    (ensure-condition invalid-argument-error
      (add-pg-record header
                     (pg-record "x"
                                :program-name "program c"
                                :program-version "2.0"
                                :command-line "run_program_c -a 1 -b 2"
                                :previous-program "!! no such PP !!")))))

(addtest (cl-sam-tests) user-header-tag/1
  (ensure (header-equal '((:HD (:VN . "1.3"))
                          (:SQ (:SN . "AL096846") (:LN . 6490)
                           (:SP . "Schizosaccharomyces pombe"))
                          (:RG (:ID . "1") (:SM . "a sample")
                           (:|Fo| . "foo")))
                        (make-sam-header "@HD	VN:1.3
@SQ	SN:AL096846	LN:6490	SP:Schizosaccharomyces pombe
@RG	ID:1	SM:a sample	Fo:foo")))
  (ensure (header-equal '((:HD (:VN . "1.3"))
                          (:SQ (:SN . "AL096846") (:LN . 6490)
                           (:SP . "Schizosaccharomyces pombe"))
                          (:RG (:ID . "1") (:SM . "a sample")
                           (:|fO| . "foo")))
                        (make-sam-header "@HD	VN:1.3
@SQ	SN:AL096846	LN:6490	SP:Schizosaccharomyces pombe
@RG	ID:1	SM:a sample	fO:foo"))))

(addtest (cl-sam-tests) user-header-tag/2
  (ensure-condition malformed-field-error
    (make-header-record "@RG	Fo:")))
