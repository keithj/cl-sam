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

(in-package :sam-test)

(defun canonical-header (header)
  (mapcar (lambda (record)
            (cons (header-type record)
                  (sort (header-tags record) #'string<
                        :key (lambda (x)
                               (symbol-name (first x))))))
          (copy-tree header)))

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

(addtest (cl-sam-tests) merge-sam-headers/1
  (let ((header '((:hd (:vn . "1.3"))
                  (:sq (:sn . "AL096846") (:ln . 6490)
                   (:sp . "Schizosaccharomyces pombe"))
                  (:pg (:id . "bwa") (:vn . "0.4.6")))))
    (ensure (equalp header (merge-sam-headers header header)))
    (ensure (equalp (canonical-header header)
                    (canonical-header
                     (merge-sam-headers
                      '((:hd (:vn . "1.3"))
                        (:sq (:sn . "AL096846") (:ln . 6490)))
                      '((:hd (:vn . "1.3"))
                        (:sq (:sn . "AL096846") (:ln . 6490)
                         (:sp . "Schizosaccharomyces pombe"))
                        (:pg (:id . "bwa") (:vn . "0.4.6")))))))))

(addtest (cl-sam-tests) merge-sam-headers/2
  (let ((header1 '((:hd (:vn . "1.3"))
                   (:sq (:sn . "AL096846") (:ln . 6490))))
        (header2 '((:hd (:vn . "1.3"))
                   (:sq (:sn . "AL096846") (:ln . 9999)
                    (:sp . "Schizosaccharomyces pombe"))
                   (:pg (:id . "bwa") (:vn . "0.4.6")))))
    (ensure-condition invalid-operation-error
      (merge-sam-headers header1 header2))))

(addtest (cl-sam-tests) make-sam-header/1
  (let ((expected '((:hd (:vn . "1.3"))
                    (:sq (:sn . "AL096846") (:ln . 6490)
                     (:sp . "Schizosaccharomyces pombe"))))
        (result (make-sam-header "@HD	VN:1.3
@SQ	SN:AL096846	LN:6490	SP:Schizosaccharomyces pombe")))
    (ensure (equalp (canonical-header expected)
                    (canonical-header result)))))

(addtest (cl-sam-tests) make-sam-header/2 ; clashes
  (ensure-condition malformed-record-error
    (make-sam-header "@SQ	SN:AL096846	LN:6490	LN:9999"))) 

(addtest (cl-sam-tests) make-sam-header/3 ; duplicates
  (ensure (equalp '((:HD (:VN . "1.3") (:SO . :coordinate)))
                  (make-sam-header
                   "@HD	VN:1.3	SO:coordinate	SO:coordinate"))))

(addtest (cl-sam-tests) make-sam-header/4 ; Comment header
  (ensure (equalp '((:HD (:VN . "1.3")) (:CO . "Test comment."))
                  (make-sam-header
                   "@HD	VN:1.3
@CO	Test comment."))))

(addtest (cl-sam-tests) make-sam-header/5 ; group order
  (ensure-condition malformed-record-error
    (make-sam-header
     "@HD	VN:1.3	SO:coordinate	GO:none")))

(addtest (cl-sam-tests) subst-sort-order/1
  (ensure (equalp (canonical-header '((:HD (:VN . "1.3") (:SO . :coordinate))))
                  (canonical-header
                   (subst-sort-order (make-sam-header "@HD	VN:1.3")
                                     :coordinate)))))

(addtest (cl-sam-tests) subst-sort-order/2
  (ensure (equalp
           (canonical-header '((:HD (:VN . "1.3") (:SO . :coordinate))))
           (canonical-header (subst-sort-order
                              (make-sam-header "@HD	VN:1.3	SO:unsorted")
                              :coordinate)))))

;; You can still manipulate headers with group order
(addtest (cl-sam-tests) subst-group-order/1
  (ensure
   (equalp (canonical-header '((:HD (:VN . "1.3") (:GO . :none))))
           (canonical-header
            (subst-group-order (make-sam-header "@HD	VN:1.3")
                               :none)))))

(addtest (cl-sam-tests) subst-group-order/2
  (ensure (equalp
           (canonical-header '((:HD (:VN . "1.3") (:GO . :reference))))
           (canonical-header (subst-group-order
                              '((:HD (:VN . "1.3") (:GO . :none)))
                              :reference)))))
