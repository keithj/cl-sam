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

(defun canonical-header (header)
  (mapcar (lambda (record)
            (cons (header-record-type record)
                  (sort (header-record-tags record) #'string<
                        :key (lambda (x)
                               (symbol-name (first x))))))
          (copy-tree header)))

(addtest (cl-sam-tests) make-header-record/1
  (ensure (equalp '(:sq (:sn . "AL096846") (:ln . 6490)
                    (:sp . "Schizosaccharomyces pombe"))
                  (make-header-record
                   "@SQ	SN:AL096846	LN:6490	SP:Schizosaccharomyces pombe"))))

(addtest (cl-sam-tests) make-header-record/2
  (ensure-condition malformed-record-error
    (make-header-record "@INVALID-TYPE	SN:AL096846")))

(addtest (cl-sam-tests) ensure-mandatory-tags/1
  (ensure (ensure-mandatory-tags '(:hd (:vn . "1.0"))))
  (ensure (ensure-mandatory-tags '(:sq (:sn . "foo") (:ln . 100))))
  (ensure (ensure-mandatory-tags '(:rg (:id . "foo") (:sm . "bar"))))
  (ensure (ensure-mandatory-tags '(:pg (:id . "foo")))))

(addtest (cl-sam-tests) ensure-mandatory-tags/2
  (ensure-condition malformed-record-error
    (ensure-mandatory-tags '(:hd nil)))
  (ensure-condition malformed-record-error
    (ensure-mandatory-tags '(:sq (:sn . "foo"))))
  (ensure-condition malformed-record-error
    (ensure-mandatory-tags '(:rg (:id . "foo"))))
  (ensure-condition malformed-record-error
    (ensure-mandatory-tags '(:pg nil))))

(addtest (cl-sam-tests) ensure-valid-tags/1
  (ensure (ensure-mandatory-tags '(:hd (:vn . "1.0")))))

(addtest (cl-sam-tests) ensure-mandatory-tags/2
  (ensure-condition malformed-record-error
    (ensure-valid-tags '(:hd (:invalid . "invalid")))))

(addtest (cl-sam-tests) merge-sam-headers/1
  (let ((header '((:hd (:vn . "1.0"))
                  (:sq (:sn . "AL096846") (:ln . 6490)
                   (:sp . "Schizosaccharomyces pombe"))
                  (:pg (:id . "bwa") (:vn . "0.4.6")))))
    (ensure (equalp header (merge-sam-headers header header)))
    (ensure (equalp (canonical-header header)
                    (canonical-header
                     (merge-sam-headers
                      '((:hd (:vn . "1.0"))
                        (:sq (:sn . "AL096846") (:ln . 6490)))
                      '((:hd (:vn . "1.0"))
                        (:sq (:sn . "AL096846") (:ln . 6490)
                         (:sp . "Schizosaccharomyces pombe"))
                        (:pg (:id . "bwa") (:vn . "0.4.6")))))))))

(addtest (cl-sam-tests) merge-sam-headers/2
  (let ((header1 '((:hd (:vn . "1.0"))
                   (:sq (:sn . "AL096846") (:ln . 6490))))
        (header2 '((:hd (:vn . "1.0"))
                   (:sq (:sn . "AL096846") (:ln . 9999)
                    (:sp . "Schizosaccharomyces pombe"))
                   (:pg (:id . "bwa") (:vn . "0.4.6")))))
    (ensure-condition invalid-operation-error
      (merge-sam-headers header1 header2))))

(addtest (cl-sam-tests) make-sam-header/1
  (let ((expected '((:hd (:vn . "1.0"))
                    (:sq (:sn . "AL096846") (:ln . 6490)
                     (:sp . "Schizosaccharomyces pombe"))))
        (result (make-sam-header "@HD	VN:1.0
@SQ	SN:AL096846	LN:6490	SP:Schizosaccharomyces pombe")))
    (ensure (equalp (canonical-header expected)
                    (canonical-header result)))))
