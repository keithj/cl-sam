;;;
;;; Copyright (c) 2009-2013 Keith James. All rights reserved.
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

(deftestsuite cl-sam-tests ()
  ())

(defmacro with-tmp-bam-file ((file) &body body)
  `(handler-bind ((error #'leave-tmp-pathname))
     (with-tmp-pathname (,file :tmpdir (merge-pathnames "data") :type "bam")
       ,@body)))

(defun make-aln-generators (num-refs read-group &rest args)
  (loop
     for ref-id from 0 below num-refs
     collect (apply #'alignment-generator ref-id read-group args)))

(defun generate-bam (filespec num-refs ref-len &rest args)
  (let* ((read-group "read_group_0")
         (gen (apply #'make-aln-generators num-refs read-group args)))
    (apply #'generate-bam-file filespec num-refs ref-len (list read-group)
           gen)))

(defun test-binary-files-eql (filespec1 filespec2)
  (with-open-file (s1 filespec1 :element-type 'octet)
    (with-open-file (s2 filespec2 :element-type 'octet)
      (ensure (loop
                 for byte1 = (read-byte s1 nil nil)
                 for byte2 = (read-byte s2 nil nil)
                 while (and byte1 byte2)
                 always (= byte1 byte2)
                 finally (cond ((null byte1)
                                (file-position s2 (1- (file-position s2))))
                               ((null byte2)
                                (file-position s1 (1- (file-position s1))))))
              :report "files were not byte identical"))))

(defun test-sam-equal (filespec1 filespec2)
  (with-bam (in1 (header1 num-refs1 ref-meta1) filespec1)
    (with-bam (in2 (header2 num-refs2 ref-meta2) filespec2)
      (loop
         with refs1 = (make-reference-table ref-meta1)
         with refs2 = (make-reference-table ref-meta2)
         while (has-more-p in1)
         for x = (next in1)
         for y = (next in2)
         do (ensure (equalp x y)
                    :report "expected ~a but found ~a"
                    :arguments ((with-output-to-string (s)
                                  (write-sam-alignment x refs1 s))
                                (with-output-to-string (s)
                                  (write-sam-alignment y refs2 s))))))))
