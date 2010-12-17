;;;
;;; Copyright (C) 2010 Keith James. All rights reserved.
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

(addtest (cl-sam-tests) read-bam-index/1
  (let* ((bam-file (generate-bam-file (tmp-pathname
                                       :tmpdir (merge-pathnames "data")
                                       :type "bam")
                                      :ref-length 100000 :read-length 50
                                      :insert-length 150 :step-size 100))
         (index (index-bam-file bam-file)))
    (ensure (probe-file bam-file))
    (let ((num-refs (length (sam::bam-index-refs index)))
          (expected 10))
      (ensure (= expected num-refs)
              :report "Expected ~d refs, but found ~d"
              :arguments (expected num-refs)))
    (delete-file bam-file)))
