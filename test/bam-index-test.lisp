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
  (let* ((num-refs 10)
         (ref-len 100000)
         (aln-sources (loop
                         for ref-id from 0 below num-refs
                         collect (alignment-generator ref-id "read_group_0"
                                                      :read-length 50
                                                      :insert-length 150
                                                      :step-size 100)))
         (bam-file (apply #'generate-bam-file (tmp-pathname
                                               :tmpdir (merge-pathnames "data")
                                               :type "bam") num-refs ref-len
                                               aln-sources))
         (index (index-bam-file bam-file)))
    (ensure (probe-file bam-file))
    (let ((expected 10))
      (ensure (= expected num-refs)
              :report "Expected ~d refs, but found ~d"
              :arguments (expected num-refs)))
    (delete-file bam-file)))

;; (let ((bam-file (tmp-pathname :tmpdir "/home/keith/" :type "bam"))
;;       (g1 (alignment-generator 0 "generated_group" :name-suffix "1"
;;                                :insert-length 100
;;                                :end (- 100000 (* 2 50) 100)))
;;       (g2 (alignment-generator 0 "generated_group"
;;                                :insert-length 1000 :name-suffix "2"
;;                                :end (- 100000 (* 2 50) 1000)))
;;       (g3 (alignment-generator 0 "generated_group"
;;                                :insert-length 10000 :name-suffix "3"
;;                                :end (- 100000 (* 2 50) 10000))))
;;   (generate-reference-file (merge-pathnames (make-pathname :type "fasta")
;;                                             bam-file)
;;                            "ref_0" 100000 (lambda ()
;;                                            #\a))
;;   (generate-bam-file bam-file 1 100000 g1 g2 g3))
