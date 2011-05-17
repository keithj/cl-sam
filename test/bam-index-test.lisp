;;;
;;; Copyright (c) 2010-2011 Keith James. All rights reserved.
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

;; The suite will eventually contain fixtures for generated stock BAM
;; files used across multiple tests.
(deftestsuite bam-indexing-tests (cl-sam-tests)
  ())

(addtest (bam-indexing-tests) index-bam-file/1
  (with-tmp-bam-file (file)
    (let* ((num-refs 10)
           (ref-len 1000)
           (bam-file (generate-bam file num-refs ref-len))
           (index (index-bam-file bam-file)))
      (let ((n (length (bam-index-refs index))))
        (ensure (= num-refs n)
                :report "Expected ~d refs, but found ~d"
                :arguments (num-refs n))))))

(addtest (bam-indexing-tests) read/write-bam-index/1
  (with-tmp-bam-file (file)
    (let* ((num-refs 10)
           (ref-len 100000)
           (bam-file (generate-bam file num-refs ref-len))
           (index (index-bam-file bam-file))
           (index-file (merge-pathnames (make-pathname :type "bai") bam-file)))
      (with-open-file (out index-file :direction :output :element-type 'octet)
        (ensure (write-bam-index index out)))
      (with-open-file (in index-file :element-type 'octet)
        (ensure (equalp index (read-bam-index in))))
      (delete-file index-file))))

;; last chunk ends at 162913:0

(addtest (bam-indexing-tests) read-bam-range/1
  (with-tmp-bam-file (file)
    (let* ((num-refs 1)
           (ref-len 100000)
           (bam-file (generate-bam file num-refs ref-len :end ref-len))
           (index (index-bam-file bam-file)))
      (flet ((count-in-region (&rest regions)
               (with-bam (bam () bam-file :index index :regions regions)
                 (loop
                    while (has-more-p bam)
                    count (next bam)))))
        (mapcar (lambda (region expected)
                  (let ((count (count-in-region region)))
                    (ensure (= expected count)
                            :report "expected ~d, but found ~d"
                            :arguments (expected count))))
                (loop
                   for i from 0 to 90000 by 10000
                   collect (list 0 i (+ i 9999)))
                '(1971 2009 2009 2009 2009 2009 2009 2009 2009 1968))))))

;; (with-bam-index (index "/home/keith/index_test.bam.bai")
;;   (let ((regions (mapcar (lambda (x)
;;                            (apply #'region x)) '((1 1000 1010)
;;                                                  (1 1000000 2000000)
;;                                                  (1 1000 2000)
;;                                                  (2 1000 1000000)))))
;;     (with-bam (bam ()  "/home/keith/index_test.bam")
;;       (make-bam-region-input bam index regions))))

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
