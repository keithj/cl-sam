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

;; FIXME -- add tests that seek into a different block
(addtest (cl-sam-tests) bgzf-seek/tell/1
  (with-bgzf (bgzf (merge-pathnames "data/c1215-name.bam") :direction :input)
    (ensure (bgzf-open-p bgzf) :report "expected an open handle")
    (dotimes (n 1)
      (let* ((x (random 4096))
             (y (random 4096))
             (bytes (sam::read-bytes bgzf x))
             (pos (bgzf-tell bgzf)))
        (declare (ignore bytes))
        ;; Read further
        (let ((more-bytes (sam::read-bytes bgzf y)))
          ;; Seek back to where we were
          (ensure (bgzf-seek bgzf pos))
          ;; Read those bytes again
          (ensure (= pos (bgzf-tell bgzf)))
          (ensure (equalp more-bytes (sam::read-bytes bgzf y))))))))

(addtest (cl-sam-tests) bam-open/close/1
  (let* ((filespec (merge-pathnames "data/c1215-name.bam"))
         (bgzf (bgzf-open filespec)))
    (ensure-condition file-error
      (bgzf-open "this/file/does/not/exist"))
    (ensure (bgzf-open-p bgzf)
            :report "expected an open handle")
    (ensure (bgzf-close bgzf)
            :report "failed to close handle")))

(addtest (cl-sam-tests) bam-open/close/2
  (let* ((filespec (merge-pathnames "data/c1215-name.bam"))
         (bgzfspec (bgzf-open filespec))
         (bgzf (bgzf-open bgzfspec)))
    (ensure (bgzf-open-p bgzf)
            :report "expected an open handle")
    (ensure (bgzf-close bgzf)
            :report "failed to close handle")))

(addtest (cl-sam-tests) bam-open-close/3
  (let* ((filespec (merge-pathnames "data/c1215-name.bam"))
         (bgzf1 (bgzf-open filespec))
         (halfway 10000))
    (read-bam-meta bgzf1)
    (loop
       repeat halfway
       do (read-alignment bgzf1))
    (let* ((bgzf2 (bgzf-open bgzf1)) ; reopen at arbitrary position
           (n (loop
                 for aln = (read-alignment bgzf2)
                 while aln
                 count aln)))
      (ensure (= halfway n)
              :report "Expected ~d, but found ~d"
              :arguments (halfway n)))))

(addtest (cl-sam-tests) bgzf-stream-read-byte/1
  (with-open-file (raw (merge-pathnames "data/c1215-name.dat")
                       :element-type 'octet)
    (let ((stream (bgzip-open (merge-pathnames "data/c1215-name.bam")
                              :direction :input)))
      (ensure (loop
                 for byte = (stream-read-byte stream)
                 until (eql :eof byte)
                 always (= (read-byte raw) byte)))
      (ensure (eql :eof (read-byte raw nil :eof))
              :report "expected raw data stream to be at EOF.")
      (ensure (close stream) :report "failed to close stream"))))

(addtest (cl-sam-tests) bgzf-stream-read-sequence/1
  (with-open-file (raw (merge-pathnames "data/c1215-name.dat")
                       :element-type 'octet)
    (let ((stream (bgzip-open (merge-pathnames "data/c1215-name.bam")
                              :direction :input))
          (buffer1 (make-array 100 :element-type 'octet :initial-element 0))
          (buffer2 (make-array 100 :element-type 'octet :initial-element 0)))
      (ensure (loop
                 for i from 0 below (file-length raw) by 100
                 for j = (stream-read-sequence stream buffer1 0 100)
                 for k = (read-sequence buffer2 raw :start 0 :end 100)
                 always (equalp buffer1 buffer2)))
      (ensure (close stream) :report "failed to close stream"))))

(addtest (cl-sam-tests) bam-parse/1
  (let* ((filespec (merge-pathnames "data/c1215-name.bam"))
         (bgzf (bgzf-open filespec)))
    (unwind-protect
         (progn
           (ensure (read-bam-magic bgzf))
           (ensure-null (read-bam-header bgzf)
                        :report "expected a null header") ; no header
           (ensure (= 1 (read-num-references bgzf)))
           (multiple-value-bind (ref len)
               (read-reference-meta bgzf)
             (ensure (string= "AL096846" ref)
                     :report "expected reference name \"AL096846\" but found ~s"
                     :arguments (ref))
             (ensure (= 6490 len)
                     :report "expected reference length 6490 but found ~d"
                     :arguments (len)))
           (loop
              for alignment = (read-alignment bgzf)
              while alignment
              do (progn
                   (ensure (stringp (read-name alignment)))
                   (let ((core (alignment-core alignment)))
                     (ensure (listp core))
                     (ensure (= 11 (length core))))
                   (let ((seq (seq-string alignment))
                         (qual (quality-string alignment)))
                     (ensure (stringp seq))
                     (ensure (stringp qual))
                     (ensure (= (length seq) (length qual))))
                   (ensure (listp (alignment-tag-values alignment))))
              count alignment into n
              finally (ensure (= 20000 n)
                              :report "expected ~d alignments, but read ~d"
                              :arguments (20000 n))))
      (bgzf-close bgzf))))

(addtest (cl-sam-tests) bam-parse/2
  (let ((expected (with-open-file (s (merge-pathnames
                                      "data/c1215_fixmate_10reads.sexp"))
                    (read s)))
        (found (with-bgzf (bgzf (namestring (merge-pathnames
                                             "data/c1215_fixmate.bam"))
                                :direction :input)
                 (multiple-value-bind (header num-refs ref-meta)
                     (read-bam-meta bgzf)
                   (declare (ignore header num-refs ref-meta))
                   (loop
                      repeat 10
                      for alignment = (read-alignment bgzf)
                      collect (list (alignment-core-alist alignment)
                                    (alignment-flag-alist alignment)
                                    (alignment-tag-values alignment)
                                    (alignment-cigar alignment)
                                    (read-name alignment)
                                    (seq-string alignment)
                                    (quality-string alignment)
                                    (alignment-tag-values alignment)))))))
    (dotimes (i 10)
      (let ((expected (elt expected i))
            (found (elt found i)))
        (ensure (equalp expected found)
                :report "alignment ~d: expected ~a but found ~a"
                :arguments (i expected found))))))

(addtest (cl-sam-tests) sequenced-pair-p/1
  (ensure (sequenced-pair-p #x0001)))

(addtest (cl-sam-tests) mapped-proper-pair-p/1
  (ensure (mapped-proper-pair-p #x0002)))

(addtest (cl-sam-tests) query-unmapped-p/1
  (ensure (query-unmapped-p #x0004)))

(addtest (cl-sam-tests) mate-unmapped-p/1
  (ensure (mate-unmapped-p #x0008)))

(addtest (cl-sam-tests) query-forward-p/1
  (ensure (query-forward-p #x0000))
  (ensure (not (query-forward-p #x0010))))

(addtest (cl-sam-tests) mate-forward-p/1
  (ensure (mate-forward-p #x0000))
  (ensure (not (mate-forward-p #x0020))))

(addtest (cl-sam-tests) first-in-pair-p/1
  (ensure (first-in-pair-p #x0040)))

(addtest (cl-sam-tests) second-in-pair-p/1
  (ensure (second-in-pair-p #x0080)))

(addtest (cl-sam-tests) alignment-not-primary-p/1
  (ensure (alignment-not-primary-p #x0100)))

(addtest (cl-sam-tests) fails-platform-qc-p/1
  (ensure (fails-platform-qc-p #x0200)))

(addtest (cl-sam-tests) pcr/optical-duplicate-p/1
  (ensure (pcr/optical-duplicate-p #x0400)))

(addtest (cl-sam-tests) bam-byte-round-trip/1
  (let ((in-filespec (namestring (merge-pathnames "data/c1215-name.bam")))
        (out-filespec (namestring
                       (tmp-pathname :tmpdir (merge-pathnames "data")
                                     :basename "bam-byte-round-trip-"))))
    (with-bgzf (in in-filespec :direction :input)
      (with-bgzf (out out-filespec :direction :output :if-exists :supersede)
        (multiple-value-bind (header num-refs ref-meta)
            (read-bam-meta in)
          (write-bam-meta out header num-refs ref-meta)
          (loop
             for aln = (read-alignment in)
             while aln
             do (write-alignment out aln)))))
    (test-binary-files-eql in-filespec out-filespec)
    (with-bgzf (in in-filespec :direction :input)
      (ensure (bgzf-eof-p in)))
    (with-bgzf (out out-filespec :direction :input)
      (ensure (bgzf-eof-p out)))
    (delete-file out-filespec)))

(addtest (cl-sam-tests) bam-uncompressed-round-trip/1
  (let ((in-filespec (namestring (merge-pathnames "data/c1215-name.bam")))
        (out-filespec (namestring
                       (tmp-pathname :tmpdir (merge-pathnames "data")
                                     :basename "bam-uncompressed-roundtrip-"))))
    (with-bam (in (header num-refs ref-meta) in-filespec)
      (with-bam (out (header num-refs ref-meta) out-filespec
                     :direction :output :if-exists :supersede :compression 0)
        (loop
           while (has-more-p in)
           for aln = (next in)
           do (consume out aln))))
    (test-sam-equal in-filespec out-filespec)
    (delete-file out-filespec)))

(addtest (cl-sam-tests) bam-uncompressed-header/1
  (let ((in-filespec (namestring (merge-pathnames "data/c1215-name.bam")))
        (out-filespec (namestring
                       (tmp-pathname :tmpdir (merge-pathnames "data")
                                     :basename "bam-uncompressed-header-roundtrip-"))))
    (with-bam (in (header num-refs ref-meta) in-filespec)
      (with-bam (out (header num-refs ref-meta) out-filespec
                     :direction :output :if-exists :supersede :compression 5
                     :compress-header nil)
        (loop
           while (has-more-p in)
           for aln = (next in)
           do (consume out aln))))
    (test-sam-equal in-filespec out-filespec)
    (delete-file out-filespec)))

(addtest (cl-sam-tests) bam-null-padded-header/1
  (let ((in-filespec (namestring (merge-pathnames "data/c1215-name.bam")))
        (out-filespec (namestring
                       (tmp-pathname :tmpdir (merge-pathnames "data")
                                     :basename "bam-null-padded-header-roundtrip-"))))
    (with-bam (in (header num-refs ref-meta) in-filespec)
      (with-bam (out (header num-refs ref-meta) out-filespec
                     :direction :output :if-exists :supersede :compression 5
                     :compress-header nil :pad-header 1024)
        (loop
           while (has-more-p in)
           for aln = (next in)
           do (consume out aln))))
    (test-sam-equal in-filespec out-filespec)
    (delete-file out-filespec)))

(addtest (cl-sam-tests) bam-semantic-round-trip/1
  (with-bgzf (bgzf (merge-pathnames "data/c1215-name.bam") :direction :input)
    (multiple-value-bind (header num-refs ref-meta)
        (read-bam-meta bgzf)
      (loop
         for aln = (read-alignment bgzf)
         while aln
         do (let ((aln2 (make-alignment-record
                         (read-name aln) (seq-string aln)
                         (alignment-flag aln)
                         :reference-id (reference-id aln)
                         :alignment-pos (alignment-position aln)
                         :mate-reference-id (mate-reference-id aln)
                         :mate-alignment-pos (mate-alignment-position aln)
                         :mapping-quality (mapping-quality aln)
                         :alignment-bin (alignment-bin aln)
                         :insert-length (insert-length aln)
                         :cigar (alignment-cigar aln)
                         :quality-str (quality-string aln)
                         :tag-values (alignment-tag-values aln))))
              (loop
               for i from 0 below (length aln)
               do (ensure (= (aref aln i) (aref aln2 i))
                          :report "different bytes at ~d: ~a -> ~a in ~a ~a"
                          :arguments (i (aref aln i) (aref aln2 i) aln aln2)))
              (ensure (equalp aln aln2)
                      :report "expected ~a giving ~a but found ~a giving ~a"
                      :arguments (aln (alignment-core-alist aln)
                                  aln2 (alignment-core-alist aln2))))))))

(addtest (cl-sam-tests) alignment-record</1
  ;; Unmapped (no reference) sort last
  (ensure (alignment-record<
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :first-in-pair)
                                  :alignment-pos 100
                                  :reference-id 0)
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :second-in-pair)
                                  :alignment-pos 1))))

(addtest (cl-sam-tests) alignment-record</2
  ;; By reference
  (ensure (alignment-record<
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :first-in-pair)
                                  :alignment-pos 100
                                  :reference-id 0)
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :second-in-pair)
                                  :alignment-pos 1
                                  :reference-id 1))))

(addtest (cl-sam-tests) alignment-record</3
  ;; By position
  (ensure (alignment-record<
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :first-in-pair)
                                  :reference-id 1
                                  :alignment-pos 1)
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :second-in-pair)
                                  :reference-id 1
                                  :alignment-pos 100))))

(addtest (cl-sam-tests) alignment-record</4
  ;; By strand
  (ensure (alignment-record<
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :first-in-pair
                                                          :query-forward)
                                  :reference-id 1
                                  :alignment-pos 1)
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :second-in-pair
                                                          :query-reverse)
                                  :reference-id 1
                                  :alignment-pos 1))))

(addtest (cl-sam-tests) alignment-name</1
  ;; By name
  (ensure (alignment-name<
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :first-in-pair))
           (make-alignment-record "aab" "acgt" (flag-bits 0 :sequenced-pair
                                                          :first-in-pair)))))

(addtest (cl-sam-tests) alignment-name</2
  ;; Fallback to template region
  (ensure (alignment-name<
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :first-in-pair)
                                  :alignment-pos 1)
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :second-in-pair)
                                  :alignment-pos 1))))

(addtest (cl-sam-tests) alignment-name</3
  ;; Fallback to strand
  (ensure (alignment-name<
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :first-in-pair
                                                          :query-forward)
                                  :alignment-pos 1)
           (make-alignment-record "aaa" "acgt" (flag-bits 0 :sequenced-pair
                                                          :first-in-pair
                                                          :query-reverse)
                                  :alignment-pos 1))))

;; FIXME -- the correct exhaustive sort order is still not defined by
;; the SAM spec. We will have to define our own.
(addtest (cl-sam-tests) sort-bam-file/1
  (let ((unsorted (merge-pathnames "data/c1215-name.bam"))
        (sorted (merge-pathnames "data/c1215-coordinate.bam")))
    (sort-bam-file unsorted sorted :sort-order :coordinate :buffer-size 10000)
    (ensure (fad:file-exists-p sorted))))

(addtest (cl-sam-tests) valid-read-name-p/1
  (ensure (sam::valid-read-name-p "a"))
  (ensure (not (sam::valid-read-name-p "@")))
  (ensure (not (sam::valid-read-name-p "a@")))
  (ensure (sam::valid-read-name-p "a"))
  (ensure (sam::valid-read-name-p (make-string 255 :initial-element #\a)))
  (ensure (not (sam::valid-read-name-p (make-string 256
                                                    :initial-element #\a)))))
