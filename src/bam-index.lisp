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

(defconstant +max-num-bins+ (/ (1- (ash 1 18)) 7)
  "The maximum number of bins possible in the BAM indexing scheme.")
(defconstant +samtools-kludge-bin+ 37450
  "Extra bin with different semantics added by samtools.")

(defconstant +linear-bin-size+ (expt 2 14)
  "The size in bases of intervals in the linear bin index.")

(defvar *bam-index-magic* (make-array 4 :element-type 'octet
                                      :initial-contents '(66 65 73 1))
  "The BAI index file magic header bytes.")

(defparameter *tree-deepening-boundaries* '(4681 585  73   9   1))
(defparameter *voffset-merge-distance* (expt 2 14)
  "If two index chunks are this number of bytes or closer to each
other, they should be merged.")

(defstruct bam-index
  "A BAM file index covering all the reference sequences in a BAM
file."
  (refs (make-array 0) :type simple-vector))

(defstruct (samtools-bam-index (:include bam-index))
  "A samtools-specific index containing extra, undocumented data:
 - unassigned: number of unmapped reads not assigned to a reference"
  (unassigned 0 :type (unsigned-byte 64)))

(defstruct (ref-index (:print-object print-ref-index))
  "An index for a single reference sequence. An index is composed of a
vector of bins (the binning index) and a vector of intervals (the
linear index). The bins are sorted by increasing bin number. The
binning index is hierarchical, while the linear index is flat, being a
projection of reads from all bins onto a single vector."
  (id -1 :type fixnum :read-only t)
  (bins (make-array 0 :initial-element nil) :type simple-vector)
  (intervals (make-array 0 :element-type 'fixnum)
             :type (simple-array fixnum (*))))

(defstruct (samtools-ref-index (:include ref-index)
                               (:print-object print-samtools-ref-index))
  "A samtools-specific reference index containing extra,
undocumented data:

 - start: the start offset of the reference
 - end: the end offset of the reference
 - mapped: number of reads mapped to the reference
 - unmapped: number of unmapped reads assigned to the reference by
   magic"
  (start 0 :type (unsigned-byte 64))
  (end 0 :type (unsigned-byte 64))
  (mapped 0 :type (unsigned-byte 64))
  (unmapped 0 :type (unsigned-byte 64)))

(defstruct (bin (:print-object print-bin))
  "A bin within a BAM binning index. A bin is a region on a reference
sequence that is identified by a number. The bin location and
numbering system used is the UCSC binning scheme by Richard Durbin and
Lincoln Stein. A bin contains a vector of chunks."
  (num 0 :type fixnum :read-only t)
  (chunks (make-array 0) :type simple-vector))

(defstruct (chunk (:print-object print-chunk))
  "A chunk within a BAM index bin. Chunks are groups of records that
lie within a bin and are close together within the BAM file."
  (start 0 :type (unsigned-byte 64))
  (end most-positive-fixnum :type (unsigned-byte 64)))

(defstruct (region (:constructor make-region (&key ref start end))
                   (:constructor region (ref start end))
                   (:print-object print-region))
  "A range of bases over a reference sequence."
  (ref nil :type t)
  (start 0 :type fixnum)
  (end 0 :type fixnum))

(defun print-ref-index (ref-index stream)
  "Prints a string representation of REF-INDEX to STREAM."
  (print-unreadable-object (ref-index stream :type t)
    (format stream "~d " (ref-index-id ref-index))
    (princ (remove-if #'null (ref-index-bins ref-index)) stream)
    (princ #\Space stream)
    (princ (run-length-encode (ref-index-intervals ref-index)) stream)))

(defun print-samtools-ref-index (ref-index stream)
  "Prints a string representation of REF-INDEX to STREAM."
  (print-unreadable-object (ref-index stream :type t)
    (let ((start (samtools-ref-index-start ref-index))
          (end (samtools-ref-index-end ref-index)))
      (format stream "~d ~d:~d..~d:~d ~d/~d "
              (ref-index-id ref-index)
              (bgzf-coffset start) (bgzf-uoffset start)
              (bgzf-coffset end) (bgzf-uoffset end)
              (samtools-ref-index-mapped ref-index)
              (samtools-ref-index-unmapped ref-index)))
    (princ (remove-if #'null (ref-index-bins ref-index)) stream)
    (princ #\Space stream)
    (princ (run-length-encode (ref-index-intervals ref-index)) stream)))

(defun print-bin (bin stream)
  "Prints a string representation of BIN to STREAM."
  (print-unreadable-object (bin stream :type t)
    (princ (bin-num bin) stream)
    (princ #\Space stream)
    (princ (bin-chunks bin) stream)))

(defun print-chunk (chunk stream)
  "Prints a string representation of CHUNK to STREAM."
  (print-unreadable-object (chunk stream :type nil)
    (let ((start (chunk-start chunk))
          (end (chunk-end chunk)))
      (format stream "~d:~d..~d:~d"
              (bgzf-coffset start) (bgzf-uoffset start)
              (bgzf-coffset end) (bgzf-uoffset end)))))

(defun print-region (region stream)
  "Prints a string representation of RANGE to STREAM."
  (print-unreadable-object (region stream :type t)
    (format stream "ref: ~a ~d..~d" (region-ref region)
            (region-start region) (region-end region))))

(defun ref-index (bam-index reference-id)
  "Returns the reference index for reference number REFERENCE-ID in
BAM-INDEX."
  (svref (bam-index-refs bam-index) reference-id))

(defun ref-index-bin (ref-index bin-num)
  "Returns the BIN number BIN-NUM in REF-INDEX."
  (binary-search (ref-index-bins ref-index) bin-num :key #'bin-num))

(defun empty-bin-p (bin)
  "Returns T if BIN contains no chunks, or NIL otherwise."
  (zerop (length (bin-chunks bin))))

(defun bin-chunk (bin chunk-num)
  "Returns the chunk number CHUNK-NUM in BIN."
  (svref (bin-chunks bin) chunk-num))

(defun max-bin-num (ref-length)
  (+ (1- (first *tree-deepening-boundaries*)) (ash ref-length -14)))

(defun region-to-bins (start end)
  "Returns a bit vector with the bits for relevant bins set, given a
0-based range. Each bin may span 2^29, 2^26, 2^23, 2^20, 2^17 or 2^14
bp. Bin 0 spans a 512 Mbp region, bins 1-8 span 64 Mbp, 9-72 8 Mbp,
73-584 1 Mbp, 585-4680 128 Kbp and bins 4681-37449 span 16 Kbp
regions."
  (check-arguments (<= 0 start end) (start end) "must satisfy (<= 0 start end)")
  ;; Offsets calculated thus:
  ;; (mapcar (lambda (x)
  ;;           (/ (1- (ash 1 x)) 7)) '(15 12 9 6 3))
  (let ((bins (make-array +max-num-bins+ :element-type 'bit
                          :initial-element 0)))
    (when (< start end)
      (let ((end (1- end)))
        (loop
           for  shift in '( -14 -17 -20 -23 -26)
           for offset in '(4681 585  73   9   1) ; tree deepening boundaries
           do (loop
                 for k from (+ offset (ash start -14))
                 to (+ offset (ash end shift))
                 do (setf (sbit bins k) 1)))))
    bins))

(defun region-to-bin (start end)
  "Returns a bin number given a 0-based range."
  (check-arguments (<= 0 start end) (start end) "must satisfy (<= 0 start end)")
  (macrolet ((define-bins (offsets shifts)
               `(cond ,@(mapcar (lambda (offset shift)
                                  `((= (ash start ,shift)
                                       (ash end ,shift))
                                    (+ ,(/ (1- (ash 1 offset)) 7)
                                       (ash end ,shift)))) offsets shifts)
                      (t
                       0))))
    (let ((end (1- end)))
      (define-bins
          ( 15  12   9   6   3)
          (-14 -17 -20 -23 -26)))))

(defun region-to-interval (coord)
  "Given a coordinate COORD, returns the array index of its interval
in the BAM linear index."
  (floor coord +linear-bin-size+))

(defun find-bins (index reference-id start end)
  "Returns a list of all the bins for reference REFERENCE-ID in INDEX,
between reference positions START and END."
  (let ((num-refs (length (bam-index-refs index))))
    (check-arguments (< -1 reference-id num-refs) (reference-id)
                     "must satisfy (<= 0 reference-id ~d)" (1- num-refs))
    (check-arguments (<= start end) (start end)
                     "must satisfy (<= start end)"))
 (let ((bin-nums (region-to-bins start end))
       (bins (ref-index-bins (ref-index index reference-id))))
   (loop
      for i from 0 below +max-num-bins+
      for bin = (when (plusp (sbit bin-nums i))
                  (binary-search bins i :key #'bin-num))
      when (and bin (not (empty-bin-p bin)))
      collect bin)))

(defun adjacentp (chunk1 chunk2)
  "Returns T if BAM index CHUNK1 and CHUNK2 are adjacent. In this
context this means that they are sequential in a file stream and close
enough together for it to be more efficient to continue reading across
the intervening bytes, rather than seeking."
  (check-arguments (<= (chunk-end chunk1) (chunk-start chunk2))
                   (chunk1 chunk2)
                   "chunk1 must end at or before the start of chunk2")
  (voffset-merge-p (chunk-end chunk1) (chunk-start chunk2)))

(defun merge-chunks (chunks)
  "Returns a new list of BAM index chunks created by merging members
of list CHUNKS that are {defun adjacentp} ."
  (if (endp (rest chunks))
      chunks
      (let ((c1 (copy-chunk (first chunks)))
            merged)
        (dolist (c2 (rest chunks) (nreverse (cons c1 merged)))
          (cond ((adjacentp c1 c2)
                 (setf (chunk-end c1) (chunk-end c2)))
                (t
                 (push c1 merged)
                 (setf c1 (copy-chunk c2))))))))

(defun voffset-merge-p (voffset1 voffset2)
  "Returns T if BGZF virtual offsets are close enough to save disk
seeks. Continuing to read from the current position may be faster than
seeking forward a short distance."
  (let ((coffset1 (bgzf-coffset voffset1))
        (coffset2 (bgzf-coffset voffset2)))
    (check-arguments (<= coffset1 coffset2) (voffset1 voffset2)
                     "the first virtual offset must be at or before the second")
    (>= (+ coffset1 *voffset-merge-distance*) coffset2)))

;; TODO -- allow chunk queue to be created from a list of start/end
;; pairs on a reference. Collect all the relevant bins first.
(defun make-bam-chunk-queue (index reference-id start end)
  "Returns a queue of merged chunks for reference REFERENCE-ID in INDEX,
between reference positions START and END."
  (queue-append (make-queue)
                (merge-chunks
                 (sort (mapcan (lambda (bin)
                                 (coerce (bin-chunks bin) 'list))
                               (find-bins index reference-id start end))
                       #'< :key #'chunk-start))))

;; Ranges are (reference-id start end)
(defun make-bam-chunk-queue2 (index regions)
  (let* ((bins (delete-duplicates (mapcan (lambda (region)
                                            (apply #'find-bins index region))
                                          regions)))
         (chunks (mapcan (lambda (bin)
                           (coerce (bin-chunks bin) 'list))
                         bins)))
    (queue-append (make-queue)
                  (merge-chunks (sort chunks #'< :key #'chunk-start)))))

(defun make-bam-chunks (index region)
  (let* ((bins (find-bins index (region-ref region) (region-start region)
                          (region-end region)))
         (chunks (mapcan (lambda (bin)
                           (coerce (bin-chunks bin) 'list))
                         bins)))
    (sort chunks #'< :key #'chunk-start)))

(defun make-bam-chunk-input (bam chunks end)
  "Returns a generator function that uses a queue of index CHUNKS to
iterate over records in BAM stream BAM. The standard generator
interface functions, CURRENT, NEXT and HAS-MORE-P may be used in
operations on the returned generator."
  (let ((chunk (queue-dequeue chunks))
        (current (read-alignment bam)))
    (defgenerator
        (more (not (null current)))
        (next (prog1
                  current
                (if (null chunk)
                    (setf current nil)
                    (let ((aln (read-alignment bam)))
                      (setf current (cond ((null aln)
                                           nil)
                                          ((> (alignment-position aln) end)
                                           nil)
                                          (t
                                           aln)))
                      (when (>= (bgzf-tell bam) (chunk-end chunk))
                        (let ((next (queue-dequeue chunks)))
                          (when next
                            (bgzf-seek bam (chunk-start next)))
                          (setf chunk next))))))))))

(defun run-length-encode (intervals)
  "Returns a run-length encoded list representing INTERVALS, a BGAZ
linear index. Each element of the list has a car of the run-length and
a cdr of the BGZF file position. Used in printing text representations
of the index."
  (labels ((rle (x &optional current (run-length 0))
             (cond ((and (null x) (zerop run-length))
                    nil)
                   ((null x)
                    (list (cons run-length current)))
                   ((null current)
                    (rle (rest x) (first x) 1))
                     ((eql current (first x))
                      (rle (rest x) current (1+ run-length)))
                     (t
                      (cons (cons run-length current)
                            (rle (rest x) (first x) 1))))))
    (rle (coerce intervals 'list))))
