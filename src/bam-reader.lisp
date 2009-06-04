;;;
;;; Copyright (C) 2009 Keith James. All rights reserved.
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

(defconstant +bam-magic+ #(66 65 77 1)
  "The BAM file magic header bytes.")

(defun read-bam-magic (bgzf)
  "Reads the BAM magic number from the handle BGZF and returns T if it
is valid or raises a {define-condition malformed-file-error} if not."
  (let ((magic (read-bytes bgzf 4)))
    (unless (equalp +bam-magic+ magic)
      (error 'malformed-file-error :text "invalid BAM magic number"))
    t))

(defun read-bam-header (bgzf)
  "Returns the unparsed BAM header from the handle BGZF as a Lisp
string."
  (let ((len (decode-int32le (read-bytes bgzf 4))))
    (when (plusp len)
      (read-string bgzf len))))

(defun read-num-references (bgzf)
  "Returns the number of reference sequences described in the BAM
header of handle BGZF."
  (decode-int32le (read-bytes bgzf 4)))

(defun read-reference-meta (bgzf)
  "Returns the reference sequence metadata for a single reference
sequence described in the BAM header of handle BGZF. Two values are
returned, the reference name as a string and the reference length as
an integer,"
  (let ((len (decode-int32le (read-bytes bgzf 4))))
    (values (read-string bgzf len :null-terminated t)
            (decode-int32le (read-bytes bgzf 4)))))

(defun read-alignment (bgzf)
  "Reads one alignment block from handle BGZF, returns it as a Lisp
array of unsigned-byte 8. The handle is advanced to the next
alignment. If no more alignments are available, returns NIL."
  (let ((size-bytes (read-bytes bgzf 4)))
    (if (null size-bytes)
        nil
      (let ((block-size (decode-int32le size-bytes)))
        (read-bytes bgzf block-size)))))
