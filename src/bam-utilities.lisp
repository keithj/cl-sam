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

(in-package :sam)

(defun repair-mapping-flags (aln)
  "Returns alignment ALN if it has correct mapped-proper-pair,
query-mapped and mate-mapped flags, or a modified copy with fixed
flags. Also sets the user tag ZF:I:<original flag> if the flags are
changed.

Alignments produced by the BWA have been observed to contain invalid
flags. This is caused by BWA creating mappings that overhang the end
of the reference. The underlying cause is reference concatenation in
the Burrows-Wheeler index.

This function attempts to fix invalid mapped-proper-pair flags in
cases where query-unmapped and/or mate-unmapped flags are set. Such
unmapped reads may also have reference-ids set and good mapping
scores. This function does not correct these other fields.

It is necessary to use this function in order to obtain correct flag
counts when using samtools flagstat or {defun flagstat} ."
  (let ((flag (alignment-flag aln)))
    (cond ((not (sequenced-pair-p flag))
           aln)
          ((and (mapped-proper-pair-p flag)
                (or (query-unmapped-p flag) (mate-unmapped-p flag)))
           (copy-alignment-record
            aln :alignment-flag (flag-bits flag :not-mapped-proper-pair)
            :tag-values (acons :zf flag (alignment-tag-values aln))))
          (t
           aln))))
