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

(in-package :cl-user)

(defpackage :bgzf-ffi
  (:use #:common-lisp #:cffi)
  (:export #:bgzf-open
           #:bgzf-close
           #:bgzf-read
           #:bgzf-write
           #:bgzf-seek
           #:bgzf-tell)
  (:documentation "A foreign function interface to the block gzip
  functions provided by the Broad Institute. See bgzf.h in the
  SAMtools package for function documentation."))

(defpackage :sam
  (:use #:common-lisp #:cffi #:trivial-gray-streams
        #:deoxybyte-utilities #:deoxybyte-io)
  (:export

   ;; Conditions
   #:bgzf-io-error

   ;; BGZF handle API
   #:bgzf
   #:bzgf-p
   #:bgzf-file
   #:bgzf-open-p
   #:bgzf-open
   #:bgzf-close
   #:bgzf-seek
   #:bgzf-tell
   #:with-bgzf-file
   
   ;; BGZF stream API
   #:bgzf-stream
   #:bgzf-input-stream
   #:bgzf-stream-open
   
   ;; Low-level BAM reading API
   #:read-bam-magic
   #:read-bam-header
   #:read-num-references
   #:read-reference-meta
   #:read-alignment

   ;; Low-level BAM API
   #:reference-id
   #:alignment-position
   #:read-name-length
   #:mapping-quality
   #:alignment-bin
   #:cigar-length
   #:alignment-flag
   #:read-length
   #:mate-reference-id
   #:mate-alignment-position
   #:insert-length

   #:read-name
   #:alignment-cigar
   #:seq-string
   #:qual-string
   #:alignment-tag-values

   #:sequenced-pair-p
   #:mapped-proper-pair-p
   #:query-unmapped-p
   #:mate-unmapped-p
   #:query-forward-p
   #:mate-forward-p
   #:first-in-pair-p
   #:second-in-pair-p
   #:alignment-not-primary-p
   #:fails-platform-qc-p
   #:pcr/optical-duplicate-p
   
   #:alignment-core
   #:alignment-core-alist
   #:alignment-flag-alist)
  (:documentation ""))
