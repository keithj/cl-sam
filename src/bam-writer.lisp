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

(in-package :sam)

(defun write-bam-magic (bgzf)
  "Writes the BAM magic number to the handle BGZF."
  (write-bytes bgzf +bam-magic+ (length +bam-magic+)))

(defun write-bam-header (bgzf header &optional (null-padding 0))
  "Writes the BAM header string HEADER to handle BGZF, followed by
padding of NULL-PADDING null bytes. This function also writes the
header length, including padding, in the 4 bytes preceding the
header."
  (let* ((hlen (length header))
         (buffer (make-array (+ 4 hlen null-padding)
                             :element-type '(unsigned-byte 8)
                             :initial-element +null-byte+)))
    ;; store header length, including padding
    (encode-int32le (+ hlen null-padding) buffer)
    (copy-array header 0 (1- hlen) buffer 4 #'char-code)
    (write-bytes bgzf buffer (length buffer))))

(defun write-num-references (bgzf n)
  "Writes the number of reference sequences N to handle BGZF."
  (let ((buffer (make-array 4 :element-type '(unsigned-byte 8))))
    (encode-int32le n buffer)
    (write-bytes bgzf buffer 4)))

(defun write-reference-meta (bgzf ref-name ref-length)
  "Writes the metadata for a single reference sequence named REF-NAME,
of length REF-LENGTH bases, to handle BGZF."
  (let* ((nlen (length ref-name))
         (blen (+ 4 nlen 1 4))
         (name-offset 4)
         (rlen-offset (+ name-offset nlen))
         (buffer (make-array blen :element-type '(unsigned-byte 8)
                             :initial-element +null-byte+)))
    (encode-int32le nlen buffer)
    (copy-array ref-name 0 (1- nlen) buffer name-offset #'char-code)
    (encode-int32le ref-length buffer rlen-offset)
    (write-bytes bgzf buffer blen)))
