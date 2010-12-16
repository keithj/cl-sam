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

(defun make-bam-output (bam)
  "Returns a consumer function that accepts an argument of a BAM
record and writes it to BAM output stream BAM. The standard consumer
interface function CONSUME may be used in operations on the returned
consumer."
  (lambda (aln)
    (write-alignment bam aln)))

(defun write-bam-magic (bgzf &key (compress t))
  "Writes the BAM magic number to the handle BGZF."
  (write-bytes bgzf *bam-magic* (length *bam-magic*) :compress compress))

(defun write-bam-header (bgzf header &key (compress t) (null-padding 0)
                         (mtime 0))
  "Writes the BAM header string HEADER to handle BGZF, followed by
padding of NULL-PADDING null bytes. This function also writes the
header length, including padding, in the 4 bytes preceding the
header."
  (let* ((hlen (length header))
         (buffer (make-array (+ 4 hlen null-padding) :element-type 'octet
                             :initial-element +null-byte+)))
    ;; store header length, including padding
    (encode-int32le (+ hlen null-padding) buffer)
    (when (plusp hlen)
      (let ((hcodes (make-array hlen :element-type 'octet)))
        (map-into hcodes #'char-code header)
        (replace buffer hcodes :start1 4)))
    (write-bytes bgzf buffer (length buffer) :compress compress :mtime mtime)))

(defun write-num-references (bgzf n)
  "Writes the number of reference sequences N to handle BGZF."
  (let ((buffer (make-array 4 :element-type 'octet)))
    (encode-int32le n buffer)
    (write-bytes bgzf buffer 4)))

(defun write-reference-meta (bgzf ref-name ref-length)
  "Writes the metadata for a single reference sequence named REF-NAME,
of length REF-LENGTH bases, to handle BGZF."
  (let* ((name-len (length ref-name))
         (bam-name-len (1+ name-len))
         (name-codes (make-array name-len :element-type 'octet))
         (buffer-len (+ 4 bam-name-len 4))
         (name-offset 4)
         (ref-len-offset (+ name-offset bam-name-len))
         (buffer (make-array buffer-len :element-type 'octet
                             :initial-element +null-byte+)))
    (encode-int32le bam-name-len buffer)
    (map-into name-codes #'char-code ref-name)
    (replace buffer name-codes :start1 name-offset)
    (encode-int32le ref-length buffer ref-len-offset)
    (write-bytes bgzf buffer buffer-len)))

(declaim (ftype (function (bgzf simple-octet-vector) fixnum) write-alignment))
(defun write-alignment (bgzf alignment-record)
  "Writes one ALIGNMENT-RECORD to handle BGZF and returns the number
of bytes written."
  (declare (optimize (speed 3)))
  (let ((alen (length alignment-record))
        (alen-bytes (make-array 4 :element-type 'octet :initial-element 0)))
    (encode-int32le alen alen-bytes)
    (+ (write-bytes bgzf alen-bytes 4)
       (write-bytes bgzf alignment-record alen))))

(defun write-bam-meta (bgzf header num-refs ref-meta &key (compress t)
                       (null-padding 0))
  "Writes BAM magic number and then all metadata to handle BGZF. The
metadata consist of the HEADER string, number of reference sequences
NUM-REFS and a list reference sequence metadata REF-META. The list
contains one element per reference sequence, each element being a list
of reference identifier, reference name and reference length. Returns
the number of bytes written.

Optional:

- null-padding (fixnum): A number of null bytes to be appended to the
end of the header string, as allowed by the SAM spec. This is useful
for creating slack space so that BAM headers may be edited in place."
  (+ (write-bam-magic bgzf :compress compress)
     (prog1
         (write-bam-header bgzf header :compress compress
                           :null-padding null-padding)
       ;; Unless compressing, flush here to move to the next bgz
       ;; member, leaving the magic number and header uncompressed and
       ;; in their own block.
       (unless compress
         (bgzf-flush bgzf :compress compress :append-eof nil)))
     (write-num-references bgzf num-refs)
     (loop
        for (nil ref-name ref-length) in ref-meta
        sum (write-reference-meta bgzf ref-name ref-length))))
