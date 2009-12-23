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

(defconstant +xlen+ 6)
(defconstant +sf1+ (char-code #\B))
(defconstant +sf2+ (char-code #\C))
(defconstant +slen+ 2)

(defconstant +member-header-length+ 18)
(defconstant +member-footer-length+ 8)

(defvar *bgz-read-buffer* (make-array 8 :element-type '(unsigned-byte 8)
                                      :initial-element 0)
  "The buffer used by {defun read-bgz-member} for reading block gzip
data. Rebind this per thread to make {defun read-bgz-member}
re-entrant.")

(defvar *bgz-write-buffer* (make-array 8 :element-type '(unsigned-byte 8)
                                       :initial-element 0)
  "The buffer used by {defun write-bgz-member} for writing block gzip
data. Rebind this per thread to make {defun write-bgz-member}
re-entrant.")

(defstruct (bgz-member (:include gz:gz-member))
  "A gzip member with extensions, as defined by RFC1952. The
extensions are described in RFC1952 and the SAM format specification."
  (sf1 +sf1+ :type uint8)
  (sf2 +sf2+ :type uint8)
  (slen +slen+ :type uint16)
  (bsize 0 :type uint16)
  (udata (make-array 0 :element-type '(unsigned-byte 8))
         :type (simple-array (unsigned-byte 8) (*))))

(defun read-bgz-member (stream &optional (buffer *bgz-read-buffer*))
  (flet ((decode-bytes (s n)
           (when (plusp (read-sequence buffer s :end n))
             (ecase n
               (1 (decode-uint8le buffer))
               (2 (decode-uint16le buffer))
               (4 (decode-uint32le buffer))))))
    (let ((id1 (decode-bytes stream 1)))
      (when id1
        (let* ((id2 (decode-bytes stream 1))
               (cm (decode-bytes stream 1))
               (flg (decode-bytes stream 1))
               (mtime (decode-bytes stream 4))
               (xfl (decode-bytes stream 1))
               (os (decode-bytes stream 1))
               (xlen (decode-bytes stream 2))
               (sf1 (decode-bytes stream 1))
               (sf2 (decode-bytes stream 1))
               (slen (decode-bytes stream 2))
               (bsize (decode-bytes stream 2))) ; 18 header bytes
          (declare (ignore xfl))
          (assert (and (= gz:+id1+ id1)
                       (= gz:+id2+ id2)
                       (= gz:+cm-deflate+ cm)
                       (not (zerop (logand gz:+flag-extra+ flg)))
                       (= +xlen+ xlen)
                       (= +sf1+ sf1)
                       (= +sf2+ sf2)
                       (= +slen+ slen)) () "Invalid block gzip header")
          (let* ((deflated-size (1+ (logand bsize #xffff)))
                 (cdata (make-array (- deflated-size
                                       +member-header-length+
                                       +member-footer-length+)
                                    :element-type '(unsigned-byte 8))))
            (read-sequence cdata stream)
            (let* ((crc32 (decode-bytes stream 4))
                   (isize (decode-bytes stream 4))) ; 8 footer bytes
              (make-bgz-member :mtime mtime :os os :xlen xlen
                               :bsize deflated-size
                               :cdata cdata :crc32 crc32 :isize isize))))))))

(defun write-bgz-member (bgz stream &optional (buffer *bgz-write-buffer*))
  (flet ((encode-bytes (x s)
           (let ((n (etypecase x
                      ((unsigned-byte 8) (encode-int8le x buffer) 1)
                      ((unsigned-byte 16) (encode-int16le x buffer) 2)
                      ((unsigned-byte 32) (encode-int32le x buffer) 4))))
             (write-sequence buffer s :end n))))
    ;; Header
    (encode-bytes (bgz-member-id1 bgz) stream)
    (encode-bytes (bgz-member-id2 bgz) stream)
    (encode-bytes (bgz-member-cm bgz) stream)
    (encode-bytes (bgz-member-flg bgz) stream)
    (encode-bytes (get-universal-time) stream)
    (encode-bytes (bgz-member-xfl bgz) stream)
    ;; FIXME -- I have no experience reliably detecting OS type. I'm
    ;; sure this can be improved.
    (encode-bytes  #+:unix gz:+os-unix+
                   #+:windows gz:+os-ntfs-filesystem+
                   #-(or :unix :windows) gz:+os-unknown+ stream)
    (encode-bytes (bgz-member-xlen bgz) stream)
    (encode-bytes (bgz-member-sf1 bgz) stream)
    (encode-bytes (bgz-member-sf2 bgz) stream)
    (encode-bytes (bgz-member-slen bgz) stream)
    (encode-bytes (bgz-member-bsize bgz) stream)
    ;; Payload
    (write-sequence (bgz-member-cdata bgz) stream)
    ;; Footer
    (encode-bytes (bgz-member-crc32 bgz) stream)
    (encode-bytes (bgz-member-isize bgz) stream)))


(defun test (filespec1 filespec2)
  (with-open-file (stream1 filespec1 :element-type '(unsigned-byte 8))
    (with-open-file (stream2 filespec2 :element-type '(unsigned-byte 8)
                             :direction :output :if-exists :supersede)
      (loop
         for bgz = (read-bgz-member stream1)
         while bgz
         do (let* ((cdata (bgz-member-cdata bgz))
                   (isize (bgz-member-isize bgz))
                   (data (make-array isize :element-type '(unsigned-byte 8))))
              (write-sequence
               (gz:inflate-vector cdata data
                                  :suppress-header t :window-bits 15) stream2))
         count bgz))))
