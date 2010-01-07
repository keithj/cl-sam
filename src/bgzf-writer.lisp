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

(defvar *bgz-write-buffer* (make-array 4 :element-type 'octet
                                       :initial-element 0)
  "The buffer used by {defun write-bgz-member} for writing block gzip
data. Rebind this per thread to make {defun write-bgz-member}
reentrant.")

(declaim (ftype (function (bgzf simple-octet-vector vector-index
                                &key (:compress t)) fixnum) write-bytes))
(defun write-bytes (bgzf bytes n &key (compress t))
  (declare (optimize (speed 3) (safety 1)))
  (let ((stream (bgzf-stream bgzf)))
    (declare (type simple-octet-vector bytes)
             (type vector-index n))
    (labels ((buffer-full-p ()
               (= (bgzf-pointer bgzf) (1- (length (bgzf-buffer bgzf)))))
             (store-overflow (num-deflated)
               (let* ((buf (bgzf-buffer bgzf))
                      (overflow (subseq buf num-deflated (bgzf-pointer bgzf))))
                 (replace buf overflow)
                 (setf (bgzf-pointer bgzf) (length overflow))))
             (deflate-to-bgz ()
               (let ((buf (bgzf-buffer bgzf))
                     (cdata (make-array +bgz-max-payload-length+
                                        :element-type 'octet)))
                 (multiple-value-bind (deflated bytes-in bytes-out)
                     (gz:deflate-vector buf cdata
                       :compression (if compress
                                        (bgzf-compression bgzf)
                                      0)
                       :suppress-header t :window-bits 15 :backoff 1024)
                   (declare (ignore deflated))
                   (declare (type bgz-payload-index bytes-in bytes-out))
                   (let* ((udata (subseq buf 0 bytes-in))
                          (bgz (make-bgz-member
                                :mtime (get-universal-time) :xlen +xlen+
                                :udata udata :cdata cdata :cend bytes-out
                                :bsize (+ bytes-out +member-header-length+
                                          +member-footer-length+)
                                :isize bytes-in
                                :crc32 (gz:crc32 udata))))
                     (write-bgz-member bgz stream)
                     (if (< bytes-in (length buf)) ; Didn't fit
                         (store-overflow bytes-in)
                       (setf (bgzf-pointer bgzf) 0))
                     buf)))))
      (cond ((< (+ (bgzf-pointer bgzf) n) (1- (length (bgzf-buffer bgzf))))
             (replace (bgzf-buffer bgzf) bytes :start1 (bgzf-pointer bgzf))
             (incf (bgzf-pointer bgzf) n)
             (when (buffer-full-p)
               (deflate-to-bgz))
             n)
            (t
             (loop
                with buffer = (bgzf-buffer bgzf)
                for i = 0 then (1+ i)
                while (< i n)
                do (progn
                     (setf (aref buffer (bgzf-pointer bgzf)) (aref bytes i))
                     (if (buffer-full-p)
                         (setf buffer (deflate-to-bgz))
                       (incf (bgzf-pointer bgzf))))
                finally (return n)))))))

(defun bgzf-flush (bgzf &key (compress t) (append-eof t))
  (if (plusp (bgzf-pointer bgzf))
      (let ((buffer (subseq (bgzf-buffer bgzf) 0 (bgzf-pointer bgzf)))
            (cdata (make-array +bgz-max-payload-length+ :element-type 'octet)))
        (multiple-value-bind (deflated bytes-in bytes-out)
            (gz:deflate-vector buffer cdata
              :compression (if compress
                               (bgzf-compression bgzf)
                             0)
              :suppress-header t :window-bits 15 :backoff 1024)
          (declare (ignore deflated))
          (let* ((udata (subseq buffer 0 bytes-in))
                 (bgz (make-bgz-member :mtime (get-universal-time)
                                       :xlen +xlen+ :udata udata
                                       :cdata cdata :cend bytes-out
                                       :bsize (+ bytes-out
                                                 +member-header-length+
                                                 +member-footer-length+)
                                       :isize bytes-in
                                       :crc32 (gz:crc32 udata))))
            (write-bgz-member bgz (bgzf-stream bgzf))
            (if (< bytes-in (length buffer)) ; Didn't fit
                (let ((overflow (subseq buffer bytes-in (bgzf-pointer bgzf))))
                  (replace buffer overflow)
                  (setf (bgzf-pointer bgzf) (length overflow))
                  (bgzf-flush bgzf)) ; Flush again to write the overflow
              (setf (bgzf-pointer bgzf) 0))))))
  (when append-eof
    (write-sequence *empty-bgz-record* (bgzf-stream bgzf))))

(defun write-bgz-member (bgz stream &optional (buffer *bgz-write-buffer*))
  "Writes one BGZ member to STREAM.

Arguments:

- stream (octet output-stream): An open stream.

Optional:

- buffer (simple-octet-vector): A re-usable read buffer which must be
  able to contain at least 4 bytes. Defaults to *BGZ-WRITE-BUFFER*
  . This argument is only required where the the function calls must
  be reentrant (normally achieved by rebinding *BGZ-WRITE-BUFFER* per
  thread).

Returns:

- The number of bytes written."
  (declare (optimize (speed 3)))
  (flet ((encode-bytes (x n s)
           (ecase n
             (1 (encode-int8le x buffer))
             (2 (encode-int16le x buffer))
             (4 (encode-int32le x buffer)))
           (write-sequence buffer s :end n)))
    (let ((bsize (bgz-member-bsize bgz)))
      (encode-bytes (bgz-member-id1 bgz) 1 stream)
      (encode-bytes (bgz-member-id2 bgz) 1 stream)
      (encode-bytes (bgz-member-cm bgz) 1 stream)
      (encode-bytes (bgz-member-flg bgz) 1 stream)
      (encode-bytes (get-universal-time) 4 stream)
      (encode-bytes (bgz-member-xfl bgz) 1 stream)
      (encode-bytes gz:+os-unknown+ 1 stream)
      (encode-bytes (bgz-member-xlen bgz) 2 stream)
      (encode-bytes (bgz-member-sf1 bgz) 1 stream)
      (encode-bytes (bgz-member-sf2 bgz) 1 stream)
      (encode-bytes (bgz-member-slen bgz) 2 stream)
      (encode-bytes (1- bsize) 2 stream) ; 18 header bytes
      (write-sequence (bgz-member-cdata bgz) stream
                      :end (bgz-member-cend bgz)) ; Payload
      (encode-bytes (bgz-member-crc32 bgz) 4 stream)
      (encode-bytes (bgz-member-isize bgz) 4 stream) ; 8 footer bytes
      bsize)))
