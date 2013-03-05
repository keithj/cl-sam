;;;
;;; Copyright (c) 2009-2013 Keith James. All rights reserved.
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

(defun write-chunked-bytes (bgzf bytes n &key (compress t) (mtime 0))
  (let ((chunk-size (- +bgz-max-payload-length+ 1024)))
    (if (> n chunk-size)
        (let ((num-chunks (ceiling n chunk-size)))
          (loop
             repeat num-chunks
             for start = 0 then (+ start chunk-size)
             for end = chunk-size then (min (+ end chunk-size) (length bytes))
             sum (write-bytes bgzf (subseq bytes start end) (- end start)
                              :compress compress :mtime mtime)))
        (write-bytes bgzf bytes n :compress compress :mtime mtime))))

(declaim (ftype (function (bgzf simple-octet-vector vector-index
                           &key (:compress t) (:mtime integer)) fixnum)
                write-bytes))
(defun write-bytes (bgzf bytes n &key (compress t) (mtime 0))
  "Write N elements from vector BYTES to stream BGZF.

Key:

- compress (boolean): If T, compress at current BGZF compression level
  when the current block becomes full. Defaults to T.
- mtime (uint32): The zlib mtime. Defaults to 0.

Returns:

- The number of bytes written."
  (declare (optimize (speed 3) (safety 1)))
  ;; (declare (optimize (debug 3) (safety 3)))
  (let ((stream (bgzf-stream bgzf))
        (compression (if compress
                         (bgzf-compression bgzf)
                         0))
        (deflate-backoff 1024))
    (declare (type simple-octet-vector bytes)
             (type vector-index n))
    (labels ((buffer-full-p ()
               (= (bgzf-pointer bgzf) (length (bgzf-buffer bgzf))))
             (store-overflow (num-deflated)
               (let* ((buf (bgzf-buffer bgzf))
                      (overflow (subseq buf num-deflated (bgzf-pointer bgzf))))
                 (replace buf overflow)
                 (setf (bgzf-pointer bgzf) (length overflow))))
             (deflate-to-bgz ()
               (let ((buf (bgzf-buffer bgzf))
                     (cdata (make-array +bgz-max-payload-length+
                                        :element-type 'octet
                                        :initial-element 0)))
                 (multiple-value-bind (deflated bytes-in bytes-out)
                     (gz:deflate-vector buf cdata
                       :compression compression :suppress-header t
                       :window-bits 15 :backoff deflate-backoff)
                   (declare (ignore deflated))
                   (declare (type (integer 0 65536) bytes-in bytes-out))
                   (let* ((udata (subseq buf 0 bytes-in))
                          (bgz (make-bgz-member
                                :mtime mtime :xlen +xlen+
                                :udata udata :cdata cdata
                                :cend bytes-out
                                :bsize (+ bytes-out
                                          +member-header-length+
                                          +member-footer-length+)
                                :isize bytes-in
                                :crc32 (gz:crc32 udata))))
                     (write-bgz-member bgz stream (bgzf-util-buffer bgzf))
                     (if (< bytes-in (length buf)) ; Didn't fit
                         (store-overflow bytes-in)
                         (setf (bgzf-pointer bgzf) 0))
                     buf)))))
      (cond ((< (+ (bgzf-pointer bgzf) n) (length (bgzf-buffer bgzf)))
             (replace (bgzf-buffer bgzf) bytes :start1 (bgzf-pointer bgzf))
             (incf (bgzf-pointer bgzf) n)
             (when (buffer-full-p)
               (deflate-to-bgz))
             n)
            (t
             (loop
                with buffer = (bgzf-buffer bgzf)
                for i from 0 below n
                do (progn
                     (when (buffer-full-p)
                       (setf buffer (deflate-to-bgz)))
                     (setf (aref buffer (bgzf-pointer bgzf)) (aref bytes i))
                     (incf (bgzf-pointer bgzf)))
                finally (return n)))))))

(defun bgzf-flush (bgzf &key (compress t) (append-eof t) (mtime 0))
  (if (plusp (bgzf-pointer bgzf))
      (let ((buffer (subseq (bgzf-buffer bgzf) 0 (bgzf-pointer bgzf)))
            (cdata (make-array +bgz-max-payload-length+ :element-type 'octet
                               :initial-element 0))
            (compression (if compress
                             (bgzf-compression bgzf)
                             0))
            (deflate-backoff 1024))
        (multiple-value-bind (deflated bytes-in bytes-out)
            (gz:deflate-vector buffer cdata
              :compression compression :suppress-header t
              :window-bits 15 :backoff deflate-backoff)
          (declare (ignore deflated))
          (let* ((udata (subseq buffer 0 bytes-in))
                 (bgz (make-bgz-member :mtime mtime
                                       :xlen +xlen+ :udata udata
                                       :cdata cdata :cend bytes-out
                                       :bsize (+ bytes-out
                                                 +member-header-length+
                                                 +member-footer-length+)
                                       :isize bytes-in
                                       :crc32 (gz:crc32 udata))))
            (write-bgz-member bgz (bgzf-stream bgzf) (bgzf-util-buffer bgzf))
            (if (< bytes-in (1- (length buffer))) ; Didn't fit
                (let ((overflow (subseq buffer bytes-in (bgzf-pointer bgzf))))
                  (replace buffer overflow)
                  (setf (bgzf-pointer bgzf) (length overflow))
                  (bgzf-flush bgzf)) ; Flush again to write the overflow
                (setf (bgzf-pointer bgzf) 0))))))
  (when append-eof
    (write-sequence *empty-bgz-record* (bgzf-stream bgzf))))

(defun write-bgz-member (bgz stream buffer)
  "Writes one BGZ member to STREAM.

Arguments:

- bgz (bgz member): A bgz member to write.
- stream (octet output-stream): An open stream.
- buffer (simple-octet-vector): A re-usable read buffer which must be
  able to contain at least 4 bytes.

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
      (encode-bytes (bgz-member-mtime bgz) 4 stream)         ; mtime
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
