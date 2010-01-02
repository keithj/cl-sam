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

(declaim (ftype (function (bgzf bam-alignment fixnum &key (:compression fixnum))
                          fixnum) write-bytes))

(defun write-bytes (bgzf bytes n &key (compression 5))
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
                     (cdata (make-array +bgz-max-size+ :element-type 'octet)))
                 (multiple-value-bind (deflated bytes-in bytes-out)
                     (gz:deflate-vector buf cdata :compression compression
                                        :suppress-header t :window-bits 15
                                        :backoff 1024)
                   (declare (ignore deflated))
                   (let* ((udata (subseq buf 0 bytes-in))
                          (bgz (make-bgz-member
                                :mtime (get-universal-time) :xlen +xlen+
                                :udata buf :cdata (subseq cdata 0 bytes-out)
                                :bsize (+ bytes-out +member-header-length+
                                          +member-footer-length+)
                                :isize bytes-in :crc32 (gz:crc32 udata))))
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

(defun bgzf-flush (bgzf)
  (if (plusp (bgzf-pointer bgzf))
      (let ((buffer (subseq (bgzf-buffer bgzf) 0 (bgzf-pointer bgzf)))
            (cdata (make-array +bgz-max-size+ :element-type 'octet)))
        (multiple-value-bind (deflated bytes-in bytes-out)
            (gz:deflate-vector buffer cdata :compression 5 :suppress-header t
                               :window-bits 15 :backoff 1024)
          (declare (ignore deflated))
          (let* ((udata (subseq buffer 0 bytes-in))
                 (bgz (make-bgz-member :mtime (get-universal-time)
                                       :xlen +xlen+ :udata udata
                                       :cdata (subseq cdata 0 bytes-out)
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
  ;; (write-sequence *empty-bgzf-record* (bgzf-stream bgzf))
  )
