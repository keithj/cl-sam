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

(defconstant +bgz-max-size+ (expt 2 16))

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
  (declare (optimize (speed 3)))
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
                       (plusp (logand flg gz:+flag-extra+))
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
                               :bsize deflated-size :cdata cdata
                               :crc32 crc32 :isize isize))))))))

(defun write-bgz-member (bgz stream &optional (buffer *bgz-write-buffer*))
  (declare (optimize (speed 3)))
  (flet ((encode-bytes (x n s)
           (ecase n
             (1 (encode-int8le x buffer))
             (2 (encode-int16le x buffer))
             (4 (encode-int32le x buffer)))
           (write-sequence buffer s :end n)))
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
    (encode-bytes (1- (bgz-member-bsize bgz)) 2 stream) ; 18 header bytes
    (write-sequence (bgz-member-cdata bgz) stream) ; Payload
    (encode-bytes (bgz-member-crc32 bgz) 4 stream)
    (encode-bytes (bgz-member-isize bgz) 4 stream) ; 8 footer bytes
    bgz))

(defstruct bgzf
  (stream nil :type t)
  (buffer (make-array +bgz-max-size+ :element-type '(unsigned-byte 8))
          :type (simple-array (unsigned-byte 8) (*)))
  (position 0 :type (unsigned-byte 48))
  (offset 0 :type uint16)
  (pointer 0 :type uint16)
  (init nil :type t)
  (eof nil :type t))

(defun bgzf-open (filespec &key (direction :input))
  (let ((stream (open filespec :element-type '(unsigned-byte 8)
                      :direction direction)))
    (make-bgzf :stream stream)))

(defun bgzf-close (bgzf)
  (let ((stream (bgzf-stream bgzf)))
    (when (output-stream-p stream)
       (bgzf-flush bgzf))
    (close stream)))

(defun bgzf-seek (bgzf position)
  (let ((stream (bgzf-stream bgzf))
        (pos (ash position -16))
        (off (logand position #xffff)))
    (when (file-position stream pos)
      (setf (bgzf-position bgzf) pos
            (bgzf-offset bgzf) off
            (bgzf-buffer bgzf) (bgz-member-udata (read-bgz-member stream))))))

(defun bgzf-tell (bgzf)
  (logior (ash (bgzf-position bgzf) 16) (bgzf-offset bgzf)))

;; The algorithm below relies on a fresh udata (uncompressed data)
;; vector being made for each bgz member as it is read. This vector is
;; setf'd to the buffer slot in the bgzf reader struct. The vector
;; length is used to determine where the uncompressed data end. There
;; may be scope for a speed increase by inflating directly into the
;; bgzf buffer. However, this would need the end index of the
;; uncompressed data to be tracked because the bgzf buffer would
;; always remain the same length. The current function is fast enough
;; (faster than using the samtools bgzf shared library via CFFI), so
;; I'm not going to do this yet.
(defun read-bytes (bgzf n)
  (declare (optimize (speed 3) (safety 1)))
  (declare (type array-index n))
  (let ((stream (bgzf-stream bgzf)))
    (flet ((buffer-empty-p ()
             (= (bgzf-offset bgzf) (1- (length (bgzf-buffer bgzf)))))
           (inflate-to-buffer ()
             (loop
                for bgz = (read-bgz-member stream)
                until (or (null bgz)                      ; eof
                          (plusp (bgz-member-isize bgz))) ; skip if empty
                finally (if bgz
                            (let ((udata (make-array (bgz-member-isize bgz)
                                                     :element-type
                                                     '(unsigned-byte 8)))
                                  (pos (file-position stream)))
                              (check-type pos (unsigned-byte 48))
                              (gz:inflate-vector (bgz-member-cdata bgz)
                                                 udata :suppress-header t
                                                 :window-bits 15)
                              (setf (bgzf-buffer bgzf) udata
                                    (bgzf-position bgzf) (1+ pos)
                                    (bgzf-offset bgzf) 0)
                              (return udata))
                          (setf (bgzf-eof bgzf) t)))))
      (unless (bgzf-init bgzf)
        (inflate-to-buffer)
        (setf (bgzf-init bgzf) t))
      (let ((bytes (make-array n :element-type '(unsigned-byte 8))))
        (cond ((bgzf-eof bgzf)
               nil)
              ((< (+ (bgzf-offset bgzf) n) (1- (length (bgzf-buffer bgzf))))
               (replace bytes (bgzf-buffer bgzf) :start2 (bgzf-offset bgzf))
               (if (buffer-empty-p)
                   (inflate-to-buffer)
                 (incf (bgzf-offset bgzf) n))
               bytes)
              (t
               (loop
                  with buffer = (bgzf-buffer bgzf)
                  for i = 0 then (1+ i)
                  while (< i n)
                  do (progn
                       (setf (aref bytes i) (aref buffer (bgzf-offset bgzf)))
                       (if (buffer-empty-p)
                           (setf buffer (inflate-to-buffer))
                         (incf (bgzf-offset bgzf))))
                  finally (return bytes))))))))

(defun write-bytes (bgzf bytes n)
  (let ((stream (bgzf-stream bgzf)))
    (labels ((buffer-full-p ()
               (= (bgzf-pointer bgzf) (1- (length (bgzf-buffer bgzf)))))
             (store-overflow (num-deflated)
               (let* ((buf (bgzf-buffer bgzf))
                      (overflow (subseq buf num-deflated (bgzf-pointer bgzf))))
                 (replace buf overflow)
                 (setf (bgzf-pointer bgzf) (length overflow))))
             (deflate-from-buffer ()
               (let ((buf (bgzf-buffer bgzf))
                     (cdata (make-array +bgz-max-size+ :element-type
                                        '(unsigned-byte 8))))
                 (multiple-value-bind (deflated bytes-in bytes-out)
                     (gz:deflate-vector buf cdata :compression 5
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
      (cond ((< (+ (bgzf-pointer bgzf) n)
                (1- (length (bgzf-buffer bgzf))))
             (replace (bgzf-buffer bgzf) bytes :start1 (bgzf-pointer bgzf))
             (incf (bgzf-pointer bgzf) n)
             (when (buffer-full-p)
               (deflate-from-buffer))
             n)
            (t
             (loop
                with buffer = (bgzf-buffer bgzf)
                for i = 0 then (1+ i)
                while (< i n)
                do (progn
                     (setf (aref buffer (bgzf-pointer bgzf)) (aref bytes i))
                     (if (buffer-full-p)
                         (setf buffer (deflate-from-buffer))
                       (incf (bgzf-pointer bgzf))))
                finally (return n)))))))

(defun bgzf-flush (bgzf)
  (if (plusp (bgzf-pointer bgzf))
      (let ((buffer (subseq (bgzf-buffer bgzf) 0 (bgzf-pointer bgzf)))
            (cdata (make-array +bgz-max-size+
                               :element-type '(unsigned-byte 8))))
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
  (write-sequence +empty-bgzf-record+ (bgzf-stream bgzf)))

(defun bgzf-virtual-position (bgzf &optional position)
  (if position
      (setf (bgzf-position bgzf) (ash position -16)
            (bgzf-offset bgzf) (logand position #xffff))
    (logior (ash (bgzf-position bgzf) 16) (bgzf-offset bgzf))))

;; As yet untested
(defun bgzf-eof-p (bgzf)
  (let* ((stream (bgzf-stream bgzf))
         (pos (file-position stream)))
    (unwind-protect
         (cond ((< (file-length stream) (length +empty-bgzf-record+))
                (error 'bgzf-io-error :text "incomplete file"))
               ((file-position stream (- (file-length stream)
                                         (length +empty-bgzf-record+)))
                (loop
                   for byte across +empty-bgzf-record+
                   always (equal byte (read-byte stream))))
               (t
                nil))
      (file-position stream pos))))
