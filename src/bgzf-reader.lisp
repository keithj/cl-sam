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

(defun read-bytes (bgzf n &key buffer)
  (declare (optimize (speed 3)))
  (let ((stream (bgzf-stream bgzf))
        (num-read 0))
    (declare (type vector-index n num-read))
    (labels ((buffer-empty-p ()
               (= (bgzf-offset bgzf) (1- (length (bgzf-buffer bgzf)))))
             (inflate-vec (source dest)
               (gz:inflate-vector source dest :suppress-header t
                                  :window-bits 15))
             (inflate-from-bgz ()
               (loop
                  for pos = (file-position stream)
                  for bgz = (read-bgz-member stream (bgzf-util-buffer bgzf))
                  until (or (null bgz)                      ; eof
                            (plusp (bgz-member-isize bgz))) ; not empty
                  finally (if bgz
                              (let ((udata (make-array (bgz-member-isize bgz)
                                                       :element-type 'octet
                                                       :initial-element 0)))
                                (check-type pos (unsigned-byte 48))
                                (inflate-vec (bgz-member-cdata bgz) udata)
                                (setf (bgzf-buffer bgzf) udata
                                      (bgzf-position bgzf) pos
                                      (bgzf-offset bgzf) (bgzf-load-seek bgzf)
                                      (bgzf-load-seek bgzf) 0)
                                (return udata))
                              (setf (bgzf-eof bgzf) t
                                    (bgzf-position bgzf) (file-position stream)
                                    (bgzf-offset bgzf) 0)))))
      (unless (bgzf-loaded-p bgzf)
        (inflate-from-bgz)
        (setf (bgzf-loaded-p bgzf) t))
      (let ((inflated (or buffer (make-array n :element-type 'octet
                                             :initial-element 0))))
        (declare (type simple-octet-vector inflated))
        (values
         (cond ((bgzf-eof bgzf)
                nil)
               ((< (+ (bgzf-offset bgzf) n) (1- (length (bgzf-buffer bgzf))))
                (replace inflated (bgzf-buffer bgzf) :start2 (bgzf-offset bgzf))
                (incf num-read n)
                (if (buffer-empty-p)
                    (inflate-from-bgz)
                    (incf (bgzf-offset bgzf) n))
                inflated)
               (t
                (loop
                   with buf = (bgzf-buffer bgzf)
                   for i = 0 then (1+ i)
                   while (and buf (< i n))
                   do (progn
                        (setf (aref inflated i) (aref buf (bgzf-offset bgzf)))
                        (incf num-read)
                        (if (buffer-empty-p)
                            (setf buf (inflate-from-bgz))
                            (incf (bgzf-offset bgzf))))
                   finally (return inflated))))
         num-read)))))

(defun read-string (bgzf n &key null-terminated)
  "Reads N characters from handle BGZF and returns them as a Lisp
string. The NULL-TERMINATED keyword is used to indicate whether the C
string is null-terminated so that the terminator may be consumed."
  (let ((bytes (read-bytes bgzf n)))
    (if null-terminated
        (octets-to-string bytes 0 (1- n))
        (octets-to-string bytes))))

(defun read-bgz-member (stream buffer)
  "Reads one BGZ member from STREAM, using BUFFER to hold integers as
they are read.

Arguments:

- stream (octet input-stream): An open stream.
- buffer (simple-octet-vector): A re-usable read buffer which must be
  able to contain at least 4 bytes.

Returns:

- A BGZ structure, or NIL."
  (declare (optimize (speed 3)))
  (flet ((decode-bytes (n)
           (when (plusp (read-sequence buffer stream :end n))
             (ecase n
               (1 (decode-uint8le buffer))
               (2 (decode-uint16le buffer))
               (4 (decode-uint32le buffer)))))
         (skip-string (stream)
           (loop
              for c of-type octet = (read-byte stream)
              until (zerop c))))
    (let ((fpos (file-position stream))
          (id1 (decode-bytes 1)))
      (when id1
        (let* ((id2 (decode-bytes 1))
               (cm (decode-bytes 1))
               (flg (decode-bytes 1))
               (mtime (decode-bytes 4))
               (xfl (decode-bytes 1))
               (os (decode-bytes 1))
               (xlen (decode-bytes 2))
               (sf1 (decode-bytes 1))
               (sf2 (decode-bytes 1))
               (slen (decode-bytes 2))
               (bsize (decode-bytes 2))) ; 18 header bytes
          (declare (ignore xfl))
          (unless (and (= gz:+id1+ id1)
                       (= gz:+id2+ id2)
                       (= gz:+cm-deflate+ cm)
                       (plusp (logand flg gz:+flag-extra+))
                       (= +xlen+ xlen)
                       (= +sf1+ sf1)
                       (= +sf2+ sf2)
                       (= +slen+ slen))
            (let ((tmpl #.(concatenate 'string
                                       "invalid BGZ header id1:~d id2:~d cm:~d "
                                       "flg:~d xlen:~d sf1:~d sf2:~d slen:~d "
                                       "aat ~a in ~a")))
              (error 'bgzf-io-error
                     :text (format nil tmpl id1 id2 cm flg xlen sf1 sf2 slen
                                   fpos stream))))
          (let* ((deflated-size (1+ (logand bsize #xffff)))
                 (cdata-len (- deflated-size +member-header-length+
                               +member-footer-length+))
                 (cdata (make-array cdata-len :element-type 'octet
                                    :initial-element 0)))
            ;; I haven't seen BGZF files in the wild having the
            ;; following, but potentially they may.
            ;; Skip fname if present
            (when (plusp (logand flg gz:+flag-name+))
              (skip-string stream))
            ;; skip fcomment if present
            (when (plusp (logand flg gz:+flag-comment+))
              (skip-string stream))
            ;; skip fhcrc if present
            (when (plusp (logand flg gz:+flag-fhcrc+))
              (dotimes (n 2)
                (read-byte stream)))
            (read-sequence cdata stream)
            (let* ((crc32 (decode-bytes 4))
                   (isize (decode-bytes 4))) ; 8 footer bytes
              (make-bgz-member :mtime mtime :os os :xlen xlen
                               :bsize deflated-size
                               :cdata cdata :cend cdata-len
                               :crc32 crc32 :isize isize))))))))
