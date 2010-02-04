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

(defvar *bgz-read-buffer* (make-array 4 :element-type 'octet :initial-element 0)
  "The buffer used by {defun read-bgz-member} for reading block gzip
data. Rebind this per thread to make {defun read-bgz-member}
reentrant.")

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
                  for bgz = (read-bgz-member stream)
                  until (or (null bgz)                      ; eof
                            (plusp (bgz-member-isize bgz))) ; skip if empty
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
                            (setf (bgzf-eof bgzf) t)))))
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
  (let ((len (if null-terminated
                 (1- n)
               n)))
    (let ((bytes (read-bytes bgzf n)))
      (make-sb-string bytes 0 (1- len)))))

(defun read-bgz-member (stream &optional (buffer *bgz-read-buffer*))
  "Reads one BGZ member from STREAM.

Arguments:

- stream (octet input-stream): An open stream.

Optional:

- buffer (simple-octet-vector): A re-usable read buffer which must be
  able to contain at least 4 bytes. Defaults to *BGZ-READ-BUFFER*
  . This argument is only required where the the function calls must
  be reentrant (normally achieved by rebinding *BGZ-READ-BUFFER* per
  thread).

Returns:

- A BGZ structure, or NIL."
  (declare (optimize (speed 3)))
  (flet ((decode-bytes (s n)
           (when (plusp (read-sequence buffer s :end n))
             (ecase n
               (1 (decode-uint8le buffer))
               (2 (decode-uint16le buffer))
               (4 (decode-uint32le buffer))))))
    (let ((fpos (file-position stream))
          (id1 (decode-bytes stream 1)))
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
          (unless (and (= gz:+id1+ id1)
                       (= gz:+id2+ id2)
                       (= gz:+cm-deflate+ cm)
                       (plusp (logand flg gz:+flag-extra+))
                       (= +xlen+ xlen)
                       (= +sf1+ sf1)
                       (= +sf2+ sf2)
                       (= +slen+ slen))
            (let ((tmpl #.(txt "invalid BGZ header id1:~d id2:~d cm:~d flg:~d"
                               "xlen:~d sf1:~d sf2:~d slen:~d at ~a in ~a")))
              (error 'bgzf-io-error
                     :text (format nil tmpl id1 id2 cm flg xlen sf1 sf2 slen
                                   fpos stream))))
          (let* ((deflated-size (1+ (logand bsize #xffff)))
                 (cdata-len (- deflated-size +member-header-length+
                               +member-footer-length+))
                 (cdata (make-array cdata-len :element-type 'octet
                                    :initial-element 0)))
            (read-sequence cdata stream)
            (let* ((crc32 (decode-bytes stream 4))
                   (isize (decode-bytes stream 4))) ; 8 footer bytes
              (make-bgz-member :mtime mtime :os os :xlen xlen
                               :bsize deflated-size
                               :cdata cdata :cend cdata-len
                               :crc32 crc32 :isize isize))))))))
