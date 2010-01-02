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
  (declare (optimize (speed 3) (safety 1)))
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
                  for bgz = (read-bgz-member stream)
                  until (or (null bgz)                      ; eof
                            (plusp (bgz-member-isize bgz))) ; skip if empty
                  finally (if bgz
                              (let ((udata (make-array (bgz-member-isize bgz)
                                                       :element-type 'octet))
                                    (pos (file-position stream)))
                                (check-type pos (unsigned-byte 48))
                                (inflate-vec (bgz-member-cdata bgz) udata)                                
                                (setf (bgzf-buffer bgzf) udata
                                      (bgzf-position bgzf) (1+ pos)
                                      (bgzf-offset bgzf) 0)
                                (return udata))
                            (setf (bgzf-eof bgzf) t)))))
      (unless (bgzf-init bgzf)
        (inflate-from-bgz)
        (setf (bgzf-init bgzf) t))
      (let ((inflated (or buffer (make-array n :element-type 'octet))))
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
