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

(in-package :sam-test)

(deftestsuite cl-sam-tests ()
  ())

(defparameter *bgzf-eof-bytes*
  '(#o037 #o213 #o010 #o004 #o000 #o000 #o000 #o000 #o000 #o377
    #o006 #o000 #o102 #o103 #o002 #o000 #o033 #o000 #o003 #o000
    #o000 #o000 #o000 #o000 #o000 #o000 #o000 #o000))

(defun read-raw-data (filespec)
  (with-open-file (stream filespec :direction :input
                          :element-type '(unsigned-byte 8))
    (let ((data (make-array (file-length stream))))
      (read-sequence data stream)
      data)))

(defun test-binary-files (filespec1 filespec2)
  (with-open-file (s1 filespec1 :element-type '(unsigned-byte 8))
    (with-open-file (s2 filespec2 :element-type '(unsigned-byte 8))
      (ensure (loop
                 for byte1 = (read-byte s1 nil nil)
                 for byte2 = (read-byte s2 nil nil)
                 while (and byte1 byte2)
                 always (= byte1 byte2)
                 finally (cond ((null byte1)
                                (file-position s2 (1- (file-position s2))))
                               ((null byte2)
                                (file-position s1 (1- (file-position s1))))))
              :report "files were not byte identical")
      ;; Files may have an empty BGZF record as EOF terminator
      (let ((eof1 (loop
                     for byte = (read-byte s1 nil nil)
                     while byte
                     collect byte))
            (eof2 (loop
                     for byte = (read-byte s2 nil nil)
                     while byte
                     collect byte)))
        (ensure (or (null eof1) (equalp *bgzf-eof-bytes* eof1))
                :report "expected nil or ~a but found ~a"
                :arguments (*bgzf-eof-bytes*  eof1))
        (ensure (or (null eof2) (equalp *bgzf-eof-bytes* eof2))
                :report "expected nil or ~a but found ~a"
                :arguments (*bgzf-eof-bytes*  eof2))))))

