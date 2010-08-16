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

(in-package :sam-test)

(deftestsuite cl-sam-tests ()
  ())

(defun test-binary-files (filespec1 filespec2)
  (with-open-file (s1 filespec1 :element-type 'octet)
    (with-open-file (s2 filespec2 :element-type 'octet)
      (ensure (loop
                 for byte1 = (read-byte s1 nil nil)
                 for byte2 = (read-byte s2 nil nil)
                 while (and byte1 byte2)
                 always (= byte1 byte2)
                 finally (cond ((null byte1)
                                (file-position s2 (1- (file-position s2))))
                               ((null byte2)
                                (file-position s1 (1- (file-position s1))))))
              :report "files were not byte identical"))))
