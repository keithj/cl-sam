;;;
;;; Copyright (C) 2009-2011 Keith James. All rights reserved.
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

(defsystem cl-sam-test
    :depends-on (:cl-sam
                 (:version :lift "1.7.0")
                 (:version :deoxybyte-io "0.9.6"))
    :components ((:module :cl-sam-test
                          :serial t
                          :pathname "test/"
                          :components ((:file "package")
                                       (:file "pileup-parser")
                                       (:file "cl-sam-test")
                                       (:file "sam-test")
                                       (:file "bam-test")
                                       (:file "bam-index-test")))))
