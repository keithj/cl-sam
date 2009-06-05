;;;
;;; Copyright (C) 2009 Keith James. All rights reserved.
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

(in-package :cl-user)

(eval-when (:compile-toplevel :load-toplevel :execute)
  (when (asdf:find-system :deoxybyte-systems nil)
    (asdf:operate 'asdf:load-op :deoxybyte-systems)))

(defpackage :cl-sam-system
  (:use :common-lisp :asdf :deoxybyte-systems))

(in-package :cl-sam-system)

(defsystem cl-sam
    :name "cl-sam"
    :author "Keith James"
    :licence "GPL v3"
    :depends-on (:cffi :deoxybyte-utilities :deoxybyte-io :deoxybyte-unix)
    :in-order-to ((test-op (load-op :cl-sam :cl-sam-test)))
    :components
    ((:module :cl-sam
              :serial t
              :pathname "src/"
              :components ((:file "package")
                           (:file "conditions")
                           (:file "bgzf-ffi")
                           (:file "bgzf-reader")
                           (:file "bgzf-stream")
                           (:file "bam-reader")
                           (:file "bam")))
     (:lift-test-config :lift-tests
                        :pathname "cl-sam-test.config"
                        :target-system :cl-sam)
     (:cldoc-config :cldoc-documentation
                    :pathname "doc/html/"
                    :target-system :cl-sam)))
