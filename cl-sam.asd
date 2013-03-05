;;;
;;; Copyright (c) 2009-2013 Keith James. All rights reserved.
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

(asdf:load-system :deoxybyte-systems)

(in-package :uk.co.deoxybyte-systems)

(defsystem cl-sam
    :name "cl-sam"
    :version "0.15.1"
    :author "Keith James"
    :licence "GPL v3"
    :depends-on ((:version :deoxybyte-gzip "0.5.1")
                 (:version :deoxybyte-unix "0.7.2"))
    :in-order-to ((test-op (load-op :cl-sam :cl-sam-test)))
    :components
    ((:module :cl-sam
              :serial t
              :pathname "src/"
              :components ((:file "package")
                           (:file "conditions")
                           (:file "bgzf")
                           (:file "bgzf-reader")
                           (:file "bgzf-writer")
                           (:file "bgzip-stream")
                           (:file "sam")
                           (:file "bam")
                           (:file "bam-index")
                           (:file "bam-index-reader")
                           (:file "bam-index-writer")
                           (:file "bam-indexer")
                           (:file "bam-reader")
                           (:file "bam-writer")
                           (:file "sam-writer")
                           (:file "external-bam-sort")
                           (:file "bam-utilities")
                           (:file "cl-sam")))
     (:lift-test-config :lift-tests
                        :pathname "cl-sam-test"
                        :target-system :cl-sam)
     (:cldoc-config :cldoc-documentation
                    :pathname "doc/html/"
                    :target-system :cl-sam)))
