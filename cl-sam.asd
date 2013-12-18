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
    :version "0.16.0"
    :author "Keith James"
    :licence "GPL v3"
    :depends-on ((:version :deoxybyte-systems "1.0.0")
                 (:version :deoxybyte-gzip "0.6.0")
                 (:version :deoxybyte-unix "0.8.0"))
    :in-order-to ((test-op (load-op :cl-sam :cl-sam-test))
                  (doc-op (load-op :cl-sam :cldoc)))
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
                           (:file "cl-sam"))))
    :perform (test-op :after (op c)
                      (maybe-run-lift-tests :cl-sam
                                            "cl-sam-test.config"))
    :perform (doc-op :after (op c)
                     (maybe-build-cldoc-docs :cl-sam "doc/html")))
