
# This Makefile is, for the most part, a wrapper around the Lisp build
# system ASDF. It provides targets for building non-Lisp components,
# such as texinfo documentation.


.PHONY:	all core fasl doc clean

default: all

all: fasl doc

fasl:
	sbcl --dynamic-space-size 512 --noinform --noprint \
	--eval "(progn (asdf:operate 'asdf:compile-op :cl-sam) (quit))"

doc:
	sbcl --dynamic-space-size 512 --noinform --noprint --load make-doc.lisp

test:
	sbcl --dynamic-space-size 512 --noinform --noprint \
	--eval "(progn (asdf:operate 'asdf:test-op :cl-sam) (quit))"

coverage: clean
	sbcl --dynamic-space-size 512 --load make-coverage-report.lisp

clean:
	find . -name \*.fasl -exec rm {} \;
