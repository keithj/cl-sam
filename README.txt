Introduction

cl-sam is a Common Lisp toolkit for manipulation of DNA sequence
alignment data stored in the Sequence Alignment/Map (SAM) format
<http://samtools.sourceforge.net>. It is meant to be used in
conjunction with the SAMtools C toolkit.

The SAM specficiation describes text (SAM) and binary (BAM)
formats. cl-sam originated as a small set of functions for exploring
example BAM files and therefore provides only a subset of possible
SAM/BAM read/write operations.

The operations supported in this version are:

         SAM    BAM
Read      No    Yes
Write     No    Yes

cl-sam is able to create BAM records de novo and may be used to create
a BAM file from scratch or edit a BAM stream. SAM/BAM header
manipulation functions are included so that new headers may be created
and headers from different BAM files may be systematically merged
without redundancy or avoidable conflicts. Where conflicts are
unavoidable, error conditions are raised to alert the user.

Sorting operations are available using an external merge sort that is
extensible by user-supplied sorting predicates.

BAM read performance is roughly 0.5x the speed of SAMtools C
implementation (SBCL 1.0.30 X86_64).


Installation

cl-sam uses ASDF for system definition. Copy or symlink cl-sam.asd
(and optionally cl-sam-test.asd) to your asdf:*central-registry* and
load cl-sam with the asdf:operate function:

 (asdf:operate 'asdf:load-op :cl-sam)

or with the equivalent deoxybyte-systems:load-system function:

 (dxs:load-system :cl-sam)


Tests

To run the unit and regression tests you need to have LIFT
installed. Run the tests with the asdf:operate function:

 (asdf:operate 'asdf:test-op :cl-sam)

or with the equivalent deoxybyte-systems:test-system function:

 (dxs:test-system :cl-sam)


Documentation

See the Lisp docstrings, particularly the package docstrings for an
overview. HTML documentation may be generated with the command:

 (dxs:document-system :cl-sam)

at the REPL, provided that CLDOC is installed.

The components of cl-sam are divided by file as follows:

 bam.lisp          High-level BAM record data reading functions.
 bam-reader.lisp   High-level BAM file reading functions.
 bam-writer.lisp   High-level BAM file writing functions.

 sam.lisp          High-level SAM record data reading functions.
 sam-reader.lisp   High-level SAM file reading functions.
 sam-writer.lisp   High-level SAM file writing functions.

 bgzf-reader.lisp  Low-level BGZF seek and read functions.
 bgzf-stream.lisp  Low-level BGZF Gray stream implementation. An alternative
                   to the bgzf-reader functions, but roughly 2x slower.

 bgzf-ffi.lisp     Low-level C FFI interface to BGZF files.

The test suite contains examples of use.


Dependencies

deoxybyte-systems       git://github.com/keithj/deoxybyte-systems.git
deoxybyte-utilities     git://github.com/keithj/deoxybyte-utilities.git
deoxybyte-io            git://github.com/keithj/deoxybyte-io.git
deoxybyte-unix          git://github.com/keithj/deoxybyte-unix.git

CFFI                    http://common-lisp.net/project/cffi/

SAMtools                http://samtools.sourceforge.net version >= 0.1.5

cl-sam requires a C shared library for block-gzip file reading. This
is created from the SAMtools C source using the makefile provided in
the bgzf directory of cl-sam.


Optional dependencies

LIFT                    http://common-lisp.net/project/lift/
CLDOC                   http://common-lisp.net/project/cldoc/
