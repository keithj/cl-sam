Introduction

cl-sam is a Common Lisp toolkit for manipulation of DNA sequence
alignment data stored in the Sequence Alignment/Map (SAM) format
(http://samtools.sourceforge.net). cl-sam is de novo implementation of
the SAM spec in Common Lisp, using zlib via a C foreign function
interface.

While cl-sam is slower than the C and Java implementations for some
(but not all) operations, performance is good enough for real work.
Moreover, it offers the advantages of rapid development at all levels
in the SAM/BAM stack (zlib, bgzf, BAM and SAM) not afforded by
scripting lanuage bindings of samtools.

cl-sam is able to create BAM records de novo and may be used to create
a BAM file from scratch or edit a BAM stream. SAM/BAM header
manipulation functions are included so that new headers may be created
and headers from different BAM files may be systematically merged
without redundancy or avoidable conflicts. Where conflicts are
unavoidable, error conditions are raised to alert the user.

The operations supported in this version are:

         SAM    BAM
Read      No    Yes
Write    Yes    Yes

Sorting operations are available using an external merge sort that is
extensible by user-supplied sorting predicates. Typical coordinate
sort performance on the same 350 Mb queryname-sorted BAM file, median
of 5 runs:

                                 samtools C 0.1.8 sort      185 sec

                            Picard 1.32       SortSam
             (IcedTea Java 1.6 -server -Xmx2G, 64-bit)      170 sec

                            cl-sam 0.9.5 sort-bam-file
               (Steel Bank Common Lisp 1.0.43, 64-bit)      185 sec

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
 bam-index.lisp    Structures representing a samtools BAM index.

 sam.lisp          High-level SAM record data reading functions.
 sam-reader.lisp   High-level SAM file reading functions.
 sam-writer.lisp   High-level SAM file writing functions.

 bgzf-reader.lisp  Low-level BGZF seek and read functions.
 bgzf-stream.lisp  Low-level BGZF Gray stream implementation. An alternative
                   to the bgzf-reader functions, but roughly 2x slower.

The test suite contains examples of use.


Dependencies

deoxybyte-systems       git://github.com/keithj/deoxybyte-systems.git
deoxybyte-utilities     git://github.com/keithj/deoxybyte-utilities.git
deoxybyte-io            git://github.com/keithj/deoxybyte-io.git
deoxybyte-gzip          git://github.com/keithj/deoxybyte-gzip.git

CFFI                    http://common-lisp.net/project/cffi/


Optional dependencies

LIFT                    http://common-lisp.net/project/lift/
CLDOC                   http://common-lisp.net/project/cldoc/
