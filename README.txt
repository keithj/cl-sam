Introduction

cl-sam is a Common Lisp toolkit for manipulation of DNA sequence
alignment data stored in the Sequence Alignment/Map (SAM) format
<http://samtools.sourceforge.net>. It is meant to be used in
conjunction with the SAMtools C toolkit.

The SAM specficiation describes text (SAM) and binary (BAM)
formats. cl-sam originated as a small set of functions for exploring
example BAM files and therefore provides only a subset of possible
SAM/BAM read/write operations.

The operations supported in this release are:

         SAM    BAM
Read      No    Yes
Write     No    No

BAM read performance is roughly 0.5x the speed of SAMtools C
implementation (SBCL 1.0.29.54 X86_64).


Documentation

See the Lisp docstrings, particularly the package docstrings for an
overview. HTML documentation may be generated with the command:

 (dxs:document-system :cl-sam)

at the REPL, provided that CLDOC is installed.

The components of cl-sam are divided by file as follows:

 bam.lisp          High-level BAM record data reading functions.
 bam-reader.lisp   High-level BAM file reading functions.

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

SAMtools                http://samtools.sourceforge.net

cl-sam requires a C shared library for block-gzip file reading. This
is created from the SAMtools C source using the makefile provided in
the bgzf directory of cl-sam.
