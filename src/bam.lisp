;;;
;;; Copyright (c) 2009-2013 Keith James. All rights reserved.
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

(in-package :sam)

(defconstant +tag-size+ 2
  "The size of a BAM auxilliary tag in bytes.")
(defconstant +alignment-size-tag+ 4
  "The size of the BAM alignment size indicator.")
(defconstant +null-byte+ #x00
  "The termination byte for BAM strings.")
(defconstant +unknown-quality+ #xff
  "The value for unavailable quality scores.")
(defconstant +unknown-position+ -1
  "The position value for unmapped reads.")
(defconstant +unknown-reference+ -1
  "The reference id for unmapped reads.")

(defvar *bam-magic* (make-array 4 :element-type 'octet
                                :initial-contents '(66 65 77 1))
  "The BAM file magic header bytes.")

(declaim (type simple-base-string *invalid-read-name-chars*))
(defvar *invalid-read-name-chars*
  (make-array 4 :element-type 'base-char
              :initial-contents '(#\Space #\Tab #\Linefeed #\Return)))

(defmacro with-bam ((var (&optional header num-refs ref-meta) filespec
                         &rest args &key (compress-header t)
                         (pad-header 0) index regions &allow-other-keys)
                    &body body)
  "Evaluates BODY with VAR bound to a new BAM generator function on
pathname designator FILESPEC. The direction (:input versus :output) is
determined by the stream-opening arguments in ARGS. The standard
generator interface functions NEXT and HAS-MORE-P may be used in
operations on the returned generator.

On reading, HEADER, NUM-REFS and REF-META will be automatically bound
to the BAM file metadata i.e. the metadata are read automatically and
the iterator positioned before the first alignment record.

Optionally, iteration may be restricted to a specific reference,
designated by REF-NUM and to alignments that start between reference
positions START and END. Furthermore, a BAM-INDEX object INDEX may be
provided to allow the stream to seek directly to the desired region of
the file.

On writing, HEADER, NUM-REFS and REF-META should be bound to
appropriate values for the BAM file metadata, which will be
automatically written to the underlying stream.

The COMPRESS-HEADER and PAD-HEADER keyword arguments are only
applicable on writing where they control whether the records and
header block should be compressed and whether the header string should
be padded with nulls to allow space for expansion. Writing an
uncompressed, padded header means that a header that fits in the first
BGZF block may be updated without re-writing the entire BAM file. The
BGZF COMPRESSION keyword is also accepted, permitting the Zlib
compression level to be set for the stream.

A list REGIONS may be supplied to limit the returned alignments to
specific references and reference coordinates. REGIONS may be region
objects or region designators in the form of list tuples

;;; (<reference designator> start end)

where a reference designator is either the reference name string or
its identifier number in the BAM file. REGIONS will be normalised
automatically by sorting according to reference position in the BAM
file (according to the BAM metadata) and then by start and
end. Overlapping regions will be merged.

For example:

To count all records in a BAM file:

;;; (with-bam (in () \"in.bam\")
;;;   (loop
;;;      while (has-more-p in)
;;;      count (next in)))

To copy only records with a mapping quality of >= 30 to another BAM
file:

;;; (with-bam (in (header n ref-meta) \"in.bam\")
;;;   (with-bam (out (header n rer-meta) \"out.bam\" :direction :output)
;;;     (let ((q30 (discarding-if (lambda (x)
;;;                                 (< (mapping-quality x) 30)) in)))
;;;       (loop
;;;          while (has-more-p q30)
;;;          do (consume out (next q30))))))

To count records with mapping quality of >= 30 in a set of genomic
ranges, using an index:

;;; (with-bam-index (index \"index.bai\")
;;;   (with-bam (bam () bam-file :index index
;;;                  :regions '((0 1000000 1100000) (1 1000000 1100000)))
;;;     (loop
;;;        while (has-more-p bam)
;;;        count (>= (mapping-quality (next bam)) 30))))"
  (with-gensyms (bgzf hdr nrefs rmeta)
    (let ((header (or header `,hdr))
          (num-refs (or num-refs `,nrefs))
          (ref-meta (or ref-meta `,rmeta)))
      (let ((default-header (with-output-to-string (s)
                              (write-sam-header
                               `((:HD (:VN . ,*sam-version*))) s))))
        `(with-bgzf (,bgzf ,filespec ,@(remove-key-values
                                        '(:compress-header :pad-header
                                          :index :regions) args))
           ,@(if (search '(:direction :output) args)
                 `((write-bam-meta ,bgzf
                                   ,(or header default-header)
                                   ,(or num-refs 0) ,ref-meta
                                   :compress ,compress-header
                                   :null-padding ,pad-header)
                   (let ((,var (make-bam-output ,bgzf)))
                     ,@body))
                 `((multiple-value-bind (,header ,num-refs ,ref-meta)
                       (read-bam-meta ,bgzf)
                     ,@(cond (ref-meta
                              `((declare (ignorable ,header ,num-refs))))
                             (num-refs
                              `((declare (ignorable ,header)))))
                     (let ((,var (cond ((and ,index ,regions)
                                        (make-bam-index-input
                                         ,bgzf ,index (normalise-regions
                                                       ,regions ,ref-meta)))
                                       (,regions
                                        (make-bam-scan-input
                                         ,bgzf (normalise-regions
                                                ,regions ,ref-meta)))
                                       (t
                                        (make-bam-full-scan-input ,bgzf)))))
                       ,@body)))))))))

(defmacro define-alignment-tag (tag value-type &optional docstring)
  "Defines a new alignment tag to hold a datum of a particular SAM
type.

Arguments:

- tag (symbol): The tag e.g. :rg
- value-type (symbol): The value type, one of :char , :string , :hex
:int32 or :float .

Optional:

- docstring (string): Documentation for the tag."
  (let ((encode-fn (ecase value-type
                     (:char 'encode-char-tag)
                     (:string 'encode-string-tag)
                     (:hex 'encode-hex-tag)
                     (:int32 'encode-int-tag)
                     (:float 'encode-float-tag)))
        (prefix (make-array 2 :element-type 'octet
                            :initial-contents (loop
                                                 for c across (symbol-name tag)
                                                 collect (char-code c)))))
    `(progn
       (defmethod encode-alignment-tag (value (tag (eql ,tag))
                                        aln index)
         (let ((i (+ +tag-size+ index)))
           (,encode-fn
            value (replace aln ,prefix :start1 index) i)))
       (defmethod alignment-tag-documentation ((tag (eql ,tag)))
         (declare (ignore tag))
         ,docstring))))

(defstruct (alignment (:conc-name aln-)
                      (:constructor make-aln (read-name seq-str qual-str flag
                                              ref-id mate-ref-id pos mate-pos
                                              cigar ref-len
                                              map-qual insert-length
                                              bin record))
                      (:print-object print-aln))
  (bin 0 :type uint16)
  (flag 0 :type uint16)
  (pos +unknown-position+ :type int32)
  (mate-pos +unknown-position+ :type int32)
  (ref-id +unknown-reference+ :type int32)
  (mate-ref-id +unknown-reference+ :type int32)
  (ref-len 0 :type fixnum)
  (map-qual 0 :type uint8)
  (insert-length 0 :type int32)
  (read-name "*" :type simple-string)
  (seq-str "*" :type simple-string)
  (qual-str "*" :type simple-string)
  (cigar nil :type list)
  (record (make-array 0 :element-type 'octet) :type simple-octet-vector))

(defun make-alignment (aln)
  "Returns a new ALIGNMENT structure given BAM record ALN."
  (make-aln (read-name aln) (seq-string aln) (quality-string aln)
            (alignment-flag aln) (reference-id aln) (mate-reference-id aln)
            (alignment-position aln) (mate-alignment-position aln)
            (alignment-cigar aln) (alignment-reference-length aln)
            (mapping-quality aln) (insert-length aln)
            (alignment-bin aln)
            aln))

(defun print-aln (aln stream)
  (print-unreadable-object (aln stream)
    (princ (aln-read-name aln) stream)
    (princ #\Space stream)
    (format stream "~d:~d..~d:~d ~d ~a ~a ~s ~s"
            (aln-ref-id aln) (aln-pos aln)
            (aln-mate-ref-id aln) (aln-mate-pos aln)
            (aln-map-qual aln)
            (flag-symbol (alignment-flag (aln-record aln)))
            (aln-cigar aln)
            (aln-seq-str aln)
            (aln-qual-str aln))))

(defun flag-symbol (flag)
  "Returns a 3 character string symbolising the read pairing denoted
by FLAG. The left character indicates the query strand, the middle the
alignment pairing and the right the mate strand.

Query mapped forward  >
Query mapped reverse  <
Query unmapped        .

Not paired             .
Mapped proper pair     =
Mapped pair            -
Singleton              ~

Mate mapped forward     >
Mate mapped reverse     <
Mate unmapped           ."
  (let ((str (with-output-to-string (s)
               (princ (cond ((query-unmapped-p flag)       #\.)
                            ((query-forward-p flag)        #\>)
                            (t                             #\<)) s)
               (princ (cond ((not (sequenced-pair-p flag)) #\.)
                            ((mapped-proper-pair-p flag)   #\=)
                            ((and (query-mapped-p flag)
                                  (mate-mapped-p flag))    #\-)
                            (t                             #\~)) s)
               (princ (cond ((not (sequenced-pair-p flag)) #\Space)
                            ((mate-unmapped-p flag)        #\.)
                            ((mate-forward-p flag)         #\<)
                            (t                             #\>)) s))))
    (if (first-in-pair-p flag)
        str
        (nreverse str))))

(defgeneric encode-alignment-tag (value tag vector index)
  (:documentation "Performs binary encoding of VALUE into VECTOR under
  TAG at INDEX, returning VECTOR."))

(defgeneric alignment-tag-documentation (tag)
  (:documentation "Returns the documentation for TAG or NIL if none is
  available."))

(defmethod encode-alignment-tag (value tag vector index)
  (declare (ignore value vector index))
  (error "Unknown tag ~a." tag))

(defmethod alignment-tag-documentation (tag)
  (error "Unknown tag ~a." tag))

(define-alignment-tag :am :int32
  "Smaller single-end mapping quality of the two reads in a pair.")
(define-alignment-tag :as :int32
  "Alignment score generated by aligner.")
(define-alignment-tag :bq :string
  (txt "Offset to base alignment quality (BAQ), of the same length as the"
       "read sequence. At the i-th read base, BAQi = Qi - (BQi - 64) where"
       "Qi is the i-th base quality."))
(define-alignment-tag :cc :string
  "Reference name of the next hit, \"=\" for the same chromosome.")
(define-alignment-tag :cm :int32
  "Number of color differences.")
(define-alignment-tag :cp :int32
  "Leftmost coordinate of the next hit.")
(define-alignment-tag :cq :string
  (txt "Color read quality on the same strand as the reference, encoded"
       "in the same way as the QUAL field."))
(define-alignment-tag :cs :string
  "Color read sequence on the same strand as the reference.")
(define-alignment-tag :e2 :string
  "The 2nd most likely base calls. Same encoding and same length as QUAL.")
(define-alignment-tag :fi :int32
  "The index of fragment in the template.")
(define-alignment-tag :fs :string
  "Fragment suffix.")
(define-alignment-tag :fz :string
  (txt "Flow signal intensities on the original strand of the read,"
       "stored as (uint16 t) round(value * 100.0)."))
(define-alignment-tag :gc :string
  (txt "CIGAR-like string describing the overlaps in the format of"
       "[0-9]+[SG]. For backward compatibility."))
(define-alignment-tag :gs :string
  "Sequence in the overlap. For backward compatibility.")
(define-alignment-tag :gq :string
  (txt "Quality in the overlap, encoded in the same way as the QUAL field."
       "For backward compatibility."))
(define-alignment-tag :lb :string
  (txt "Library. Value should be consistent with the header RG-LB tag if"
       "@RG is present."))
(define-alignment-tag :h0 :int32
  "Number of perfect hits.")
(define-alignment-tag :h1 :int32
  "Number of 1-difference hits (an in/del counted as a difference).")
(define-alignment-tag :h2 :int32
  "Number of 2-difference hits (an in/del counted as a difference).")
(define-alignment-tag :hi :int32
  (txt "Query hit index, indicating the alignment record is the i-th one"
       "stored in SAM."))
(define-alignment-tag :ih :int32
  (txt "Number of stored alignments in SAM that contains the query in the"
       "current record."))
(define-alignment-tag :md :string
  (txt "String for mismatching positions in the format of "
       "[0-9]+(([ACGTN]|\^[ACGTN]+)[0-9]+)*"))
(define-alignment-tag :mf :int32
  "MAQ pair flag (MAQ specific). For backward compatibility.")
(define-alignment-tag :mq :int32
  "The mapping quality score the mate alignment.")
(define-alignment-tag :nh :int32
  (txt "Number of reported alignments that contains the query in the"
       "current record."))
(define-alignment-tag :nm :int32
  "Number of nucleotide differences.")
(define-alignment-tag :pg :string
  (txt "Program that generates the alignment; match the header PG-ID tag"
       "if @PG is present."))
(define-alignment-tag :pq :int32
  (txt "Phred likelihood of the read pair, conditional on both the mapping"
       "locations being correct."))
(define-alignment-tag :pu :string
  (txt "Platform unit. Value should be consistent with the header RG-PU"
       "tag if @RG is present."))
(define-alignment-tag :q2 :string
  "Phred quality for the mate, encoded is the same as the QUAL field.")
(define-alignment-tag :r2 :string
  "Sequence of the mate.")
(define-alignment-tag :rg :string
  (txt "Read group. Value matches the header RG-ID tag if @RG is present"
       "in the header."))
(define-alignment-tag :s2 :string
  (txt "Encoded base probabilities for the other 3 bases for the"
       "mate-pair read, encoded in the same way as the SQ field."
       "For backward compatibility."))
(define-alignment-tag :sm :int32
  (txt "Mapping quality if the read is mapped as a single read rather"
       "than as a read pair."))
(define-alignment-tag :sq :string
  (txt "Encoded base probabilities for the suboptimal bases at each position."
       "For backward compatibility."))
(define-alignment-tag :u2 :int32
  (txt "Phred probility of the 2nd call being wrong conditional on"
       "the best being wrong. The same encoding as QUAL."))
(define-alignment-tag :uq :int32
  (txt "Phred likelihood of the read sequence, conditional on the mapping"
       "location being correct."))

;;; Various user tag extensions
;;; cl-sam
(define-alignment-tag :za :string "cl-sam: annotation")
(define-alignment-tag :zf :int32 "cl-sam: original flag before correction")

;;; BWA
(define-alignment-tag :x0 :int32 "BWA: number of best hits")
(define-alignment-tag :x1 :int32 "BWA: number of suboptimal hits")
(define-alignment-tag :xa :string "BWA: alternative hits (chr,pos,CIGAR,NM;)*")
(define-alignment-tag :xe :int32 "BWA: numnber of supporting reads")
(define-alignment-tag :xf :int32
  "BWA: support from forward or reverse alignment")
(define-alignment-tag :xg :int32 "BWA: number of gap extensions")
(define-alignment-tag :xm :int32 "BWA: number of mismatches in alignment")
(define-alignment-tag :xo :int32 "BWA: number of gap opens")
(define-alignment-tag :xs :int32 "BWA: suboptimal alignment score")
(define-alignment-tag :xt :char "BWA: Type: Unique/Repeat/N/Mate-sw")

(defun make-reference-table (ref-meta-list)
  "Returns a hash-table mapping reference identifiers to reference
names for the reference data in REF-META-LIST."
  (let ((ref-table (make-hash-table :size (length ref-meta-list))))
    (dolist (ref-meta ref-meta-list ref-table)
      (setf (gethash (first ref-meta) ref-table) (second ref-meta)))))

(defun make-alignment-record (read-name seq-str alignment-flag
                              &key (reference-id +unknown-reference+)
                              (alignment-pos +unknown-position+)
                              (mate-reference-id +unknown-reference+)
                              (mate-alignment-pos +unknown-position+)
                              (mapping-quality 0) (alignment-bin 0)
                              (insert-length 0)
                              cigar quality-str tag-values)
  "Returns a new alignment record array.

Arguments:

- read-name (string): The read name.
- seq-str (string): The read sequence.
- alignment-flag (integer): The binary alignment flag.

Key:

- reference-id (integer): The reference identifier, defaults to -1
- alignment-pos (integer): The 1-based alignment position, defaults to -1.
- mate-reference-id (integer): The reference identifier of the mate,
  defaults to -1.
- mate-alignment-pos (integer): The 1-based alignment position of the mate.
- mapping-quality (integer): The mapping quality, defaults to 0.
- alignment-bin (integer): The alignment bin, defaults to 0.
- insert-length (integer): The insert size, defaults to 0.
- cigar (alist): The cigar represented as an alist of operations e.g.

;;; '((:M . 9) (:I . 1) (:M . 25))

- quality-str (string): The read quality string.
- tag-values (alist): The alignment tags represented as an alist e.g.

;;; '((:XT . #\U) (:NM . 1) (:X0 . 1) (:X1 . 0)
;;;   (:XM . 1) (:XO . 0) (:XG . 0) (:MD . \"3T31\"))

The tags must have been defined with {defmacro define-alignment-tag} .

Returns:

- A vector of '(unsigned-byte 8)."
  (let ((slen (length seq-str))
        (qlen (length quality-str)))
    (check-arguments (or (null quality-str) (= slen qlen))
                     (seq-str quality-str)
                     (txt "read sequence and quality strings were not"
                          "the same length"))
    (let* ((i 32)
           (j (+ i (1+ (length read-name))))
           (k (+ j (* 4 (length cigar))))
           (m (+ k (ceiling slen 2)))
           (n (+ m (if (null quality-str)
                       slen
                       qlen)))
           (sizes (loop
                     for (nil . value) in tag-values
                     collect (alignment-tag-bytes value)))
           (aln (make-array (+ n (apply #'+ sizes)) :element-type 'octet
                            :initial-element 0)))
      (encode-int32le reference-id aln 0)
      (encode-int32le alignment-pos aln 4)
      (encode-int8le (1+ (length read-name)) aln 8)
      (encode-int8le mapping-quality aln 9)
      (encode-int16le alignment-bin aln 10)
      (encode-int16le (length cigar) aln 12)
      (encode-int16le alignment-flag aln 14)
      (encode-int16le slen aln 16)
      (encode-int32le mate-reference-id aln 20)
      (encode-int32le mate-alignment-pos aln 24)
      (encode-int32le insert-length aln 28)
      (encode-read-name read-name aln i)
      (encode-cigar cigar aln j)
      (encode-seq-string seq-str aln k)
      (if (null quality-str)
          (encode-unknown-quality slen aln m)
          (encode-quality-string quality-str aln m))
      (loop
         for (tag . value) in tag-values
         for size in sizes
         with offset = n
         do (progn
              (encode-alignment-tag value tag aln offset)
              (incf offset size))
         finally (return aln)))))

(defun copy-alignment-record (aln &key alignment-flag
                              (reference-id (reference-id aln))
                              (alignment-pos (alignment-position aln))
                              (mate-reference-id (mate-reference-id aln))
                              (mate-alignment-pos (mate-alignment-position aln))
                              (mapping-quality (mapping-quality aln))
                              (alignment-bin (alignment-bin aln))
                              (insert-length (insert-length aln))
                              (cigar (alignment-cigar aln))
                              (tag-values (alignment-tag-values aln)))
  "Returns a copy of ALN, optionally setting some fields to new values."
  (let ((alignment-flag (or alignment-flag (alignment-flag aln))))
    (make-alignment-record (read-name aln) (seq-string aln) alignment-flag
                           :reference-id reference-id
                           :alignment-pos alignment-pos
                           :mate-reference-id mate-reference-id
                           :mate-alignment-pos mate-alignment-pos
                           :mapping-quality mapping-quality
                           :alignment-bin alignment-bin
                           :insert-length insert-length :cigar cigar
                           :quality-str (quality-string aln)
                           :tag-values tag-values)))

(defun flag-bits (flag &rest bit-names)
  "Returns an integer FLAG that had BAM flag bits named by symbols
BIT-NAMES set.

Arguments:

-  flag (unsigned-byte 8): a BAM alignment flag.

Rest:

- bit-names (symbols): Any number of valid bit flag names:

;;; :sequenced-pair
;;; :mapped-proper-pair
;;; :query-mapped , :query-unmapped
;;; :mate-mapped , :mate-unmapped
;;; :query-forward , :query-reverse
;;; :mate-forward , :mate-reverse
;;; :first-in-pair , :second-in-pair
;;; :alignment-primary , :alignment-not-primary
;;; :fails-platform-qc
;;; :pcr/optical-duplicate

Returns:

- An (unsigned-byte 8)"
  (let ((f flag))
    (dolist (name bit-names (ensure-valid-flag f))
      (destructuring-bind (bit value)
          (ecase name
            (:not-sequenced-pair      '( 0 0))
            (:sequenced-pair          '( 0 1))
            (:not-mapped-proper-pair  '( 1 0))
            (:mapped-proper-pair      '( 1 1))
            (:query-mapped            '( 2 0))
            (:query-unmapped          '( 2 1))
            (:mate-mapped             '( 3 0))
            (:mate-unmapped           '( 3 1))
            (:query-forward           '( 4 0))
            (:query-reverse           '( 4 1))
            (:mate-forward            '( 5 0))
            (:mate-reverse            '( 5 1))
            (:first-in-pair           '( 6 1))
            (:second-in-pair          '( 7 1))
            (:alignment-primary       '( 8 0))
            (:alignment-not-primary   '( 8 1))
            (:fails-platform-qc       '( 9 1))
            (:pcr/optical-duplicate   '(10 1)))
        (setf (ldb (byte 1 bit) f) value)))))

(declaim (ftype (function (simple-octet-vector) (signed-byte 32))
                reference-id))
(defun reference-id (aln)
  "Returns the reference sequence identifier of ALN. This is an
integer locally assigned to a reference sequence within the context of
a BAM file."
  (declare (optimize (speed 3)))
  (decode-int32le aln 0))

(declaim (ftype (function (simple-octet-vector) (signed-byte 32))
                alignment-position))
(defun alignment-position (aln)
  "Returns the 0-based sequence coordinate of ALN in the reference
sequence of the first base of the clipped read."
  (declare (optimize (speed 3)))
  (decode-int32le aln 4))

(defun alignment-read-length (aln)
  "Returns the length of the alignment on the read."
  (loop
     for (op . len) in (alignment-cigar aln)
     when (member op '(:i :m :s := :x)) sum len))

(defun alignment-reference-length (aln)
  "Returns the length of the alignment on the reference."
  (loop
     for (op . len) in (alignment-cigar aln)
     when (member op '(:d :m :n := :x)) sum len))

(defun alignment-reference-end (aln)
  "Returns the end position of the alignment on the
reference. Requires decoding the CIGAR string."
  (1- (+ (alignment-position aln) (alignment-reference-length aln))))

(declaim (inline read-name-length))
(defun read-name-length (aln)
  "Returns the length in ASCII characters of the read name of ALN."
  (decode-uint8le aln 8))

(defun mapping-quality (aln)
  "Returns the integer mapping quality of ALN."
  (decode-uint8le aln 9))

(defun alignment-bin (aln)
  "Returns an integer that indicates the alignment bin to which
ALN has been assigned."
  (decode-uint16le aln 10))

(defun cigar-length (aln)
  "Returns the number of CIGAR operations in ALN."
  (decode-uint16le aln 12))

(declaim (ftype (function ((simple-octet-vector) &key (:validate t))
                          (unsigned-byte 16))
                alignment-flag))
(defun alignment-flag (aln &key validate)
  "Returns an integer whose bits are flags that describe properties of
the ALN. If the VALIDATE key is T (the default) the flag's bits are
checked for internal consistency."
  (let ((flag (decode-uint16le aln 14)))
    (when validate
      (ensure-valid-flag flag aln))
    flag))

(defun sequenced-pair-p (flag)
  "Returns T if FLAG indicates that the read was sequenced as a member
of a pair, or NIL otherwise."
  (logbitp 0 flag))

(defun mapped-proper-pair-p (flag)
  "Returns T if FLAG indicates that the read was mapped as a member of
a properly oriented read-pair, or NIL otherwise."
  (logbitp 1 flag))

(defun query-unmapped-p (flag)
  "Returns T if FLAG indicates that the read was not mapped to a
reference, or NIL otherwise."
  (logbitp 2 flag))

(defun query-mapped-p (flag)
  "Returns T if FLAG indicates that the read's mate was mapped to a
reference, or NIL otherwise."
  (not (query-unmapped-p flag)))

(defun mate-unmapped-p (flag)
  "Returns T if FLAG indicates that the read's mate was not mapped to
a reference, or NIL otherwise."
  (logbitp 3 flag))

(defun mate-mapped-p (flag)
  "Returns T if FLAG indicates that the read's mate was mapped to a
reference, or NIL otherwise."
  (not (mate-unmapped-p flag)))

(declaim (inline query-forward-p))
(defun query-forward-p (flag)
  "Returns T if FLAG indicates that the read was mapped to the forward
strand of a reference, or NIL if it was mapped to the reverse strand."
  (not (query-reverse-p flag)))

(declaim (inline query-reverse-p))
(defun query-reverse-p (flag)
  "Returns T if FLAG indicates that the read was mapped to the reverse
strand of a reference, or NIL if it was mapped to the forward strand."
  (declare (type uint16 flag))
  (logbitp 4 flag))

(defun mate-forward-p (flag)
  "Returns T if FLAG indicates that the read's mate was mapped to the
forward, or NIL if it was mapped to the reverse strand."
  (not (mate-reverse-p flag)))

(defun mate-reverse-p (flag)
  "Returns T if FLAG indicates that the read's mate was mapped to the
reverse, or NIL if it was mapped to the forward strand."
  (logbitp 5 flag))

(defun first-in-pair-p (flag)
  "Returns T if FLAG indicates that the read was the first in a pair
of reads from one template, or NIL otherwise."
  (logbitp 6 flag))

(defun second-in-pair-p (flag)
  "Returns T if FLAG indicates that the read was the second in a pair
of reads from one template, or NIL otherwise."
  (logbitp 7 flag))

(defun alignment-not-primary-p (flag)
  "Returns T if FLAG indicates that the read mapping was not the
primary mapping to a reference, or NIL otherwise."
  (logbitp 8 flag))

(defun alignment-primary-p (flag)
  "Returns T if FLAG indicates that the read mapping was the primary
mapping to a reference, or NIL otherwise."
  (not (alignment-not-primary-p flag)))

(defun fails-platform-qc-p (flag)
  "Returns T if FLAG indicates that the read failed plaform quality
control, or NIL otherwise."
  (logbitp 9 flag))

(defun pcr/optical-duplicate-p (flag)
  "Returns T if FLAG indicates that the read is a PCR or optical
duplicate, or NIL otherwise."
  (logbitp 10 flag))

(setf (fdefinition 'multiple-frags-p) #'sequenced-pair-p)
(setf (fdefinition 'proper-aligned-frags-p) #'mapped-proper-pair-p)
(setf (fdefinition 'frag-unmapped-p) #'query-unmapped-p)
(setf (fdefinition 'frag-mapped-p) #'query-mapped-p)
(setf (fdefinition 'next-unmapped-p) #'mate-unmapped-p)
(setf (fdefinition 'next-mapped-p) #'mate-mapped-p)
(setf (fdefinition 'frag-forward-p) #'query-forward-p)
(setf (fdefinition 'frag-reverse-p) #'query-reverse-p)
(setf (fdefinition 'next-forward-p) #'mate-forward-p)
(setf (fdefinition 'next-reverse-p) #'mate-reverse-p)
(setf (fdefinition 'first-frag-p) #'first-in-pair-p)
(setf (fdefinition 'last-frag-p) #'second-in-pair-p)

(defun medial-frag-p (flag)
  "Returns T if FLAG indicates that the read was sequenced as a
non-terminal part of a linear template, or NIL otherwise."
  (and (logbitp 6 flag) (logbitp 7 flag)))

(defun unknown-frag-p (flag)
  "Returns T if FLAG indicates that the read was sequenced as part of
a non-linear template or as part of a linear template where the
position has been lost, or NIL otherwise."
  (not (or (logbitp 6 flag) (logbitp 7 flag))))

(defun valid-mapped-pair-p (flag)
  "Returns T if FLAG indicates valid mapping states for a pair of
mapped reads, that is both must be mapped, or NIL otherwise."
  (and (query-mapped-p flag) (mate-mapped-p flag)))

(defun valid-mapped-proper-pair-p (flag)
  "Returns T if FLAG indicates valid proper mapping states for a pair
of mapped reads, that is both must be mapped and on opposite strands,
or NIL otherwise."
  (and (valid-mapped-pair-p flag)
       (not (eql (query-forward-p flag) (mate-forward-p flag)))))

(defun valid-flag-p (flag)
  "Returns T if the paired-read-specific bits of FLAG are internally
consistent."
  (cond ((and (sequenced-pair-p flag) (mapped-proper-pair-p flag))
         (valid-mapped-proper-pair-p flag))
        ((not (sequenced-pair-p flag))
         (not (or (mate-unmapped-p flag)
                  (first-in-pair-p flag)
                  (second-in-pair-p flag))))
        (t
         t)))

(defun valid-read-name-p (str)
  "Returns T if STR is a valid read name matching the regex
[!-?A-~]1,255 , or NIL otherwise. The length limit is implicit in the
BAM format, there being 8 bits to store the read name length. This
function is not used as it is not clear that this check is accepted by
other implementations."
  (flet ((name-char-p (c)
           (and (char/= #\@ c) (< 32 (char-code c) 127))))
    (let ((len (length str)))
      (and (< 0 len 256) (loop
                            for c across str
                            always (name-char-p c))))))

(defun read-length (aln)
  "Returns the length of the read described by ALN."
  (decode-int32le aln 16))

(defun mate-reference-id (aln)
  "Returns the integer reference ID of ALN."
  (declare (optimize (speed 3)))
  (decode-int32le aln 20))

(defun mate-alignment-position (aln)
  "Returns the 0-based sequence position of the read mate's alignment
described by ALN."
  (decode-int32le aln 24))

(defun insert-length (aln)
  "Returns the insert length described by ALN."
  (decode-int32le aln 28))

(defun read-name (aln)
  "Returns the read name string described by ALN."
  (decode-read-name aln 32 (read-name-length aln)))

(defun alignment-cigar (aln)
  "Returns the CIGAR record list of the alignment described by
ALN. CIGAR operations are given as a list, each member being a list of
a CIGAR operation keyword and an integer operation length."
  (let* ((name-len (read-name-length aln))
         (cigar-index (+ 32 name-len))
         (cigar-bytes (* 4 (cigar-length aln))))
    (decode-cigar aln cigar-index cigar-bytes)))

(defun seq-string (aln)
  "Returns the sequence string described by ALN."
  (multiple-value-bind (read-len cigar-index cigar-bytes seq-index
                        qual-index tag-index)
      (alignment-indices aln)
    (declare (ignore cigar-index cigar-bytes qual-index tag-index))
    (decode-seq-string aln seq-index read-len)))

(defun quality-string (aln)
  "Returns the sequence quality string described by ALN."
  (multiple-value-bind (read-len cigar-index cigar-bytes seq-index
                        qual-index tag-index)
      (alignment-indices aln)
    (declare (ignore cigar-index cigar-bytes seq-index tag-index))
    (decode-quality-string aln qual-index read-len)))

(defun alignment-tag-values (aln)
  "Returns an alist of tag and values described by ALN."
  (multiple-value-bind (read-len cigar-index cigar-bytes seq-index
                        qual-index tag-index)
      (alignment-indices aln)
    (declare (ignore read-len cigar-index cigar-bytes qual-index seq-index))
    (decode-tag-values aln tag-index)))

(defun alignment-core (aln &key validate)
  "Returns a list of the core data described by ALN. The list elements
are comprised of reference-id, alignment-position, read-name length,
mapping-quality alignment-bin, cigar length, alignment flag, read
length, mate reference-id, mate alignment-position and insert length."
  (list (reference-id aln)
        (alignment-position aln)
        (read-name-length aln)
        (mapping-quality aln)
        (alignment-bin aln)
        (cigar-length aln)
        (alignment-flag aln :validate validate)
        (read-length aln)
        (mate-reference-id aln)
        (mate-alignment-position aln)
        (insert-length aln)))

(defun alignment-core-alist (aln &key validate)
  "Returns the same data as {defun alignment-core} in the form of an
alist."
  (pairlis '(:reference-id :alignment-pos :read-name-length
             :mapping-quality :alignment-bin :cigar-length
             :alignment-flag :read-length :mate-reference-id
             :mate-alignment-position :insert-length)
           (alignment-core aln :validate validate)))

(defun alignment-flag-alist (aln &key validate)
  "Returns the bitwise flags of ALN in the form of an alist. The
primary purpose of this function is debugging."
  (let ((flag (alignment-flag aln :validate validate)))
    (pairlis '(:sequenced-pair :mapped-proper-pair :query-unmapped
               :mate-unmapped :query-forward :mate-forward :first-in-pair
               :second-in-pair :alignment-not-primary :fails-platform-qc
               :pcr/optical-duplicate)
             (list (sequenced-pair-p flag)
                   (mapped-proper-pair-p flag)
                   (query-unmapped-p flag)
                   (mate-unmapped-p flag)
                   (query-forward-p flag)
                   (mate-forward-p flag)
                   (first-in-pair-p flag)
                   (second-in-pair-p flag)
                   (alignment-not-primary-p flag)
                   (fails-platform-qc-p flag)
                   (pcr/optical-duplicate-p flag)))))

(defun alignment-indices (aln)
  "Returns 7 integer values which are byte-offsets within ALN at which
the various core data lie. See the SAM spec."
  (let* ((read-len (read-length aln))
         (name-len (read-name-length aln))
         (cigar-index (+ 32 name-len))
         (cigar-bytes (* 4 (cigar-length aln)))
         (seq-index (+ cigar-index cigar-bytes))
         (seq-bytes (ceiling read-len 2))
         (qual-index (+ seq-index seq-bytes))
         (tag-index (+ qual-index read-len)))
    (values read-len cigar-index cigar-bytes seq-index qual-index tag-index)))

(defun decode-read-name (aln index num-bytes)
  "Returns a string containing the template/read name of length
NUM-BYTES, encoded at byte INDEX in ALN."
  ;; The read name is null terminated and the terminator is included
  ;; in the name length
  (octets-to-string aln index (1- (+ index num-bytes))))

(defun encode-read-name (read-name aln index)
  "Returns ALN having encoded READ-NAME into it, starting at INDEX."
  (loop
     for i from 0 below (length read-name)
     for j = index then (1+ j)
     do (setf (aref aln j) (char-code (char read-name i)))
     finally (setf (aref aln (1+ j)) +null-byte+))
  aln)

(defun decode-seq-string (aln index num-bytes)
  "Returns a string containing the alignment query sequence of length
NUM-BYTES. The sequence must be present in ALN at INDEX."
  (declare (optimize (speed 3)))
  (declare (type simple-octet-vector aln)
           (type (unsigned-byte 32) index num-bytes))
  (flet ((decode-base (nibble)
           (ecase nibble
             (0 #\=)
             (1 #\A)
             (2 #\C)
             (4 #\G)
             (8 #\T)
             (15 #\N))))
    (loop
       with seq = (make-array num-bytes :element-type 'base-char
                              :initial-element #\Nul)
       for i from 0 below num-bytes
       for j of-type uint32 = (+ index (floor i 2))
       do (setf (char seq i)
                (decode-base (if (evenp i)
                                 (ldb (byte 4 4) (aref aln j))
                                 (ldb (byte 4 0) (aref aln j)))))
       finally (return seq))))

(defun encode-seq-string (str aln index)
  "Returns ALN having encoded STR into it, starting at INDEX."
  (flet ((encode-base (char)
           (ecase (char-upcase char)
             (#\= 0)
             (#\A 1)
             (#\C 2)
             (#\G 4)
             (#\T 8)
             (#\N 15))))
    (loop
       for i from 0 below (length str)
       for j = (+ index (floor i 2))
       for nibble = (encode-base (char str i))
       do (if (evenp i)
              (setf (ldb (byte 4 4) (aref aln j)) nibble)
              (setf (ldb (byte 4 0) (aref aln j)) nibble))
       finally (return aln))))

(defun decode-quality-string (aln index num-bytes)
  "Returns a string containing the alignment query sequence of length
NUM-BYTES. The sequence must be present in ALN at INDEX. The SAM spec
states that quality data are optional, with absence indicated by
0xff. If the first byte of quality data is 0xff, NIL is returned."
  (flet ((encode-phred (x)
           (code-char (+ 33 (min 93 x)))))
    (if (= +unknown-quality+ (aref aln index))
        nil
        (loop
           with str = (make-array num-bytes :element-type 'base-char
                                  :initial-element #\Nul)
           for i from 0 below num-bytes
           for j from index below (+ index num-bytes)
           do (setf (char str i)
                    (encode-phred (decode-uint8le aln j)))
           finally (return str)))))

(defun encode-quality-string (str aln index)
  "Returns ALN having encoded quality string STR into it, starting at
INDEX."
  (flet ((decode-phred (x)
           (- (char-code x) 33)))
    (loop
       for i from 0 below (length str)
       for j = (+ index i)
       do (encode-int8le (decode-phred (char str i)) aln j)))
  aln)

(defun encode-unknown-quality (n aln index)
  "Returns ALN having encoded N unknown qualities into it, starting at
INDEX."
  (loop
     for i from index below (+ index n)
     do (setf (aref aln i) +unknown-quality+))
  aln)

(defun decode-cigar (aln index num-bytes)
  "Returns an alist of CIGAR operations from NUM-BYTES bytes within
ALN, starting at INDEX. The decoding of the = and X operations is not
documented in the spec."
  (flet ((decode-len (uint32)
           (ash uint32 -4))
         (decode-op (uint32)
           (ecase (ldb (byte 4 0) uint32)
             (0 :m)
             (1 :i)
             (2 :d)
             (3 :n)
             (4 :s)
             (5 :h)
             (6 :p)
             (7 :=)
             (8 :x))))
    (loop
       for i from index below (1- (+ index num-bytes)) by 4
       collect (let ((x (decode-uint32le aln i)))
                 (cons (decode-op x) (decode-len x))))))

(defun encode-cigar (cigar aln index)
  "Returns ALN having encoded alist CIGAR into it, starting at INDEX."
  (flet ((encode-op-len (op len)
           (let ((uint32 (ash len 4)))
             (setf (ldb (byte 4 0) uint32) (ecase op
                                             (:m 0)
                                             (:i 1)
                                             (:d 2)
                                             (:n 3)
                                             (:s 4)
                                             (:h 5)
                                             (:p 6)
                                             (:= 7)
                                             (:x 8)))
             uint32)))
    (loop
       for (op . length) in cigar
       for i = index then (+ 4 i)
       do (encode-int32le (encode-op-len op length) aln i)
       finally (return aln))))

(defun decode-tag-values (aln index)
  "Returns a list of auxilliary data from ALN at INDEX. The BAM
two-letter data keys are transformed to Lisp keywords."
  (declare (optimize (speed 3) (safety 0)))
  (declare (type simple-octet-vector aln)
           (type vector-index index))
  (loop
     with tag-index of-type vector-index = index
     while (< tag-index (length aln))
     collect (let* ((type-index (+ tag-index +tag-size+))
                    (type-code (code-char (aref aln type-index)))
                    (tag (intern (octets-to-string aln tag-index
                                                   (+ 2 tag-index)) 'keyword))
                    (val-index (1+ type-index)))
               (declare (type vector-index val-index))
               (let  ((val (ecase type-code
                             (#\A         ; A printable character
                              (setf tag-index (+ val-index 1))
                              (code-char (aref aln val-index)))
                             (#\C         ; C unsigned 8-bit integer
                              (setf tag-index (+ val-index 1))
                              (decode-uint8le aln val-index))
                             ((#\H #\Z) ; H hex string, Z printable string
                              (let ((end (position +null-byte+ aln
                                                   :start val-index)))
                                (setf tag-index (1+ end))
                                (octets-to-string aln val-index end)))
                             (#\I         ; I unsigned 32-bit integer
                              (setf tag-index (+ val-index 4))
                              (decode-uint32le aln val-index))
                             (#\S         ; S unsigned short
                              (setf tag-index (+ val-index 2))
                              (decode-uint16le aln val-index))
                             (#\c         ; c signed 8-bit integer
                              (setf tag-index (+ val-index 1))
                              (decode-int8le aln val-index))
                             (#\f         ; f single-precision float
                              (setf tag-index (+ val-index 4))
                              (decode-float32le aln val-index))
                             (#\i         ; i signed 32-bit integer
                              (setf tag-index (+ val-index 4))
                              (decode-int32le aln val-index))
                             (#\s         ; s signed short
                              (setf tag-index (+ val-index 2))
                              (decode-int16le aln val-index)))))
                 (cons tag val)))))

(defun encode-int-tag (value aln index)
  "Returns ALN having encoded integer VALUE into it, starting at
INDEX. BAM format is permitted to use more compact integer storage
where possible."
  (destructuring-bind (type-char encoder)
      (etypecase value
        ((integer 0 255)                  '(#\C encode-int8le))
        ((integer -128 127)               '(#\c encode-int8le))
        ((integer 0 65535)                '(#\S encode-int16le))
        ((integer -32768 32767)           '(#\s encode-int16le))
        ((integer -2147483648 2147483647) '(#\I encode-int32le))
        ((integer 0 4294967295)           '(#\i encode-int32le)))
    (setf (aref aln index) (char-code type-char))
    (funcall encoder value aln (1+ index))))

(defun encode-float-tag (value aln index)
  "Returns ALN having encoded float VALUE into it,
starting at INDEX."
  (setf (aref aln index) (char-code #\f))
  (encode-float32le value aln (1+ index)))

(defun encode-char-tag (value aln index)
  "Returns ALN having encoded character VALUE into it,
starting at INDEX."
  (setf (aref aln index) (char-code #\A))
  (encode-int8le (char-code value) aln (1+ index)))

(defun encode-hex-tag (value aln index)
  "Returns ALN having encoded hex string VALUE into it,
starting at INDEX."
  (when (parse-integer value :radix 16)
    (setf (aref aln index) (char-code #\H))
    (%encode-string-tag value aln (1+ index)) ))

(defun encode-string-tag (value aln index)
  "Returns ALN having encoded string VALUE into it,
starting at INDEX."
  (setf (aref aln index) (char-code #\Z))
  (%encode-string-tag value aln (1+ index)))

(declaim (inline %encode-string-tag))
(defun %encode-string-tag (value aln index)
  (let* ((len (length value))
         (term-index (+ index len)))
    (loop
       for i from 0 below len
       for j = index then (1+ j)
       do (setf (aref aln j) (char-code (char value i)))
       finally (setf (aref aln term-index) +null-byte+))
    aln))

(defun alignment-tag-bytes (value)
  "Returns the number of bytes required to encode VALUE."
  (etypecase value
    (character 4)
    (string (+ 4 (length value)))       ; includes null byte
    (single-float 7)
    ((integer 0 255) 4)
    ((integer -128 127) 4)
    ((integer 0 65535) 5)
    ((integer -32768 32767) 5)
    ((integer -2147483648 2147483647) 7)
    ((integer 0 4294967295) 7)))

(defun ensure-valid-reference-name (str)
  "Returns STR if it is a valid reference sequence name, or raises a
{define-condition malformed-field-error} if not."
  (declare (optimize (speed 3)))
  (and (check-field (valid-reference-name-p str) nil str
                    "invalid reference name")
       str))

(defun ensure-valid-read-name (str)
  (declare (optimize (speed 3)))
  (declare (type simple-string str))
  (loop
     for c across str
     do (check-field (not (find c *invalid-read-name-chars* :test #'char=))
                     nil str "invalid character ~c in read name" c)
     finally (return str)))

(defun ensure-valid-flag (flag &optional aln)
  (cond ((mapped-proper-pair-p flag)
         (cond ((not (sequenced-pair-p flag))
                (flag-validation-error
                 flag (txt "the sequenced-pair bit was not set in a mapped"
                           "proper pair flag") aln))
               ((not (valid-mapped-pair-p flag))
                (flag-validation-error
                 flag (txt "one read was marked as unmapped in a mapped"
                           "proper pair flag") aln))
               ((not (valid-mapped-proper-pair-p flag))
                (flag-validation-error
                 flag (txt "reads were not mapped to opposite strands in a"
                           "mapped proper pair flag") aln))
               (t
                flag)))
        ((not (sequenced-pair-p flag))
         (cond ((mate-reverse-p flag)
                (flag-validation-error
                 flag "the mate-reverse bit was set in an unpaired read"
                 aln))
               ((mate-unmapped-p flag)
                (flag-validation-error
                 flag "the mate-unmapped bit was set in an unpaired read"
                 aln))
               ((first-in-pair-p flag)
                (flag-validation-error
                 flag "the first-in-pair bit was set in an unpaired read"
                 aln))
               ((second-in-pair-p flag)
                (flag-validation-error
                 flag "the second-in-pair bit was set in an unpaired read"
                 aln))
               (t
                flag)))
        (t
         flag)))

(defun flag-validation-error (flag message &optional aln)
  "Raised a {define-condition malformed-field-error} for alignment
FLAG in ALN, with MESSAGE."
  (if aln
      (let ((reference-id (reference-id aln))
            (read-name (read-name aln))
            (pos (alignment-position aln)))
        (error 'malformed-field-error :record aln :field flag
               :format-control (txt "invalid flag ~b set for read ~s at ~a"
                                    "in reference ~d: ~a")
               :format-arguments (list flag read-name pos reference-id
                                       message)))
      (error 'malformed-field-error
             :field flag
             :format-control "invalid flag ~b set: ~a"
             :format-arguments (list flag message))))

(defun bam-sort-error (previous-ref previous-pos ref pos &optional message
                       &rest message-arguments)
  (error 'bam-sort-error
         :prev-reference previous-ref :reference ref
         :prev-position previous-pos :position pos
         :format-control message
         :format-arguments message-arguments))
