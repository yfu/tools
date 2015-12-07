#!/usr/bin/env python

# Part of BP discovery pipeline
# For each line in the bam file, count the number of matches (not mismatches, since the CIGAR M flag may represent match or mismatch)
# Output each seqname and the corresponding # matches
# Note that the input sam/bam file is generated using bowtie2 local mode, so NM tag is not useful for counting matches.

import pysam as ps
from Bio.Seq import Seq
import sys
import argparse
import sys

DEBUG = False
parser = argparse.ArgumentParser(description='A tools to output the lengths of the aligned portions of reads in a SAM/BAM file')
# parser.add_argument('-s', '--strand-ref', help='enforcing the reads to be on the reference strand', required=False, dest="sr", default=False, action="store_true")
parser.add_argument('-f', '--file', help='the input file', required=True)
parser.add_argument('-d', '--debug', help='debug mode', required=False, action="store_true", default=False)
args = parser.parse_args()
# True: output the sequences on the reference strand
# False: output the sequences in the orignal direction
# sr = args.sr
infile = args.file
DEBUG = args.debug
bam = ps.AlignmentFile(infile, "rb")
nas = 0
ns = 0

def rc(s):
    return str(Seq.reverse_complement(Seq(s)))

c = 0
for read in bam.fetch():
    read_alignment_len = read.query_alignment_length
    # M 0
    # I 1
    # D 2
    # S 4
    if DEBUG:
        print >>sys.stderr, read.cigartuples
    indel_len = 0
    m_mm_len = 0
    for operation, length in read.cigartuples:
        if operation == 1 or operation == 2:
            indel_len += length
        if operation == 0:
            m_mm_len += length
    if DEBUG:
        print >>sys.stderr, "Total length of indel " + str(indel_len)
        print >>sys.stderr, "Total length of match and mismatch " + str(m_mm_len)
    ed = 0
    for tag, val in read.get_tags():
        if tag == "NM":
            ed = val
            break
    if DEBUG:
        print >>sys.stderr, "Edit distance " + str(ed)
    # number of mismatches
    n_mm = ed - indel_len
    if DEBUG:
        print >>sys.stderr, "Length of mismatches " + str(n_mm)
    n_m = m_mm_len - n_mm
    # c += 1
    # l = read.query_alignment_end - read.query_alignment_start
    # print read.query_name + "\t" + str(l)
    print read.query_name + "\t" + str(n_m)
print >> sys.stderr, "Number of alignments processed " + str(c)
