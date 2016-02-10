#!/usr/bin/env python

# This is for searching potential targets in DEG
# Given a bam file from BWA with many mismatches allowed (3, 4, 5, 6, 7), this script reports
# the 5' position of the reads, and normalized read counts

import pysam as ps
from Bio.Seq import Seq
import sys
import argparse
import sys

DEBUG = False
parser = argparse.ArgumentParser(description='A tools to parse a bam file from BAM for obtaining potential target sites')
# parser.add_argument('-s', '--strand-ref', help='enforcing the reads to be on the reference strand', required=False, dest="sr", default=False, action="store_true")
parser.add_argument('-f', '--file', help='the input BAM file from BWA', required=True)
parser.add_argument('-u', '--unique', help='consider unique mappers only', action="store_true")
parser.add_argument('-d', '--debug', help='debug mode', required=False, action="store_true", default=False)
args = parser.parse_args()
# True: output the sequences on the reference strand
# False: output the sequences in the orignal direction
# sr = args.sr
infile = args.file
DEBUG = args.debug
require_unique = args.unique
bam = ps.AlignmentFile(infile, "rb")
nas = 0
ns = 0

def rc(s):
    return str(Seq.reverse_complement(Seq(s)))

# Count the number of species and reads at each postion
pos2species = {}
pos2read = {}
for read in bam.fetch():
    read_alignment_len = read.query_alignment_length
    # M 0
    # I 1
    # D 2
    # S 4
    if DEBUG:
        print >>sys.stderr, read.cigartuples
    nm_tag = read.get_tag("NM")
    if read.is_reverse:
        fiveprime = read.reference_end - 1
    else:
        fiveprime = read.reference_start
    if require_unique:
        if nm_tag == 1
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
