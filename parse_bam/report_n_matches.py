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

parser = argparse.ArgumentParser(description='A tools to output the lengths of the aligned portions of reads in a SAM/BAM file')
# parser.add_argument('-s', '--strand-ref', help='enforcing the reads to be on the reference strand', required=False, dest="sr", default=False, action="store_true")
parser.add_argument('-f', '--file', help='the input file', required=True)
args = parser.parse_args()
# True: output the sequences on the reference strand
# False: output the sequences in the orignal direction
# sr = args.sr
infile = args.file

bam = ps.AlignmentFile(infile, "rb")
nas = 0
ns = 0

def rc(s):
    return str(Seq.reverse_complement(Seq(s)))

c = 0
for read in bam.fetch():
    c += 1
    l = read.query_alignment_end - read.query_alignment_start
    print read.query_name + "\t" + str(l)
print >> sys.stderr, "Number of alignments processed " + str(c)
