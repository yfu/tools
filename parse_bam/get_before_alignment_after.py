#!/usr/bin/env python

# Part of BP discovery pipeline
# 1st col is before; 2nd column is the aligned part; 3rd column is after
# By default, the reported sequences are in the same direction of the original reads
# adding options -s to enforce the output to be on the reference strand, that is, if a
# seq maps to the other strand, this seq will be reverse complemented and reported

# Notice that bowtie and bowtie2 output sams files in which the seq field is always in the same direction with the reference
import pysam as ps
from Bio.Seq import Seq
import sys
import argparse
import sys

parser = argparse.ArgumentParser(description='A tools to output the soft-clipped sequences on both ends')
parser.add_argument('-s', '--strand-ref', help='enforcing the reads to be on the reference strand', required=False, dest="sr", default=False, action="store_true")
parser.add_argument('-f', '--file', help='the input file', required=True)
args = parser.parse_args()
# True: output the sequences on the reference strand
# False: output the sequences in the orignal direction
sr = args.sr
infile = args.file

bam = ps.AlignmentFile(infile, "rb")
nas = 0
ns = 0

def rc(s):
    return str(Seq.reverse_complement(Seq(s)))

for read in bam.fetch():
    strand = "s"
    if read.is_reverse:
        nas += 1
        strand = "as"
    else:
        ns += 1
    us = read.query_sequence[0:read.query_alignment_start]
    ds = read.query_sequence[read.query_alignment_end:]
    al = read.query_alignment_sequence
    # Print the mapping and the query information
    print bam.getrname(read.reference_id) + "|" + str(read.reference_start) + "|" + str(read.reference_end) + "|" + strand + "|" + read.query_name,    
    if sr:
        print us + "\t",
        print al + "\t",
        print ds + "\t"
    else:
        print rc(ds) + "\t",
        print rc(al) + "\t",
        print rc(us) + "\t"
print >> sys.stderr, "Number of reads mapping to the sense strand " + str(ns)
print >> sys.stderr, "Number of reads mapping to the antisense strand " + str(nas)
