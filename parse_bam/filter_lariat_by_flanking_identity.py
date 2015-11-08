#!/usr/bin/env python

# Part of BP discovery pipeline
# This script is to eliminate spurious lariats by checking their similarity to the genome sequence
# First, all putative lariats are mapped to the genome to get a bam file.
# Then, given a bam file, this script will output the length of matching portion on the left and on the right
# lariat_name	80	121

# Usage: python ~/repo/tools/parse_bam/filter_lariat_by_flanking_identity.py -a 100 -f Mattick.RNaseR.HeLa.BP.rep1.fq.gz.map_to_put_lar.bam > good.bam 2>lariat_support.stats
#
# Genome:         --------------------------------------------------A-------------------------------------------------
# One good lariat -------------------------------------MMMMMMMMMMMMMAMM-----------------------------------------------
# One bad lariat  -------------------------------------MMMMMMMMMMMMMAMMMMMMMMMMMMM------------------------------------
# Output (note that the BP is included in the upstream (left)):
# lariat1	14	2
# lariat2	14	13

import pysam as ps
from Bio.Seq import Seq
import sys
import argparse
import sys

parser = argparse.ArgumentParser(description='A script to output the start and end position of the aligned reads')
# parser.add_argument('-a', '--anchor', help='(1-based) position of the anchor point on the read (for example, in bp discover, the parameter should be set to 100 ([0:100) is the lariat and upstream and [100:200) is 5\' intron ', required=True, dest="ap", type=int)
parser.add_argument('-f', '--file', help='the input file', required=True)
# parser.add_argument('-o', '--output', help='the output file', required=True)
args = parser.parse_args()
# True: output the sequences on the reference strand
# False: output the sequences in the orignal direction
# ap  = args.ap - 1                         # 0-based postiion

infile = args.file
bam = ps.AlignmentFile(infile, "rb")
# out = ps.AlignmentFile('-', "wb", template=bam)
out = ps.AlignmentFile('-', "wb", template=bam)
nas = 0
ns = 0

# names of each supporting reads
readname_sup = {}
# seq
seq_sup = {}
for r in bam.references:
    readname_sup[r] = []
    seq_sup[r] = []

def rc(s):
    return str(Seq.reverse_complement(Seq(s)))

def print_sam_line(r):
    print 
for read in bam.fetch():
    ref_name = bam.getrname(read.reference_id)
    read_name = read.query_name
    read_seq = read.seq

    if read.is_reverse:
        strand = "-"
    else:
        strand = "+"
    # Left coordinate
    # start index of the aligned query portion of the sequence (0-based, inclusive).
    lc = read.query_alignment_start
    # end index of the aligned query portion of the sequence (0-based, exclusive)
    rc = read.query_alignment_end - 1     # (0-based, inclusive)
    print read_name + "\t" + strand + "\t" + str(lc) + "\t" + str(rc)
