#!/usr/bin/env python

# Part of BP discovery pipeline
# Given a bam file, this script will output the files that overlap the branchpoint.

# Usage: python ~/repo/tools/parse_bam/filter_ol_bp.py -a 100 -l 10 -r 10 -f Mattick.RNaseR.HeLa.BP.rep1.fq.gz.map_to_put_lar.bam
#
# Lariat:    --------------------------------------------------A-------------------------------------------------
# Good read: -------------------------------------MMMMMMMMMMMMMAMMMMMMMMMMMMM------------------------------------
# Bad read:  -------------------------------MMMMMMMMMMMMMMMMMMMAM------------------------------------------------
# Bad read:  --------------------------------------------------AMMMMMMMMMMMMMMMMMMMM-----------------------------
# Bad read:  --------------------MMMMMMMMMMMMMMMMMMMMM-----------------------------------------------------------

import pysam as ps
from Bio.Seq import Seq
import sys
import argparse
import sys

parser = argparse.ArgumentParser(description='A script to ')
parser.add_argument('-a', '--anchor', help='(1-based) position of the anchor point on the read (for example, in bp discover, the parameter should be set to 100 ([0:100) is the lariat and upstream and [100:200) is 5\' intron ', required=True, dest="ap", type=int)
parser.add_argument('-l', '--left-overlap', help='enforcing length of alignment (matches/mismatches) on the left', required=True, dest="lol", type=int)
parser.add_argument('-r', '--right-overlap', help='enforcing length of alignment (matches/mismatches) on the right', required=True, dest="rol", type=int)
parser.add_argument('-f', '--file', help='the input file', required=True)
# parser.add_argument('-o', '--output', help='the output file', required=True)
args = parser.parse_args()
# True: output the sequences on the reference strand
# False: output the sequences in the orignal direction
ap  = args.ap - 1                         # 0-based postiion
lol = args.lol
rol = args.rol

infile = args.file
bam = ps.AlignmentFile(infile, "rb")
# out = ps.AlignmentFile('-', "wb", template=bam)
out = ps.AlignmentFile('-', "wb", template=bam)
nas = 0
ns = 0

def rc(s):
    return str(Seq.reverse_complement(Seq(s)))

def print_sam_line(r):
    print 
for read in bam.fetch():
    strand = "s"
    if read.is_reverse:
        nas += 1
        strand = "as"
    else:
        ns += 1
    # Left coordinate
    # start index of the aligned query portion of the sequence (0-based, inclusive).
    lc = read.reference_start
    # end index of the aligned query portion of the sequence (0-based, exclusive)
    rc = read.reference_end - 1     # (0-based, inclusive)
    # print "%d\t%d" % (lol, rol)
    if ap - lc >= lol and rc - ap >= rol:
        # print str(lc) + "\t" + str(rc) + "\t" + read.cigarstring
        # print read.tostring(bam)
        out.write(read)
    # al = read.query_alignment_sequence
