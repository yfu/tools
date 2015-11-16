#!/usr/bin/env python

# Usage: python mark_duplicates_umi.py -f test.sorted.bam >test.dup_marked.bam 2>test.dup_marked.log

DEBUG = False
# A FASTQ file is first processed by reformat_umi_fastq.py to put the
# barcode info in the fasta header. Then it is aligned by aligner such as
# STAR. The bam file can then be processed by this script.
# It will go through the bam file and search for r1 reads that align to the same
# position and 
import pysam as ps
from Bio.Seq import Seq
import sys
import argparse
import sys

parser = argparse.ArgumentParser(description='A pair of FASTQ files are first reformatted using reformat_umi_fastq.py and then is aligned to get the bam file. This script can parse the barcode in the bam file to mark duplicates. For now, it is assumed that any r1 reads with the same 5\' ends and the same barcode are duplicates and only one of them is kept.')
parser.add_argument('-f', '--file', help='the input bam file', required=True)
parser.add_argument('-d', '--debug', help='turn on debug mode', action="store_true")
# parser.add_argument('-o', '--output', help='the output file', required=True)
args = parser.parse_args()
# True: output the sequences on the reference strand
# False: output the sequences in the orignal direction

infile = args.file
bam = ps.AlignmentFile(infile, "rb")
out = ps.AlignmentFile('-', "wb", template=bam)
DEBUG = args.debug

# Stores the reads to be marked as duplicates (only the 2nd and later ones are marked; the first one is not)
dup_read_names = []

# Position + barcode as the key; count as the value
posbc = {}

for read in bam.fetch():
    # Ignore r2 in the first pass
    if read.is_read2:
        continue
    elif read.is_paired == False:
        sys.exit("This is not paired end data. Exiting...")
    # out.write(read)
    read_n = read.query_name
    read_bc = read_n.split("_")[0]
    read_chr = read.reference_id
    if read.is_reverse:
        # read_5 is the position of the 5'end of the read on the reference
        strand = "-"
        read_5 = read.reference_end - 1
    else:
        strand = "+"
        read_5 = read.reference_start
    if DEBUG:
        print >>sys.stderr, "5'end: " + str(read_5)
        print >>sys.stderr, "Reference ID: " + str(read.reference_id)
        print >>sys.stderr, "Barcode: " + read_bc
        print >>sys.stderr, "Read: " + str(read)
    read_id = str(read_chr) + "|" + str(read_5) + "|" + read_bc + "|" + strand
    if read_id in posbc:
        posbc[read_id] += 1
        if DEBUG:
            print >>sys.stderr, "Duplicates found:" + str(read)        
        dup_read_names.append(read_n)
    else:
        posbc[read_id] = 1

bam.close()
# The 2nd pass:
bam = ps.AlignmentFile(infile, "rb")
for read in bam.fetch():
    if read.query_name in dup_read_names:
        # Add 0x400 flags for PCR duplicates
        read.flag = read.flag + 0x400
        out.write(read)
    else:
        out.write(read)
