#!/usr/bin/env python

# Given a bam file, this script outputs the multimappers based on the NH tag
import argparse
import pysam as ps

parser = argparse.ArgumentParser(description='A script to output the multimappers in a bam file')
parser.add_argument('-f', '--file', help='the input file', required=True)
# parser.add_argument('-o', '--output', help='the output file', required=True)
args = parser.parse_args()

infile = args.file
bam = ps.AlignmentFile(infile, "rb")
# out = ps.AlignmentFile('-', "wb", template=bam)
out = ps.AlignmentFile('-', "wb", template=bam)

for read in bam.fetch():
    tag = read.get_tag("NH")
    if(tag != 1):
        out.write(read)
