#!/usr/bin/env python

# Given a bed2 file from STDIN, this script will take the # of copies and the NTMs into consideration 
# and generate the coverage for all chromosomes (or, contigs, in my case)...

import argparse
import sys

parser = argparse.ArgumentParser(description='A tools to calculate bed2 coverage')

parser.add_argument('-s', '--strand', help='The fasta file', required=True)

args = parser.parse_args()
requested_strand = args.strand
coverages = {}

for line in sys.stdin.readlines():
    line = line.strip()
    line = line.split()
    chrom = line[0]
    start = int(line[1])
    end = int(line[2])
    copy = float(line[3])
    ntm = float(line[4])
    strand = line[5]
    seq = line[6]
    
    effective_copy = copy / ntm
    if chrom not in coverages.keys():
        coverages[chrom] = {}
    coverages[chrom][1] = 2

    for i in range(start, end):
        # The following three lines will assign 0 coverage to the regions that are covered by the other strand
        # For example, if a read covers 1000-1200 on the forward strand, this script will record 1000-1200 on the reverse
        # strand as 0.
        if(strand != requested_strand):
            if i not in coverages[chrom].keys():
                 coverages[chrom][i] = 0
        else:
            if i in coverages[chrom].keys():
                coverages[chrom][i] += effective_copy
            else: 
                coverages[chrom][i] = effective_copy           

for i in sorted(coverages.keys()):
    cur_dict = coverages[i]
    m = max(cur_dict.keys())
    for j in range(0, m+1):
        if j not in cur_dict.keys():
            print str(i) + "\t" + str(j) + "\t" + str(0)
        else:
            print str(i) + "\t" + str(j) + "\t" + str(cur_dict[j])


