#!/usr/bin/env python

# Given a file in bed2 format, this script will calculate the number of reads in non-overlapping 
# windows and spit out a bed format defining piRNA clusters.

# genome.length has one line for each chromosome: the first column is the name and the second column is the length of this chromosom.
# Usage python pirna_cluster_caller.py genome.length your.bed2> pirna_definition.bed
# Check if it is correct by awk '{ d[$1] += $4/$5 } END{ for(i in d) {print i, d[i]} }' test.1000 and see if the signals match.

import sys
import argparse

win_size = 1000
# Threshold: 5 reads every 1000 bp.
threshold = 5
parser = argparse.ArgumentParser(description='Call piRNA clusters')

parser.add_argument('-f', '--file', help='bed2 file', required=True)
parser.add_argument('-g', '--genome', help='A genome file containing chromosomes as the first column and the lenght of that chromosome as the second column', required=True)
parser.add_argument('-u', '--unique', help='Use unique mappers only', required=False, action="store_true")
args = parser.parse_args()

chrom_len_fn = args.genome
bed2_fn = args.file

chrom_len = {}

with open(chrom_len_fn, "r") as f:
    for line in f.readlines():
        line = line.strip()
        [chrom, length] = line.split()
        if chrom in chrom_len:
            print >>sys.stderr, "You have duplicate chromosomes. Exiting..."
            sys.exit(1)
        chrom_len[chrom] = length
# print chrom_len
# all_signals is a dict. Its keys are chromosome and the value is a list. The list stores the signals in each window.
all_signals = {}
# Allocate the keys:
for c in chrom_len.keys():
    tmp_len = chrom_len[c]
    tmp_n_win = float(tmp_len) / win_size
    # Take care of both cases: say a chromosome is 5000bp long, you should divide it into 5 parts. If the chromosome 
    # is 5001bp long, you should divide it into 6 windows.
    if (tmp_n_win.is_integer()):
        all_signals[c] = [0] * int(tmp_n_win)
    else:
        all_signals[c] = [0] * (int(tmp_n_win) + 1)
print >>sys.stderr,  "Finish reading the genome file."

f = open(bed2_fn, "r")
for line in f.readlines():
    line = line.strip()
    e = line.split()
    chrom = e[0] 
    start = int(e[1])
    end = int(e[2])
    copy = int(e[3])
    ntm = int(e[4])
    strand = e[5]
    seq = e[6]
    # print ", " . join(e)
    if ntm != 1 and args.unique == True:
        continue
    signal = float(copy) / ntm
    # Only consider start of the read to make the calculation simple
    win_n = start / win_size
    all_signals[chrom][win_n] += signal

# Now output those windows with more than n reads
for k in all_signals.keys():
    cur_chrom = all_signals[k]
    for i in range(len(cur_chrom)):
        if cur_chrom[i] > threshold:
            print "\t" . join( [k, str(i * win_size), str((i + 1) * win_size), ".", str(cur_chrom[i]), "."] )
