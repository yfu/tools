#!/usr/bin/env python

# Convert a bed2 file to a bedGraph file.
# Only 5' end is considered
# This script requires the bed2 file to be sorted by chromosome coordinate
# Usage: bed2_to_bedgraph.py Zamore.SRA.total.00dpp.testis.trimmed.x_rRNA.x_hairpin.mm9v1.unique.siRNA.sorted.bed2 watson.bedgraph crick.bedgraph

import sys

fn = sys.argv[1]
fh = open(fn)

output_watson = open(sys.argv[2], "w")
output_crick = open(sys.argv[3], "w")

cur_chrom_watson = ""
cur_p_watson = ""
cur_signal_watson = 0

cur_chrom_crick = ""
cur_p_crick = ""
cur_signal_crick = 0

for line in fh:
    line = line.strip()
    col = line.split()
    chrom, start, end, copy, ntm, strand, seq = col
    start = int(start)
    end = int(end)
    copy = int(copy)
    ntm = int(ntm)
    if strand == "+":
        if cur_p_watson == "":
            cur_chrom_watson = chrom
            cur_p_watson = start
            cur_signal_watson += float(copy) / ntm
        elif cur_p_watson == start and cur_chrom_watson == chrom:
            cur_signal_watson += copy/ntm
        else:
            print >> output_watson, "\t".join([cur_chrom_watson, str(cur_p_watson), str(cur_p_watson + 1), str(cur_signal_watson)]) 
            cur_p_watson = start
            cur_chrom_watson = chrom
            cur_signal_watson = float(copy) / ntm
    else:
        if cur_p_crick == "":
            cur_chrom_crick = chrom
            cur_p_crick = start
            cur_signal_crick += float(copy) / ntm
        if cur_p_crick == start and cur_chrom_crick == chrom:
            cur_signal_crick += copy/ntm
        else:
            print >> output_crick, "\t".join([cur_chrom_crick, str(cur_p_crick), str(cur_p_crick + 1), str(-cur_signal_crick)])
            cur_p_crick = start
            cur_chrom_crick = chrom
            cur_signal_crick = float(copy) / ntm
