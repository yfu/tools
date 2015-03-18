#!/usr/bin/env python

# Convert a bed6 file into a bedgraph file. It works (almost) in the same way as bedtools bamtobed, but it also takes multimappers into consideration
# This input SAM must have the NH tag reported by the aligner

# NOTE: You need to first convert a BAM file to a bed12 file, then reverse one of the reads in a pair (for dUTP, reverse the r1), then convert such a bed12 file to a bed6 file:
# END_TO_REVERSE_STRAND=1
# bedtools bamtobed -bed12 -tag NH -i ${input} | awk -v strand=$END_TO_REVERSE_STRAND 'BEGIN{FS=OFS="\t"}{e=substr($4,length($4)); if (e==strand) $6=($6=="+"?"-":"+"); print $0; }' 
# Usage: bed6_to_bedgraph.py ....

import HTSeq
import sys

gl_fn = sys.argv[1]
gl = {}

with open(gl_fn) as fh:
    for line in fh.xreadlines():
        line = line.strip()
        col = line.split()
        gl[col[0]] = int(col[1])
print >>sys.stderr, "Chromosome lengths have been loaded"

# print >>sys.stderr, "Start loading the SAM file..."
# for line in fileinput.input():
fn = sys.argv[2]
b = HTSeq.BED_Reader(fn)
# Use ndarray since the data will be constantly changing, unlike genome annotations
ga = HTSeq.GenomicArray(chroms=gl, stranded=True, typecode="d", storage="ndarray")

# Output prefix
output_pf = sys.argv[3]
# Output filename prefix
print >>sys.stderr, "Start loading the input bed file...."
for aln in b:
    iv = aln.iv
    if iv.strand != "+" and iv.strand != "-":
        print >>sys.stderr, "Found an alignment w/o strand information. I will ignore it..."
        continue
    copy = 1
    ntm = int(aln.score)
    impact = float(copy) / ntm
    ga[iv] += impact

print >>sys.stderr, "Start writing bedgraph files"
watson_bg = output_pf +  ".Watson.bedGraph"
crick_bg = output_pf + ".Crick.bedGraph"
ga.write_bedgraph_file(watson_bg, strand="+")
ga.write_bedgraph_file(crick_bg, strand="-")
