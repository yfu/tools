#!/usr/bin/env python

# usage: python assign_pirna_to_genes.py SRAtoTranscripts.all.ox.bed2 > SRAtoTranscripts.all.ox.count
# Given a bed2 file (small RNAs mapped to transcriptome), this script outputs the number of
# reads each Trinity gene has (unique reads and nonunique reads)

import sys
import re
from collections import defaultdict

fn = sys.argv[1]

# key is the read sequence; value is the copy number
read2copy = {}
# key is the read sequence; value is the list of gene name(s)
read2gene = {}
# all gene names
genes = {}
with open(fn) as fh:
    for line in fh.xreadlines():
        line = line.strip()
        col = line.split()
        transcript, start, end, copy, ntm, strand, seq = col
        start = int(start)
        end = int(end)
        copy = int(copy)
        ntm = int(ntm)
        mat = re.match("(\S+)_(i[0-9]+)", transcript)
        if mat == None:
            print >>sys.stderr, "This file does not seem have \
            transcript names as the first column"
            print >>sys.stderr, line
            sys.exit(1)
        gene = mat.group(1)
        isoform = mat.group(2)
        genes[gene] = 1
        ## print gene + "\t" + isoform + "\t" + seq
        read2copy[seq] = copy
        if seq not in read2gene:
            read2gene[seq] = set([gene, ])
        else:
            read2gene[seq] = read2gene[seq] | set([gene, ])

uniq_counts = defaultdict(int)
nonuniq_counts = defaultdict(float)

for read in read2gene:
    # print read + "\t" + str(read2gene[read]) + "\t" + str(read2copy[read])
    g = list(read2gene[read])
    c = read2copy[read]
    if len(g) == 1:
        uniq_counts[g[0]] += c
    else:
        for i in g:
            nonuniq_counts[i] += float(c) / len(g)
print "#gene\tuniq\tnon-uniq"
for g in genes:
    print g + "\t" + str(uniq_counts[g]) + "\t" + str(nonuniq_counts[g])
