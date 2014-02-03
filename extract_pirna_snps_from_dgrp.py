#!/usr/bin/env python

# Extract the SNPs inside the piRNA clsuters
# From SNP files to bed files
import sys
pirna_loc_fn = '/data/fuy2/flycommon/Brennecke.pirnaCluster.bed'

pirna_loc = []
fh = open(pirna_loc_fn)
for line in fh.readlines():
    pirna_loc.append(line.strip().split())
fh.close()

chrom = ['2L', '2R', '3L', '3R', 'X']
# chrom = ['2L']
pre = "Variants_Sparse_"
suf = ".sample_swap_fixed.txt"
for i in chrom:
    fh = open(pre + i + suf)
    out = open(i + "_snps.bed", 'w')
    fh.readline() # Ignore the first line
    for line in fh.readlines():
        e = line.strip().split(',')
        print >> out, "chr" + i + "\t" + e[0] + "\t" + e[0] + '\t' + ','.join(e[1:-1])
        

