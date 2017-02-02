#!/usr/bin/env python

# Report all siRNAs that have its pair present. The pair should
# have a 2-nt overhang at 3' ends

import sys

bed2_fn = sys.argv[1]
overhang = 2

sirna_min = 20
sirna_max = 22

# Stores all the potential good 5' positions of siRNA pairs
d = {}
with open(bed2_fn) as fh:
    for line in fh:
        line = line.strip()
        col = line.split()
        chrom, start, end, copy, ntm, strand = col[0:6]
        start = int(start)
        end = int(end)
        copy = float(copy)
        ntm = float(ntm)

        l = end - start

        if l >= sirna_min and l <= sirna_max:
            if strand == "+":
                p = end - 1 - 2
                d[chrom + ":" + str(start-overhang) + ":" + str(end-overhang) + "-"] = line
            else:
                p = start + 2
                d[chrom + ":" + str(start+overhang) + ":" + str(end+overhang) + "+"] = line
            

with open(bed2_fn) as fh:
    for line in fh:
        line = line.strip()
        col = line.split()
        chrom, start, end, copy, ntm, strand = col[0:6]
        start = int(start)
        end = int(end)
        copy = float(copy)
        ntm = float(ntm)

        k = chrom + ":" + str(start) + ":" + str(end) + strand

        if k in d:
            # print line + "\t" + d[k]
            print line
            
