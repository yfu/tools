#!/usr/bin/env python

# Extend (or shrink) 3' ends of the ranges defined in a 
# bed file (and take strandedness into consideration), such that
# the ranges are of the same length


import sys

# Desired length
d_len = int(sys.argv[1])
for line in sys.stdin:
    line = line.strip()
    e = line.split()
    chr = e[0]
    start = int(e[1])
    end = int(e[2])
    name = e[3]
    val = e[4]
    strand = e[5]
    if(len(e)>6):
        rest = e[6:]

    if(strand == "+"):
        end = start + d_len
    else:
        start = end - d_len
    if(start<0):
        continue
    if(len(e) <=6):
        print "\t".join([chr, str(start), str(end), name, val, strand])
    else:
        print "\t".join([chr, str(start), str(end), name, val, strand] + rest)
