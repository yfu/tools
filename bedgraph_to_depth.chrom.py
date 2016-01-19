#!/usr/bin/env python

# The same with bedgraph_to_depth.py except that this one also accept
# Given one bedgraph file, this script output a vector-like file describing depth.
# This script glues the bedgraph file with simple plotting functions in R.
# start and end parameters are like those in bed files, i.e. start from 0 and half open
# half closed
# Usage: bedgraph_to_depth.py my.bedGraph chrom 0 1000

import sys
import re

pat = re.compile(r"^track ")

if(len(sys.argv) !=5):
    print "Usage: %s my.bedgraph 0 1234" % sys.argv[0]
    print "Note that start and end are like those in bed files: start from 0 and half open half closed"
    exit(1)
bg_fn = sys.argv[1]
bg = open(bg_fn, "r")
chrom = sys.argv[2]
start = int(sys.argv[3])
end   = int(sys.argv[4])

i = start
# chrom = ""
for line in bg.readlines():
    if pat.match(line) is not None:
        print >>sys.stderr, "Found a track line in the bedGraph file. Skipping:"
        print >>sys.stderr, line
        continue
    if i > end:
        break
    # print line
    col = line.split()
    my_chrom = col[0]
    s = int(col[1])
    e = int(col[2])
    sig = float(col[3])
    if my_chrom == chrom:
        while(i<s and i<end):
            print chrom + "\t" + str(i) + "\t" + "0"
            i += 1
        while(i<e and i<end):
            # print i
            print chrom + "\t" + str(i) + "\t" + str(sig)        
            i += 1
# For the last part of the gene (which might not have signals), we need to output this part
# as well
while(i < end):
    print chrom + "\t" + str(i) + "\t" + "0"
    i += 1
    

