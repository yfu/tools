#!/usr/bin/env python

# Parse the output from gtf_db_to_length.py
import sys
import numpy

d = {}
with open(sys.argv[1]) as f:
    for line in f.readlines():
        line = line.split()
        g, _, l = line
        l = float(l)
        if g not in d:
            d[g] = [l]
        else:
            d[g].append(l)
for g in d:
    print "%s\t%.2f" % (g, numpy.median(d[g]))
