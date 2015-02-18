#!/usr/bin/env python

# Parse the output from pps_simple.py

import sys
import numpy

fn = sys.argv[1]
fh = open(fn)

lst = []
tenth_count = 0.0

highest_pos = -1
highest_sig = -1.0

for line in fh.readlines():
    line = line.strip()
    col = line.split()
    pos = int(col[0])
    sig = float(col[1])

    if pos == 10:
        tenth_count = sig
    else:
        lst.append(sig)
    if sig > highest_sig:
        highest_sig = sig
        highest_pos = pos

if numpy.std(lst) != 0:
    print "Raw counts at 10th" + "\t" + str(tenth_count)
    pps = (tenth_count - numpy.mean(lst)) / numpy.std(lst)
    print "Ping-Pong score" + "\t" + str(pps)
    
    print "Highest position" + "\t" + str(highest_pos)
    print "Highest signal" + "\t" + str(highest_sig)
else:
    print "Cannot calculate Ping-Pong score: standard deviation is 0"
