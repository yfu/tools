#!/usr/bin/env python

# Generalized Ping-Pong calculator

# Author: Yu Fu

# This script is designed to be able to handle 5'-5' (Ping-Pong), 5'-3' (phasing), 3'-3' distances given a bed2 file.
# By default, it calculates the distance using all reads.
# You can restrict it by only considering reads from the same strand, or reads from different strands.

import numpy as np
import sys
# Load the bed2 file as a
# Load a ChromeInfo file to set up the numpy array

cl_fn = sys.argv[1]
cl = {}
with open(cl_fn) as fh:
    for line in fh.xreadlines():
        line = line.strip()
        col = line.split()
        cl[col[0]] = int(col[1])

# Each chromosome is an ndarray
sig = {}
for i in cl.keys():
    # 2-D array: rows are positions and 2 columns are strands (watson, then crick)
    sig[i] = np.ndarray(shape=[ cl[i], 2 ], dtype=float)
  
print >>sys.stderr, "Chromosome lengths have been loaded"

for i in sig:
    print i + "\t",
    print sig[i].shape

# Read the bed2 file
fn = sys.argv[2]
chrom, start, end, copy, ntm, strand, seq 

