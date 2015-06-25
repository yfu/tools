#!/usr/bin/env python

# Generalized Ping-Pong calculator

# Author: Yu Fu

# This script is designed to be able to handle 5'-5' (Ping-Pong), 5'-3' (phasing), 3'-3' distances given a bed2 file.
# By default, it calculates the distance using all reads.
# You can restrict it by only considering reads from the same strand, or reads from different strands.

# Global variable to set how wide the range should be. e.g. if you wants to analyze -50 to +49 Ping-Pong, set it to 50
RANGE=100
# When index is 0, there is just one nt overlap
hist = range(-RANGE, RANGE)
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
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
# Stores signals of 5' ends and 3' end separately
sig5 = {}
sig3 = {}
for i in cl.keys():
    # 2-D array: rows are positions and 2 columns are strands
    # column 0 for the Watson strand
    # column 1 for the Crick strand
    # ndarray solution
    sig5[i] = np.zeros(shape=[ cl[i], 2 ], dtype=float)
    sig3[i] = np.zeros(shape=[ cl[i], 2 ], dtype=float)    
    # Compressed Sparse Column matrix (CSC) solution or lil
    # sig5[i] = csr_matrix((cl[i], 2), dtype=float)
    # sig3[i] = csr_matrix((cl[i], 2), dtype=float)    
print >>sys.stderr, "Chromosome lengths have been loaded"

# for i in sig:
#     print i + "\t",
#     print sig[i].shape

# Read the bed2 file
fn = sys.argv[2]
fh = open(fn, "r")
counter = 0
for line in fh.xreadlines():
    counter += 1
    if (counter % 1000 == 0):
        print >>sys.stderr, "Processed %d lines" % counter
    chrom, start, end, copy, ntm, strand, seq = line.split()
    start = int(start)
    end = int(end)
    copy = int(copy)
    ntm = int(ntm)
    cur_sig = float(copy)/ntm
    if (strand == "+"):
        sig5[chrom][start, 0] += cur_sig
        sig3[chrom][end-1, 0] += cur_sig
    elif (strand == "-"):
        sig5[chrom][end-1, 1] += cur_sig
        sig3[chrom][start, 1] += cur_sig

        
def ftfd(s5, chrom, flag="same_strand"):
    '''5' to 5' distance
    flag can be "same strand", "different strand", or "mixed strand"
    '''
    ret = range(-RANGE, RANGE)
    chrom_l = cl[chrom]
    if (flag=="same_strand"):
        for i in range(s5.shape[1]):
            sig_i_w = s5[i, 0]
            sig_i_c = s5[i, 1]
            for j in range(-RANGE, RANGE):
                p = i + j
                if p >= 0 or p < chrom_l:
                    ret[j] += sig_i_w * s5[p, 0]
                    ret[j] += sig_i_c * s5[p, 1]


