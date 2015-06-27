#!/usr/bin/env python

# Generalized Ping-Pong calculator

# Author: Yu Fu
# 
# This script is designed to be able to handle 5'-5' (Ping-Pong), 5'-3' (phasing), 3'-3' distances given a bed2 file.
# By default, it calculates the distance using all reads.
# You can restrict it by only considering reads from the same strand, or reads from different strands.

# Global variable to set how wide the range should be. e.g. if you wants to analyze -50 to +50 Ping-Pong, set it to 50.
# 0 means there is exactly 1nt overlap
RANGE=50
# When index is 0, there is just one nt overlap
hist = range(-RANGE, RANGE+1)
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import sys
import multiprocessing
# Load the bed2 file as a
# Load a ChromeInfo file to set up the numpy array

def d1(s, flag="ds"):
    '''5' to 5' distance or 3' to 3' distance:
The sign of the distance is based on the reference read (the reference reads always goes frmo the left to the right)
requiring all signals on the same chromosome
flag can be "ds" (different strand), "ss" (same strand), "all" (no restriction on strand)
"ds" is for Ping-Pong.
    '''
    ret = {}
    for i in range(-RANGE, RANGE+1):
        ret[i] = 0.0
    chrom_l = cl[chrom]
    count = 0
    if flag == "ds":
        for i in range(s.shape[0]):
            count += 1
            if(count % 10000 == 0):
                print "Processing signals: " + str(count)
            # Target strand
            cur_sig_w = s[i, 0]
            cur_sig_c = s[i, 1]
            if cur_sig_w != 0:
                for j in range(-RANGE, RANGE+1):
                    # Find signal on the different strand
                    ret[j] += s[i+j, 1] * cur_sig_w
            if cur_sig_c != 0 :
                for j in range(-RANGE, RANGE+1):
                    # NOTICE the negative sign here
                    ret[-j] += s[i+j, 0] * cur_sig_c
    elif flag == "ss":
        for i in range(s.shape[0]):
            count += 1
            if(count % 10000 == 0):
                print "Processing signals: " + str(count)
            # Target strand
            cur_sig_w = s[i, 0]
            cur_sig_c = s[i, 1]
            if cur_sig_w != 0:
                for j in range(-RANGE, RANGE+1):
                    # Find signal on the different strand
                    ret[j] += s[i+j, 0] * cur_sig_w
            if cur_sig_c != 0 :
                for j in range(-RANGE, RANGE+1):
                    # NOTICE the negative sign here
                    ret[-j] += s[i+j, 1] * cur_sig_c
    elif flag == "all":
        for i in range(s.shape[0]):
            count += 1
            if(count % 10000 == 0):
                print "Processing signals: " + str(count)
            # Target strand
            cur_sig_w = s[i, 0]
            cur_sig_c = s[i, 1]
            if cur_sig_w != 0:
                for j in range(-RANGE, RANGE+1):
                    # Find signal on the different strand
                    ret[j] += s[i+j, 0] * cur_sig_w
            if cur_sig_c != 0 :
                for j in range(-RANGE, RANGE+1):
                    # NOTICE the negative sign here
                    ret[-j] += s[i+j, 1] * cur_sig_c
    return ret

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

ret = d1(sig5["TAS"], "ss")
for i in range(-RANGE, RANGE+1):
    print str(i) + "\t" + str(ret[i])
        


