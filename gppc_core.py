#!/usr/bin/env python

# Generalized Ping-Pong calculator

# Author: Yu Fu
# 
# This script defines the core functions.
# It is designed to be able to handle 5'-5' (Ping-Pong), 5'-3' (phasing), 3'-3' distances given a bed2 file.
# By default, it calculates the distance using all reads.
# You can restrict it by only considering reads from the same strand, or reads from different strands.

# Global variable to set how wide the range should be. e.g. if you wants to analyze -50 to +50 Ping-Pong, set it to 50.
# 0 means there is exactly 1nt overlap
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import sys
from multiprocessing import Process, Queue
# Load the bed2 file as a
# Load a ChromeInfo file to set up the numpy array
BLOCK = 10000
#     q.put( d1(s, rang, flag) )

def d2(sa, sb, rang):
    '''d2 calculates distances between sa and sb within the range defined by rang. Named d2 because it accepts two matrices...
'''
    ret = {}
    for i in range(-rang, rang+1):
        ret[i] = 0.0
    count = 0
    for i in range(sa.shape[0]):
        cur = sa[i, 0]
        count += 1
        if(count % BLOCK == 0):
            print "Processing signals: " + str(count)
        # Signals on the other bed2 file
        for j in range(-rang, rang+1):
            # Find signal on the different strand
            try:
                ret[j] += sb[i+j, 0] * cur
            except IndexError:
                print >>sys.stderr, "Encounter choromosome boundaries"
                # pass
    return ret
        
def d2_worker(sa, sb, rang, q):
    '''Just a wrapper around d2
'''
    q.put( d2(sa, sb, rang) )

def bed2_reader(fn, s):
    '''Just a bed2 reader that reads in a bed2 file and output the np matrix'''
    fh = open(fn, "r")
    counter = 0
    for line in fh.xreadlines():
        counter += 1
        if (counter % BLOCK == 0):
            print >>sys.stderr, "Processed %d lines" % counter
        chrom, start, end, copy, ntm, strand, seq = line.split()
        start = int(start)
        end = int(end)
        copy = int(copy)
        ntm = int(ntm)
        cur_sig = float(copy)/ntm
        s[chrom][start, 0] += cur_sig
    

# def d1(s, rang, flag):
#     '''5' to 5' distance or 3' to 3' distance:
# The sign of the distance is based on the reference read (the reference reads always goes from the left to the right)
# requiring all signals on the same chromosome
# flag can be "ds" (different strand), "ss" (same strand). If you want to get the signal from ss and ds, you can simply add them up
# "ds" is for Ping-Pong.
#     '''
#     ret = {}
#     for i in range(-rang, rang+1):
#         ret[i] = 0.0
#     count = 0
#     if flag == "ds":
#         for i in range(s.shape[0]):
#             count += 1
#             if(count % BLOCK == 0):
#                 print "Processing signals: " + str(count)
#             # Target strand
#             cur_sig_w = s[i, 0]
#             cur_sig_c = s[i, 1]
#             if cur_sig_w != 0:
#                 for j in range(-rang, rang+1):
#                     # Find signal on the different strand
#                     try:
#                         ret[j] += s[i+j, 1] * cur_sig_w
#                     except IndexError:
#                         print >>sys.stderr, "Encounter choromosome boundaries"
#                         # pass
#             if cur_sig_c != 0 :
#                 for j in range(-rang, rang+1):
#                     # NOTICE the negative sign here
#                     try:
#                         ret[-j] += s[i+j, 0] * cur_sig_c
#                     except IndexError:
#                         print >>sys.stderr, "Encounter choromosome boundaries"
#                         # pass
#     elif flag == "ss":
#         for i in range(s.shape[0]):
#             count += 1
#             if(count % BLOCK == 0):
#                 print "Processing signals: " + str(count)
#             # Target strand
#             cur_sig_w = s[i, 0]
#             cur_sig_c = s[i, 1]
#             if cur_sig_w != 0:
#                 for j in range(-rang, rang+1):
#                     # Find signal on the different strand
#                     try:
#                         ret[j] += s[i+j, 0] * cur_sig_w
#                     except IndexError:
#                         pass
#             if cur_sig_c != 0 :
#                 for j in range(-rang, rang+1):
#                     # NOTICE the negative sign here
#                     try:
#                         ret[-j] += s[i+j, 1] * cur_sig_c
#                     except IndexError:
#                         pass
#     else:
#         print >>sys.stderr, "Unrecognized parameter: " + flag
#     # Each pair is counted twice...
#     for i in range(-rang, rang+1):
#         ret[i] = ret[i] / 2        
#     return ret

# def d1_worker(s, rang, flag, q):
#     ''' Just a wrapper around d1() so that I can collect the return value from each process
# '''
