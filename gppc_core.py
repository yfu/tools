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

def d1(s, rang, flag):
    '''5' to 5' distance or 3' to 3' distance:
The sign of the distance is based on the reference read (the reference reads always goes frmo the left to the right)
requiring all signals on the same chromosome
flag can be "ds" (different strand), "ss" (same strand), "all" (no restriction on strand)
"ds" is for Ping-Pong.
    '''
    ret = {}
    for i in range(-rang, rang+1):
        ret[i] = 0.0
    count = 0
    if flag == "ds":
        for i in range(s.shape[0]):
            count += 1
            if(count % 1000000 == 0):
                print "Processing signals: " + str(count)
            # Target strand
            cur_sig_w = s[i, 0]
            cur_sig_c = s[i, 1]
            if cur_sig_w != 0:
                for j in range(-rang, rang+1):
                    # Find signal on the different strand
                    try:
                        ret[j] += s[i+j, 1] * cur_sig_w
                    except IndexError:
                        # print >>sys.stderr, "Encounter choromosome boundaries"
                        pass
            if cur_sig_c != 0 :
                for j in range(-rang, rang+1):
                    # NOTICE the negative sign here
                    try:
                        ret[-j] += s[i+j, 0] * cur_sig_c
                    except IndexError:
                        pass
                        # print >>sys.stderr, "Encounter choromosome boundaries"
    elif flag == "ss":
        for i in range(s.shape[0]):
            count += 1
            if(count % 1000000 == 0):
                print "Processing signals: " + str(count)
            # Target strand
            cur_sig_w = s[i, 0]
            cur_sig_c = s[i, 1]
            if cur_sig_w != 0:
                for j in range(-rang, rang+1):
                    # Find signal on the different strand
                    try:
                        ret[j] += s[i+j, 0] * cur_sig_w
                    except IndexError:
                        pass
            if cur_sig_c != 0 :
                for j in range(-rang, rang+1):
                    # NOTICE the negative sign here
                    try:
                        ret[-j] += s[i+j, 1] * cur_sig_c
                    except IndexError:
                        pass
    elif flag == "all":
        for i in range(s.shape[0]):
            count += 1
            if(count % 1000000 == 0):
                print "Processing signals: " + str(count)
            # Target strand
            cur_sig_w = s[i, 0]
            cur_sig_c = s[i, 1]
            if cur_sig_w != 0:
                for j in range(-rang, rang+1):
                    # Find signal on the different strand
                    try:
                        ret[j] += s[i+j, 0] * cur_sig_w
                    except IndexError:
                        pass
            if cur_sig_c != 0 :
                for j in range(-rang, rang+1):
                    # NOTICE the negative sign here
                    try:
                        ret[-j] += s[i+j, 1] * cur_sig_c
                    except IndexError:
                        pass
    return ret

def d1_worker(s, rang, flag, q):
    ''' Just a wrapper around d1() so that I can collect the return value from each process
'''
    q.put( d1(s, rang, flag) )

def d2(s5, s3, rang, flag):
    '''d2 calculates 5' to 3' distances. Named d2 because it accepts two matrices...
'''
    
