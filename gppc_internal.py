#!/usr/bin/env python

from gppc_core import d2, d2_worker
from gppc_core import bed2_reader
import sys
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
from multiprocessing import Process, Queue
import argparse
DEBUG = True

parser = argparse.ArgumentParser(description="A simple distance calculator. You need to preprocess the two bed2 files so that all lengths of all ranges in the bed2 files are 1's. Also, you need to make sure that you only include signals on the strand(s) you are interested in. You probably want to use the wrapper script instead of this scripts to calculate Ping-Pong signals and phasing patterns.\nThis script only uses the start position of each entry.")
parser.add_argument('-a', type=str, dest="a", required=True, metavar="a.bed2", help="The bed2 file containing 5' or 3' ends of small RNA reads on the plus strand")
parser.add_argument('-b', type=str, dest="b", required=True, metavar="b.bed2", help="The bed2 file containing 5' or 3' ends of small RNA reads on the minus strand")
parser.add_argument('-c', type=str, dest="c", required=True, metavar="XX.ChromInfo.txt", help="A tab-delimited file with chromosome names and lengths")
parser.add_argument('-d', dest="DEBUG", action="store_true", help="Turn on Debugging mode")
parser.add_argument('-r', '--range', dest="rang", default=50, type=int, help="The range you are interested in")
parser.add_argument('-v', '--version', action='version', version='0.1')
args = vars(parser.parse_args())
print args
if args["DEBUG"] is True:
    print "Degbugging mode is on."
    
rang = args["rang"]

# When index is 0, there is just one nt overlap
hist = {}
for i in range(-rang, rang+1):
    hist[i] = 0.0

if __name__ == "__main__":
    cl_fn = args['c']
    cl = {}
    with open(cl_fn) as fh:
        for line in fh.xreadlines():
            line = line.strip()
            col = line.split()
            cl[col[0]] = int(col[1])

    # Each chromosome is an ndarray and I store signals from the two files
    sa = {}
    sb = {}
    for i in cl.keys():
        # 2-D array: rows are positions and 2 columns are strands
        # column 0 for the Watson strand
        # column 1 for the Crick strand
        # ndarray solution
        sa[i] = np.zeros(shape=[ cl[i], 1 ], dtype=float)
        sb[i] = np.zeros(shape=[ cl[i], 1 ], dtype=float)    
        # Compressed Sparse Column matrix (CSC) solution or lil
        # sa[i] = csr_matrix((cl[i], 2), dtype=float)
        # sb[i] = csr_matrix((cl[i], 2), dtype=float)
    if DEBUG is True:
        for i in cl.keys():
            print "DEBUG: A file" + str(sa[i].shape)
            print "DEBUG: B file" + str(sb[i].shape)
    print >>sys.stderr, "Chromosome lengths have been loaded."
    
    # Read the bed2 file
    bed2_reader(args['a'], sa)
    bed2_reader(args['b'], sb)
    print >>sys.stderr, "Signals have been loaded."
    # print sa
    q = Queue()
    threads = []

    if DEBUG == True:
        print "Multithreading begins..."
    # for i in ("X-TAS",):
    for i in ("chr2LHet",):
    # for i in cl.keys():
        p = Process(target=d2_worker, args=(sa[i], sb[i], rang, q))
        p.start()
        threads.append(p)
    for t in threads:
        t.join()
    print q.qsize()
    while(q.empty() == False):
        r = q.get()
        print r
        for i in range(-rang, rang+1):
            hist[i] += r[i]
    print hist
        
    for i in range(-rang, rang+1):
        print str(i) + "\t" + str(hist[i])



    


