#!/usr/bin/env python

# For test data, see gppc_test_data directory
from gppc_core import d2
from gppc_core import bed2_reader
import sys
import numpy as np
# from scipy.sparse import csc_matrix
# from scipy.sparse import csr_matrix
# from scipy.sparse import lil_matrix
from multiprocessing import Process, Queue, Pool
import argparse
DEBUG = False
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="READ CAREFULLY BEFORE YOU USE IT!\nA simple distance calculator. You need to preprocess the two bed2 files so that all lengths of all ranges in the bed2 files are 1's. Also, you need to make sure that you only include signals on the strand(s) you are interested in. You probably want to use the wrapper script instead of this scripts to calculate Ping-Pong signals and phasing patterns.\nThis script only uses the start position of each entry.")
    parser.add_argument('-a', type=str, dest="a", required=True, metavar="a.bed2", help="The bed2 file containing 5' or 3' ends of small RNA reads on the plus strand")
    parser.add_argument('-b', type=str, dest="b", required=True, metavar="b.bed2", help="The bed2 file containing 5' or 3' ends of small RNA reads on the minus strand")
    parser.add_argument('-c', type=str, dest="c", required=True, metavar="XX.ChromInfo.txt", help="A tab-delimited file with chromosome names and lengths")
    parser.add_argument('-d', dest="DEBUG", action="store_true", help="Turn on Debugging mode")
    parser.add_argument('-p', dest="cpu", type=int, default=4, help="Number of threads")    
    parser.add_argument('-r', '--range', dest="rang", default=50, type=int, help="The range you are interested in (the default is 50)")
    parser.add_argument('-v', '--version', action='version', version='0.1')
    args = vars(parser.parse_args())
    DEBUG = args['DEBUG']
    cpu = args["cpu"]
    # print args
    if DEBUG:
        print >>sys.stderr, "Degbugging mode is on."
        print >>sys.stderr, "A file: " + str(args['a'])
        print >>sys.stderr, "B file: " + str(args['b'])

    rang = args["rang"]

    # When index is 0, there is just one nt overlap
    hist = {}
    for i in range(-rang, rang+1):
        hist[i] = 0.0

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
    if DEBUG:
        for i in cl.keys():
            print >>sys.stderr, "DEBUG: A file" + str(sa[i].shape)
            print >>sys.stderr, "DEBUG: B file" + str(sb[i].shape)
    print >>sys.stderr, "Chromosome lengths have been loaded."

    # Read the bed2 file
    bed2_reader(args['a'], sa)
    bed2_reader(args['b'], sb)
    print >>sys.stderr, "Signals have been loaded."
    if DEBUG:
        print >>sys.stderr, "sa"
        print >>sys.stderr, sa
        print >>sys.stderr, "sb"
        print >>sys.stderr, sb
    q = Queue()
    threads = []

    print >>sys.stderr, "Starting calculating signals..."
    # chroms = ("chr2LHet", "chr2RHet")
    # chroms = ("chr2LHet", )    
    chroms = cl.keys()

    def d2_worker_handler(i):
        print >>sys.stderr, "Handling chromosome: " + i
        return d2(sa[i], sb[i], rang)
    print >>sys.stderr, "Number of threads: " + str(cpu)

    p = Pool(cpu)
    res = p.map(d2_worker_handler, chroms)
    p.close()
    p.join()
    for hist_i in res:
        if DEBUG:
            print >>sys.stderr, hist_i
        for j in range(-rang, rang+1):
            hist[j] += hist_i[j]
    # for i in ("chr2LHet",):
    # # for i in cl.keys():
    #     p = Process(target=d2_worker, args=(sa[i], sb[i], rang, q))
    #     p.start()
    #     threads.append(p)
    # for t in threads:
    #     t.join()
    # print "Queue size: " + str(q.qsize())
    # while(q.empty() == False):
    #     r = q.get()
    #     print r
    #     for i in range(-rang, rang+1):
    #         hist[i] += r[i]
    # print hist
    for i in range(-rang, rang+1):
        print str(i) + "\t" + str(hist[i])
