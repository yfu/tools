#!/usr/bin/env python

from gppc_core import d1, d1_worker, d2, d2_worker
import sys
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
from multiprocessing import Process, Queue

if (len(sys.argv)<5):
    print >>sys.stderr, "Usage: gppc1.py asm.ChromInfo my.bed2 [55|33|53] [ss|ds|all]"
    print >>sys.stderr, "asm.ChromInfo contains the chromosome length, 55|33|53: 5' to 5', 3' to 3', 5' to 3' distance, ss|ds: consider reads on the same strand, different strand"
    print >>sys.stderr, "If you ever want to calculate distances considering reads on the same strand and the different strand, consider run ss and ds and then add them together"    
    exit(1)
rang =50
# When index is 0, there is just one nt overlap
hist = {}
for i in range(-rang, rang+1):
    hist[i] = 0.0

if __name__ == "__main__":
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
        if (counter % 1000000 == 0):
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
    q = Queue()

    ori = sys.argv[4]
    if sys.argv[3] == "55":
        s = sig5
    elif sys.argv[3] == "33":
        s = sig3
    elif sys.argv[3] == "53":
        s5 = sig5
        s3 = sig3
    threads = []

    if sys.argv[3] == "55" or sys.argv[3] == "33":
        for i in ("TAS",):        
            p = Process(target=d1_worker, args=(s5[i], s3[i], rang, ori, q))
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
    elif sys.argv[3] == "53":
        for i in ("TAS",):        
            p = Process(target=d2_worker, args=(s5[i], s3[i], rang, ori, q))
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



    


