#!/usr/bin/env python

## Assign multimappers to the windows using expectation-maximization
## Step 0: use unique mappers to determine the initial state
## Step E: assign multimappers to windows based on probabilies linearly scaled by #reads in windows
## Step M: adjust the #reads in each window
## Repeat Step # and Step M till convergence

import gzip
import sys

## Stores all the reads: sequence as the key, list of winid's from all mapping positions form the value
read_winid = {}
## Stores how many times we see the read in the library
read_count = {}
## How many reads are assigned to different windows
win_count = {}


## Intersection of bed2 and 5k windows
input_fn = sys.argv[1]
output_fn = sys.argv[2]

with gzip.open(input_fn) as fh:
    for line in fh:
        line = line.strip()
        ## print line
        col = line.split()
        winid = "\t".join([col[7], col[8], col[9], col[12]])
        seq = col[6]
        ## print winid
        win_count[winid] = 1        
        read_count[seq] = int(col[3])
        if seq in read_winid:
            read_winid[seq].append(winid)
        else:
            read_winid[seq] = [winid, ]
            
total_n_read = 0.0
for k in read_count:
    total_n_read += read_count[k]

def iteration(win_count, read_winid):
    new_win_count = {}
    for win in win_count:
        new_win_count[win] = 0
        
    for read in read_winid:
        ids = read_winid[read]
        tmp = [win_count[id] for id in ids]
        prob = [float(i) / sum(tmp) for i in tmp]
        ## print prob
        for p in zip(ids, prob):
            new_win_count[p[0]] += read_count[read] * p[1]

    return new_win_count


def print_top_win(w):
    max_line = 10
    c = 0
    for i in sorted(w, key=w.get, reverse=True):
        c += 1
        if(c == max_line):
            break
        print >>sys.stderr, w[i], i


def output_final_count(w, fn):
    with open(fn, "w") as fh:
        for i in w:
            print >>fh, i + "\t" + str(w[i])

            
def print_total_n_read(w):
    s = 0.0
    for i in w:
        s += w[i]
    print >>sys.stderr, "Total number of reads: %f" % s


def delta(w1, w2):
    '''It calculates the difference between two rounds of results
'''
    ret = 0.0
    for i in w1:
        ret += abs(w1[i] - w2[i])
    return ret

    
max_iter = 100
n_iter = 0
# epsilon = 0.001
# To make it finish earlier
epsilon = 0.05

while(True):
    n_iter += 1
    if(n_iter > max_iter):
        print >>sys.stderr, "Reaching maximum number of iterations: %d" % (max_iter)
        output_final_count(new_win_count, output_fn)
        break
    new_win_count = iteration(win_count, read_winid)
    print >>sys.stderr, "Iteration %d" % n_iter
    print_top_win(win_count)
    print >>sys.stderr, "#" * 72
    print_top_win(new_win_count)
    print_total_n_read(new_win_count)
    diff = delta(win_count, new_win_count)
    print >>sys.stderr, "Delta between two iterations: %f" % diff
    if diff <= total_n_read * epsilon:
        print >>sys.stderr, "Reached convergence!"
        output_final_count(new_win_count, output_fn)
        break
    win_count = new_win_count
    print >>sys.stderr, ""



