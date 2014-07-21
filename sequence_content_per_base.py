#!/usr/bin/env python

import sys

fn = sys.argv[1];

d = [];
expected_len=150
for i in range(expected_len):
    d.append( {'A':0, 'T':0, 'C':0, 'G':0, 'N':0 } )


with open(fn) as f:
    for line in f.readlines():
        line = line.strip();
        e = line.split();
        seq = e[0]
        for i in range(len(seq)):
            d[i][seq[i]] += int(e[1])

print "\t".join([ 'pos', 'A', 'C', 'G', 'T', 'N' ])
for i in range(expected_len):
    print str(i+1) + "\t",
    print str(d[i]['A']) + "\t",
    print str(d[i]['C']) + "\t",
    print str(d[i]['G']) + "\t",
    print str(d[i]['T']) + "\t",
    print str(d[i]['N'])

    
