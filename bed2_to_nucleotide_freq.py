#!/usr/bin/env python

# Given a bed2 file as input and output the position specific matrix.
# This takes NTMs into account
# Author: Yu Fu

import fileinput

# PWM: 4 rows representing ATCG and 50 columns recording 50 positions.
# For piRNA, probably 23-35nt reads are used.
width = 50
m = []
for i in range(width):
    m.append({"A":0, "T":0, "C":0, "G":0})

for line in fileinput.input():
    line = line.strip()
    col = line.split()
    _, _, _, copy, ntm, _, seq = col
    copy = int(copy)
    ntm = int(ntm)
    for i in range(len(seq)):
        m[i][seq[i]] += float(copy) / ntm
print "\t" + "\t".join([str(i) for i in range(width)])
for i in ["A", "T", "C", "G"]:
    print i + "\t",
    for j in range(width):
        print str(m[j][i]) + "\t",
    print ""
    
    
