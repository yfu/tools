#!/usr/bin/env python

# Author: Yu Fu (yfu@yfu.me)

from Bio import SeqIO
from itertools import *
import sys

fh = open(sys.argv[1])

# The total length of all scaffolds with N's included
total_scaf_len = 0
# A list containing the length of every scaffold
scaf_len = []

# The total length of all N's
total_gap = 0
for record in SeqIO.parse(fh, "fasta"):
    seq = str(record.seq)
    one_scaf_len = len(seq)
    total_scaf_len += one_scaf_len
    scaf_len.append(one_scaf_len)
    total_gap += seq.count('N') + seq.count('n')
    

print "Total scaffold length" + "\t" + str(total_scaf_len)
print "Total length of gaps" + "\t" + str(total_gap)

avg_scaf_len = sum(scaf_len) / len(scaf_len)
print "The average scaffold length" + "\t" + str(avg_scaf_len)
scaf_len.sort(reverse=True)

# Accumulted scaffold length, which is used to calculate N10, N20...
scaf_len_accum = 0
nxx_cutoff = map(lambda x: float(x) / 10 * total_scaf_len, range(1, 11)) # Include N100...
print nxx_cutoff
nxx_flag = 1

for i in scaf_len:
    scaf_len_accum += i
    if scaf_len_accum >= nxx_cutoff[nxx_flag - 1]:
        nxx_flag += 1
        print "N" + str(nxx_flag - 1) + '0' + "\t" + str(i)
