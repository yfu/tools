#!/usr/bin/env python

# Count those reads that are uniquely assignable to a transposon family

import sys
import collections

fn = sys.argv[1]

read2tpn = {}
read2weight = collections.defaultdict(float)

# Find those reads that can only be assigned to one tpn family
with open(fn) as fh:
    for line in fh:
        col = line.strip().split()
        chrom, start, end, copy, ntm, strand, seq = col[0:7]
        start = int(start)
        end = int(end)
        copy = float(copy)
        ntm = float(ntm)
        # tpn_id, family, class, divergence_rate, unique_id
        tpn_id, tpn_fam, tpn_cla, tpn_div, _ = chrom.split(",")
        # print tpn_id, tpn_fam, tpn_cla, tpn_div
        if seq in read2tpn:
            read2tpn[seq].add(tpn_fam)
        else:
            read2tpn[seq] = set([tpn_fam, ])
        read2weight[seq] += copy / ntm

# Read counts for tpn
tpn_count = collections.defaultdict(float)

for read in read2tpn:
    if len(read2tpn[read]) == 1:
        tpn = read2tpn[read].pop()
        tpn_count[tpn] += read2weight[read]

for tpn in tpn_count:
    print "%s\t%f" % (tpn, tpn_count[tpn])
