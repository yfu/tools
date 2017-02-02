#!/usr/bin/env python

# For each siRNA 5' ends, calc the read length distribution. If
# reads between 20-22nt is more than 66%, consider it as a siRNA loci

# stdout contains the good siRNA reads
# stderr contains the bad ones (the ones where piRNA reads dominate)
# Usage: python pick_sirna.py 23-30nt.piRNAs.bed2 20-22nt.siRNAs.bed2

import sys
import collections

CUTOFF = 2.0 / 3

fn1 = sys.argv[1]
fn2 = sys.argv[2]


def bed2key(chrom, start, end, strand):
    if(strand == "+"):
        fivep = start
    else:
        fivep = end - 1
    return chrom + ":" + str(fivep) + ":" + strand
    

def store_counts(fn):
    # Return a dict: key is the five prime position and value is the read count
    # multimappers are discarded
    ret = collections.defaultdict(float)
    with open(fn) as fh:
        i = 0
        for line in fh:
            i += 1
            # if i % 1000000 == 0:
            #     print "Processed %d lines..." % i
            col = line.strip().split()
            chrom, start, end, copy, ntm, strand, seq = col[0:7]
            start = int(start)
            end = int(end)
            copy = float(copy)
            ntm = float(ntm)
            if ntm != 1:
                continue

            k = bed2key(chrom, start, end, strand)
            ret[k] += float(copy) / ntm
    return ret

# position to piRNA read counts
pos2pi = store_counts(fn1)
# position to siRNA read counts
pos2si = store_counts(fn2)

# pos2ratio = {}
# for i in pos2si:
#     si_over_total = pos2si[i] / (pos2si[i] + pos2pi[i])
#     pos2ratio[i] = si_over_total

with open(fn2) as fh:
    i = 0
    for line in fh:
        i += 1
        col = line.strip().split()
        chrom, start, end, copy, ntm, strand, seq = col[0:7]
        start = int(start)
        end = int(end)
        copy = float(copy)
        ntm = float(ntm)
        if ntm != 1:
            continue

        if(strand == "+"):
            fivep = start
        else:
            fivep = end - 1
        
        k = bed2key(chrom, start, end, strand)
        total_pi = 0.0
        total_si = 0.0
        for i in range(-5, 6):
            total_pi += pos2pi[bed2key(chrom, fivep + i, fivep+i+1, strand)]
            total_si += pos2si[bed2key(chrom, fivep + i, fivep+i+1, strand)]

        ratio = total_si / (total_si + total_pi)
        if ratio > CUTOFF:
            print line.strip() + "\t"+ str(ratio)
        else:
            print >>sys.stderr, line.strip() + "\t" + str(ratio)
