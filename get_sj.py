#!/usr/bin/env python

# Author: Yu Fu

# Get the splice junction positions and the # reads supporting each of the junctions.
# Usage: cat reads.cluster.pairs.support.junction.noheader.sam | python run5-get-sj.py > sj_depth.bed
# Then do
# bedtools intersect -a sj_depth.bed -b ~/data/piPipes/common/dm3/Brennecke.piRNAcluster.bed6.gz -wo | awk 'BEGIN{OFS="\t"} { $4=$9; $6="."; print}' | cut -f1-6 | sort -k1,1 -k2,2n > sj_depth_cluster.bed 

import fileinput
import re

pat_cigar = re.compile(r"[0-9]+[MIDNSHP=X]")
pat_one_cigar = re.compile(r"([0-9]+)([MIDNSHP=X])")

def sj_pos(pos, cigar):
    # Given a (1-based) start position and a cigar string. Return all the splice junction positions
    # inferred from this read
    sj = []
    # print pos, cigar
    mat = pat_cigar.findall(cigar)
    cur_pos = int(pos)
    for m in mat:
        # print m
        mm = pat_one_cigar.match(m)
        num = int(mm.group(1))
        s = mm.group(2)
        if s == "M":
            cur_pos += num - 1
        elif s == "N":
            sj.append(cur_pos)
            cur_pos += num + 1
            sj.append(cur_pos)
        elif s == "S":
            pass
        elif s == "D":
            cur_pos += num + 1
        # When the cigar is "I", the position on the genome does not move
        elif s == "I":
            pass
        elif s == "H":
            pass
        # P, = and X do not appear in the sam file
        elif s == "P":
            pass
        elif s == "=":
            pass
        elif s == "X":
            pass
        else:
            print >> sys.stderr, "I found a CIGAR that I cannot recognize:", num, s
            
        # print sj
        
    # print "="*32
    return sj

if __name__ == "__main__":
    all_sj = {}
    chrom = ""
    for line in fileinput.input():
        line = line.strip()
        col = line.split()
        chrom = col[2]
        pos = col[3]
        cigar = col[5]
        sj = sj_pos(pos, cigar)
        for i in sj:
            # i is a 1-based coordinate
            # k is "chr2L	100	101", which is bed coordinates
            k = chrom + "\t" + str(i-1) + "\t" + str(i)
            if k in all_sj:
                all_sj[k] += 1
            else:
                all_sj[k] = 1
    for k in all_sj.keys():
        print k +"\t" + "." + "\t" + str(all_sj[k])
