#!/usr/bin/env python

# A simple Ping-Pong score calculator
# Given a sorted bed2 file, this script will output the nucleotide offset in
# the first column and number of read pairs having such an offset in 
# the second column
#
# Test data:
# cat ~/data/prepachytene/results/2014-11-21/SRA_mm10g/00dpp/genome_mapping/Zamore.SRA.total.00dpp.testis.trimmed.x_rRNA.x_hairpin.mm10gv1.all.bed2 | head -n1000000 > test.bed2
#
# Usage: pps_simple.py ~/data/shared/mm10/mm10.chrom.sizes ../2014-11-21/SRA_mm10g/00dpp/genome_mapping/Zamore.SRA.total.00dpp.testis.trimmed.x_rRNA.x_hairpin.mm10gv1.all.bed2

import HTSeq
import sys

gl_fn = sys.argv[1]
gl = {}
# This script will look for -50 to +50 (inclusive, so it is 101 in total) of every piRNA for its PP partner
# The two variables have to be the same
left = right = 50

hist = {}
for i in range(-left, right+1):
    hist[i] = 0
    
with open(gl_fn) as fh:
    for line in fh.xreadlines():
        line = line.strip()
        col = line.split()
        gl[col[0]] = int(col[1])
print >>sys.stderr, "Chromosome lengths have been loaded"
# for c in gl:
#     print c + "\t" + gl[c]
    
fn = sys.argv[2]
b = HTSeq.BED_Reader(fn)
# Use ndarray since the data will be constantly changing, unlike genome annotations
ga = HTSeq.GenomicArray(chroms=gl, stranded=True, typecode="d", storage="ndarray")

for aln in b:
    iv = aln.iv
    if iv.strand == "+":
        iv.end = iv.start + 1
    elif iv.strand == "-":
        iv.start = iv.end - 1
    else:
        print >>sys.stderr, "Found an alignment w/o strand information. I will ignore it..."
        continue
    copy = int(aln.name)
    ntm = int(aln.score)
    impact = float(copy) / ntm
    ga[iv] += impact

# This is the second time the same file is read
# This time, for every alignment I get, I will calculate how many
# reads have n nt overlap and store it in the dictionary.

b = HTSeq.BED_Reader(fn)
for aln in b:
    iv = aln.iv
    copy = int(aln.name)
    ntm = int(aln.score)
    impact = float(copy) / ntm
    strand = iv.strand

    # target_strand="dummy"
    if strand == "+":
        target_strand = "-"
    elif strand == "-":
        target_strand = "+"
    else:
        print >>sys.stderr, "Found an alignment w/o strand information. I will ignore it..."
        continue
    five_prime_pos = iv.start_d
    # print "Five primer position: ", five_prime_pos
    # 
    # When one piRNA species has multiple potential partners, assign the weights according to the
    # abundances (impact) of all potential partners. For example, one species have an impact of 5
    # it has two partners: 2nt-overlap species of impact 2 and 3nt-overlap species of impact 4
    # Then 2nt gains (5 * 2/(2+4)) * 2
    # 3nt gains (5*4/(2+4)) * 4
    partners = {}
    for i in range(-left, right+1):
        partners[i] = 0
    partners_sum = 0
    for i in range(-left, right+1):
        target_chrom = iv.chrom
        ## TODO: determine if five_prime_pos + i is <0 or >chromosome length
        target_pos = five_prime_pos + i
        if target_pos < 0 or target_pos >= gl[target_chrom]:
            print >>sys.stderr, "Found one alignment that is too close to the chromosome border. Skipping: " + str(iv)
            continue
        p = HTSeq.GenomicPosition(target_chrom, target_pos, target_strand)
        # print i, p
        # print ga[ p ]
        # Target
        t = ga[ p ]
        if strand == "+":
            partners[i] += t
            partners_sum += t
        elif strand == "-":
            partners[-i] += t
            partners_sum += t
    # Orphan piRNA which does not have Ping-Pong partners
    # print partners_sum
    if partners_sum == 0:
        continue
    for k in partners:
        # print impact * partners[k] / partners_sum * partners[k]
        hist[k] += impact * (partners[k] / partners_sum) * partners[k]

for i in range(-left, right+1):
    print str(i) + "\t" + str(hist[i])
