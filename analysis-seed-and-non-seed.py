#!/usr/bin/env python

import pylev
from Bio.Seq import Seq
from Bio import pairwise2
import sys

# first param is SRA fasta
# second param is DEG tab

DEBUG = False
# Human readable
# FORMAT = "r"
# table
FORMAT = "t"
# seq = Seq("TCGGGCCC")
# print seq.reverse_complement()

# Find all SRA-DEG read pairs such that g2-10 are perfect matches

# $ head SRA/Zamore.SRA.pi2H_R7.48dpp.pi17_reads.5plus.fa
# >chr17_27355053_27355076_5_1_+_TGAGGTTTCAAGGGTTGCTAGGG
# TGAGGTTTCAAGGGTTGCTAGGG
# >chr17_27361580_27361605_5_1_+_TAGCTGTCACATTTCTCCCTCTGCT
# TAGCTGTCACATTTCTCCCTCTGCT

fasta_fn = "SRA/Zamore.SRA.pi2H_R7.48dpp.pi17_reads.5plus.fa"
fasta_fn = sys.argv[1]
def rev_comp(x):
    """Just a wrapper around Bio.Seq.reverse_complement()
"""
    return str(Seq(x).reverse_complement())
# def rev_comp(x):
#     for i in range(len(x)):
#         tmp = x[i]
#         if tmp == "A":
#             x[i] = "T"
#         elif tmp == "C":
#             x[i] = "G"
#         elif tmp == "G":
#             x[i] = "C"
#         elif tmp == "T":
#             x[i] = "A"
#         elif tmp == "N":
#             x[i] = "N"
#     return x
def show_alignment(align1, align2, score, begin, end):
## """format_alignment(align1, align2, score, begin, end) -> string 
## Format the alignment prettily into a string. """ 
    s = []
    mat_str = ""
    for i in range(len(align1)):
        if align1[i] == align2[i]:
            mat_str += "|"
        else:
            mat_str += " "
    s.append("%s\n" % align1)
    s.append(mat_str + "\n")
    s.append("%s\n" % align2)
    s.append("  Score=%g\n" % score)
    return ''.join(s) 
class Bed_record:
    """Stores info in one bed record
"""
    def __init__(self, chrom, start, end, copy, ntm, strand, seq):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.copy = copy
        self.ntm = ntm
        self.strand = strand
        self.seq = seq
    def __repr__(self):
        my = ("bed_record", self.chrom, str(self.start), str(self.end), str(self.copy), str(self.ntm), str(self.strand), str(self.seq))
        return "\t".join(my)
    def output_compact(self):
        my = (self.chrom, str(self.start), str(self.end), str(self.copy), str(self.ntm), str(self.strand), str(self.seq))
        return "_".join(my)
# consider g2-g10
# [1, 10)
g_start = 1
g_end = 10
seeds = {}
with open(fasta_fn) as fh:
    for line in fh.readlines():
        line = line.strip()
        if line[0] == ">":
            header = line[1:len(line)]
            col = header.split("_")
            record = Bed_record(col[0], int(col[1]), int(col[2]), int(col[3]), int(col[4]), col[5], col[6])
        else:
            # Only consider single line fasta
            seq = line
            # print record
            myseed = record.seq[g_start: g_end]
            if myseed in seeds:
                seeds[ myseed ].append(record)
            else:
                seeds[ myseed ] = [record, ]
            
def print_seeds(seeds):
    for i in seeds:
        print "seed:" + i
        print "n_species:" + str( len(seeds[i]) )
        for j in seeds[i]:
            print j
        print 
# print_seeds(seeds)

# SRA              5' -------->--------->---------- 3'
# DEG 3' -------<---------------==============<========== 5'
# Cleavage                      |  
#

class Deg_record:
    """Stores infor in one deg record
"""
    def __init__(self, deg_id, deg_signal, deg_up, deg_site, deg_down):
        self.deg_id = deg_id
        self.deg_signal = deg_signal
        self.deg_up = deg_up
        self.deg_site = deg_site
        self.deg_down = deg_down
    def __repr__(self):
        my = ("deg_record", self.deg_id, str(self.deg_signal), self.deg_up, self.deg_site, self.deg_down)
        return "\t".join(my)
    def output_compact(self):
        my = ("degrecord", self.deg_id, str(self.deg_signal), self.deg_up, self.deg_site, self.deg_down)
        return "_".join(my)
    
tab_fn = "WT.putative_sites.x_piGene.l100_1_r100.tab"
tab_fn = sys.argv[2]
# deg reads: 2-10 or [1, 10)
# deg_start = 1
# deg_end = 10

deg_seeds = {}
with open(tab_fn) as fh:
    for line in fh.readlines():
        line = line.strip()
        col = line.split()
        deg_name, deg_sig, deg_up, deg_site, deg_down = col
        deg_sig = float(deg_sig)
        # print deg_name, deg_sig, deg_up, deg_site, deg_down
        tmp = deg_site + deg_down
        deg_seed = tmp[0: 9]
        deg_record = Deg_record(deg_name, deg_sig, deg_up, deg_site, deg_down)
        if deg_seed in deg_seeds:
            deg_seeds[ deg_seed ].append(deg_record)
        else:
            deg_seeds[ deg_seed ] = [deg_record, ]

# print_seeds(deg_seeds)

for i in seeds:
    rc_seed = rev_comp(i)
    if DEBUG:
        print "SRA seed: " + i
        print "rev_comp: " + rc_seed
    # See if the same seed (rc'ed) exists in the DEG library
    if rc_seed in deg_seeds:
        for j in deg_seeds[rc_seed]:
            a = j.deg_up[-30:]
            # b = (j.deg_site + j.deg_down)[deg_start: deg_end]
            # TODO replace 10 with a variable name
            b = (j.deg_site + j.deg_down)[0: 10]
            if DEBUG:
                print "#" + str(j)
                print "DEG read: " + a + " " + b
                print rev_comp(a+b)
                print "SRA"
            # Given a seed, for every seed-matching pair of DEG and SRA, do pairwise alignment 
            for k in seeds[i]:
                if DEBUG:
                    print k.seq
                l = len(k.seq)
                # s1: seq from around DEG sites
                ab = a + b
                s1 = rev_comp(ab)[0: l]
                # s2: seq from SRA (excluding the first pos, i.e. the 1st position does not matter)
                s2 = k.seq
                # print "s1: " + s1
                # print "s2: " + s2
                ed = pylev.levenshtein(s1[10:], s2[10:])
                if ed <= 100:
                    if DEBUG:
                        print "DEG:" + a + b
                        print "s1 from DEG: " + s1[0] + "|" + s1[1:10] + "|" + s1[10:]
                        print "s2 from SRA: " + s2[0] + "|" + s2[1:10] + "|" + s2[10:]
                        print "ed_x_pos1: " + str(ed)
                    # Only do alignment for the rest of the seq (i.e. ignore the 1st position and the seed region ( 2-10, or [1, 10) )
                    # since they are supposed to be perfectly matched
                    # alignments = pairwise2.align.globalxx(s1[10:], s2[10:])
                    # emphasize g10-g21
                    alignments = pairwise2.align.globalxx(s1[10:], s2[10:])
                    if len(alignments) >= 1:
                        aa = alignments[0]
                        # "('GGTTTTCTA-CTTTC----TA-', '--C-T-C-AGCAC-CAGGGTGG', 6.0, 0, 22)"
                        score = aa[2]
                        if score >= 1:
                            if FORMAT == "h":
                                print "DEG:" + a + b
                                print "s1: " + s1[0] + "|" + s1[1:10] + "|" + s1[10:]
                                print "s2: " + s2[0] + "|" + s2[1:10] + "|" + s2[10:]
                                print "ed_x_pos1: " + str(ed)
                                print(show_alignment(*aa))
                                print ""
                            elif FORMAT == "t":
                                # s1(DEG), s2(SRA), ed(x_1st_and_seed), s1_aln, s2_aln, ed, aln_score
                                # print show_alignment(*aa)
                                print "\t".join([ s1[0] + "|" + s1[1:10] + "|" + s1[10:], s2[0] + "|" + s2[1:10] + "|" + s2[10:], aa[0], aa[1], str(ed), str(aa[2]) ]), 
                                print j.output_compact(),
                                print k.output_compact(),
                                print 
    if DEBUG:
        print "-" * 72
