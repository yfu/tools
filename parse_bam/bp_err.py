#!/usr/bin/env python

# Part of BP discovery pipeline
# Given a bam file (reads are mapped to putative lariats), this script figures out the
# the mismatch at the branchpoint position.

import pysam
import sys
import argparse

DEBUG=False
parser = argparse.ArgumentParser(description='A script to collect stats on the errors (SNPs) on lariat-mapping reads')
parser.add_argument('-f', '--file', help='the input file', required=True)
parser.add_argument('-r', '--reference', help='the referece file (indexed fasta file)', required=True)
# parser.add_argument('-o', '--output', help='the output file', required=True)
args = parser.parse_args()
f = args.file
fa_fn = args.reference
fa = pysam.FastaFile(fa_fn)

samfile = pysam.AlignmentFile(f, "rb" )
chroms = samfile.references

# pc: pileupcolumn
print "\t". join( ("#chrom", "pos", "ref", "A", "C", "G", "T", "N") )
for pc in samfile.pileup():
    # print ("\ncoverage at base %s = %s" % (pc.pos, pc.n))
    rpos = pc.reference_pos
    c = chroms[pc.reference_id]
    # print "Inspecting lariat %s at position %d" % (c, rpos)
    d = {}
    for i in ("A", "C", "G", "T", "N"):
        d[i] = 0
    for pileupread in pc.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:  # query position is None if is_del or is_refskip is set.
            qpos = pileupread.query_position
            read_base = pileupread.alignment.query_sequence[qpos]
            # print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))
            # print "Read %s at read position %d with read base %s" % (pileupread.alignment.query_name, qpos, read_base)
            qn = pileupread.alignment.qname
            if DEBUG:
                print "#", qpos, qn, "#"
            d[ read_base.upper() ] += 1
    ref_base = fa.fetch(c, rpos, rpos+1)
    print "\t" . join((c, str(rpos+1), ref_base, str(d["A"]), str(d["C"]), str(d["G"]), str(d["T"]), str(d["N"])))
samfile.close()
                            
