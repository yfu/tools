#!/usr/bin/env python

# Given a annotation file in GTF format and a bed file from a 5' RACE data, this script generates the normalized count of 5' RACE tags on each gene
# Note that this script only considers those mapping to one feature. For example, if an alignment matches two or more genes, this read is ignored
# Test: python count_reads_for_gene_ntm.py mm9.genes.chr1.gtf 1
# Example: python count_reads_for_5race.py TODO TODO
# Author: Yu Fu

import HTSeq
import gzip
import sys
import itertools

annotation = sys.argv[1]
bed2 = sys.argv[2]
if (len(sys.argv)>3):
    feature_type = sys.argv[3]
else:
    feature_type = "exon"

gtf_file = HTSeq.GFF_Reader( annotation, end_included=True )
# ga = HTSeq.GenomicArrayOfSets("auto", stranded=True)
features = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
counts = {}

# This script will also count the number of reads falling on the antisense strand
counts_as = {}
# for feature in itertools.islice( gtf_file, 10 ):
for feature in gtf_file:
    if feature.type == feature_type:
        # feature name is the first attribute in the last column of a GTF file
        features[ feature.iv ] += feature.name
        counts[ feature.name ] = 0
        counts_as[ feature.name ] = 0

fh = open(bed2)
for line in fh:
    line = line.strip()
    col = line.split()
    chrom, start, end, copy, ntm, strand, seq = col
    start = int(start)
    end = int(end)
    copy = int(copy)
    ntm = int(ntm)

    iv = HTSeq.GenomicInterval( chrom, start, end, strand )
    strand_as = "."
    if strand == "+":
        strand_as = "-"
    elif strand == "-":
        strand_as = "+"
    else:
        print >>sys.stderr, "Found one entry without strand infomation. Skipping..."
        continue
    iv_as = HTSeq.GenomicInterval( chrom, start, end, strand_as )
    iset = None
    for iv, step_set in features[iv].steps():
        # print iv, step_set
        if iset == None:
            iset = step_set.copy()
        else:
            # Strict: this requires that the read maps to one gene. If it ever maps to two or more genes,
            # the length of the set will be longer than 2 and will be dropped.
            # See http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html?highlight=intersection_update#counting-reads-by-genes for details
            iset.intersection_update(step_set)
            # iset.update(step_set)
    l = len(iset)
    if l == 1:    
        for i in list(iset):
            counts[ i ] += float(copy) / ntm / l
    # If the read overlaps with certain features, do not flip the strand
    if (l > 0):
        continue
    iset = None
    for iv, step_set in features[iv_as].steps():
        # print iv, step_set
        if iset == None:
            iset = step_set.copy()
        else:
            iset.intersection_update(step_set)
    l = len(iset)
    if l == 1:
        for i in list(iset):
            counts_as[ i ] += float(copy) / ntm / l
for g in counts:
    print g + "\t" + str(counts[g]) + "\t" + str(counts_as[g])
    
