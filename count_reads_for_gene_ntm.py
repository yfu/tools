#!/usr/bin/env python

# Given a annotation file in GTF format and a bed2 format, this script outputs the piRNA counts for each gene
# Test: python count_reads_for_gene_ntm.py mm9.genes.chr1.gtf 1
# Example: python count_reads_for_gene_ntm.py count_reads_for_gene_ntm.mm9.genes.chr1.gtf ~/data/prepachytene/results/2014-11-21/SRA/00dpp/genome_mapping/Zamore.SRA.total.00dpp.testis.trimmed.x_rRNA.x_hairpin.mm9v1.unique.bed2
# Author: Yu Fu

import HTSeq
import gzip
import sys
import itertools

annotation = sys.argv[1]
bed2 = sys.argv[2]

gtf_file = HTSeq.GFF_Reader( annotation, end_included=True )
# ga = HTSeq.GenomicArrayOfSets("auto", stranded=True)
exons = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
counts = {}

# for feature in itertools.islice( gtf_file, 10 ):
for feature in gtf_file:
    if feature.type == "exon":
        # feature name is the first attribute in the last column of a GTF file
        exons[ feature.iv ] += feature.name
        counts[ feature.name ] = 0

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
    iset = None
    for iv, step_set in exons[iv].steps():
        # print iv, step_set
        if iset == None:
            iset = step_set.copy()
        else:
            # Strict: this requires that the read maps to one gene. If it ever maps to two or more genes,
            # the length of the set will be longer than 2 and will be dropped.
            iset.intersection_update(step_set)
            # iset.update(step_set)
    if len(iset) == 1:
        # print copy/ntm
        counts[ list(iset)[0] ] += copy / ntm
for g in counts:
    print g + "\t" + str(counts[g])
    
