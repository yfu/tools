#!/usr/bin/env python

# Given a annotation file in bed format and a bed2 format, this script outputs the data for the meta analysis
# Note that some records in the bed file have the same names. For these, I only keep the first entry and throw away the rest
# Author: Yu Fu

import HTSeq
import gzip
import sys

# Normalize features to the length of 100
meta_size = 100
meta = {}
for i in range(meta_size):
    meta[i] = 0.0
annotation = sys.argv[1]
bed2 = sys.argv[2]

if annotation[-3:] == ".gz" or annotation[-5:] == ".gzip":
    f = gzip.open(annotation)
else:
    f = open(annotation)

ga = HTSeq.GenomicArrayOfSets("auto", stranded=True)
feature_dict = {}
for line in f.readlines():
    line = line.strip()
    col = line.split()
    if (len(col) > 6):
        col = col[0:6]
    chrom, start, end, name, score, strand = col
    start = int(start)
    end = int(end)
    iv = HTSeq.GenomicInterval(chrom, start, end, strand)
    # k = name + "," + chrom + ":" + str(start) + "-" + str(end) + strand
    k = name
    if k in feature_dict:
        print >> sys.stderr, "Warning: Found features with the same name for this line. I will ignore this record: \"" + line + "\""
    else:
        feature_dict[k] = iv
        ga[ iv ] += name

# for c in ga.chrom_vectors:
#     print c + "\t" + str(ga.chrom_vectors[c])

f.close()
f = open(bed2)
for line in f.readlines():
    line = line.strip()
    col = line.split()
    # print col
    chrom, start, end, copy, ntm, strand, seq = col
    start = int(start)
    end = int(end)
    copy = int(copy)
    ntm = int(ntm)
    if strand == "+":
        start = start
        end = start + 1
        iv = HTSeq.GenomicInterval(chrom, start, start+1, strand)
    elif strand == "-":
        start = end - 1
        end = end
        iv = HTSeq.GenomicInterval(chrom, end-1, end, strand)
    else:
        print >> sys.stderr, "Error: Found one alignment without strand information"
    for iv, value in ga[iv].steps():
        # One particular alignment may map to multiple overlapping features: ntm_ntm
        features = value
        ntm_ntm = len(features)
        for f in features:
            iv = feature_dict[f]
            f_start = iv.start
            f_end = iv.end
            f_strand = iv.strand
            prop = float(start - f_start) / (f_end - f_start)
            if f_strand == "-":
                prop = 1 - prop
            i = int(prop * meta_size)
            meta[i] += float(copy) / ntm / ntm_ntm
            # print float(copy) / ntm / ntm_ntm
for i in range(meta_size):
    print str(i) + "\t" + str(meta[i])
