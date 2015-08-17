#!/usr/bin/env python

import gffutils
import sys
import numpy

db = gffutils.FeatureDB(sys.argv[1], keep_order=True)
# genes = db.all_features(limit="chr1:1-100000", featuretype="gene")
genes = db.all_features(featuretype="gene")

VERBOSE=False
for g in genes:
    # print g
    l = []
    if len(g["gene_id"]) > 1:
        print >>sys.stderr, "Warning: duplicate gene_id!"
    if VERBOSE == True:
        print "Inspecting gene %s: %s:%d-%d|%s" % ( str(g["gene_id"][0]), g.chrom, g.start, g.end, g.strand ) 
    transcripts = db.children(g, featuretype="transcript")
    for t in transcripts:
        t_len = 0
        if VERBOSE == True:
            print "\tInspecting transcript %s: %s:%d-%d|%s" % ( str(t["transcript_id"][0]), t.chrom, t.start, t.end, t.strand )
        exons = db.children(t, featuretype="exon")
        for e in exons:
            if VERBOSE == True:
                print "\t\tInspecting exon %s: %s:%d-%d|%s" % ( str(e["transcript_id"][0]), e.chrom, e.start, e.end, e.strand )
                print "\t\tExon length: %d" % (e.end - e.start + 1)
            t_len += e.end - e.start + 1
        if VERBOSE == True:
            print "Transcript length: %d" % t_len
        gid = str(g["gene_id"][0])
        tid = str(t["transcript_id"][0])
        print gid + "\t" + tid + "\t" + str(t_len)
        l.append(t_len)
