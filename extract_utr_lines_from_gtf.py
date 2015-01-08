#!/usr/bin/env python
# 
import sys, re
fn = sys.argv[1]
fh = open(fn)

pat = re.compile(r"transcript_id \"([^\"]+)\";")
# pat = re.compile(r"transcript_id")
pat2 = re.compile(r"gene_id \"([^\"]+)\";")
pat3 = re.compile(r"exon_number \"([^\"]+)\";")

transcripts = {}

three_utr = {}
five_utr = {}
n_exons_3utr = {}
n_exons_5utr = {}

class UTR():
    def __init__(self, transcript_id, chrom, start, end, strand, exon_number):
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.exon_number = exon_number
    def __repr__(self):
        return "\t".join([chrom, str(start), str(end), transcript_id, str(exon_number), strand])
        
for line in fh:
    line = line.strip()
    col = line.split("\t")
    chrom = col[0]
    category = col[2]
    start = int(col[3]) - 1
    end = int(col[4])
    strand = col[6]
    info = col[8]
    m = pat.search(info)
    if m is None:
        print >> sys.stderr, "Wrong transcript id!"
        print >> sys.stderr, "Error: malformatted line: " + line + "\n"
        exit(2)
    transcript_id = m.group(1)

    m2 = pat2.search(info)
    if m2 is None:
        print >> sys.stderr, "Wrong gene id!"
        print >> sys.stderr, "Error: malformatted line: " + line + "\n"
        exit(2)
    gene_id = m2.group(1)

    m3 = pat3.search(info)
    if m3 is None:
        print >> sys.stderr, "Wrong gene id!"
        print >> sys.stderr, "Error: malformatted line: " + line + "\n"
        exit(2)
    exon_number = m3.group(1)

    if category == "3UTR":
        if transcript_id not in small_n_exons_3utr:
            small_n_exons_3utr[transcript_id] = exon_number
        else:
            if exon_number < small_n_exons_3utr[transcript_id]:
                small_n_exons_3utr[transcript_id] = exon_number
    if category == "5UTR":
        if transcript_id not in large_n_exons_5utr:
            large_n_exons_5utr[transcript_id] = exon_number
        else:
            if exon_number > large_n_exons_5utr[transcript_id]:
                large_n_exons_5utr[transcript_id] = exon_number

    if category == "3UTR":
        if transcript_id not in three_utr:
            three_utr[transcript_id] = []
        three_utr[transcript_id].append(UTR(transcript_id, chrom, start, end, strand, exon_number))
    elif category == "5UTR":
        if transcript_id not in five_utr:    
            five_utr[transcript_id] = []
        five_utr[transcript_id].append(UTR(transcript_id, chrom, start, end, strand, exon_number))

print three_utr
for u in three_utr:
    if len(three_utr[u]) > 1:
        for i in three_utr[u]:
            if i.exon_number > small_n_exons_3utr[i.transcript_id]:
                print i

