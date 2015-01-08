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
def lonely_utr(utrs, end):
    # Given a set of utrs for one transcript, output those UTRs which occupy whole exons, i.e. exclude the case in which a UTR and a CDS share an exon
    if len(three_utr) <=1:
        print >>sys.stderr("Error: lonely_utr() should be used only for those cases when a gene has multiple 3UTRs or multiple 5UTRs")
        exit(2)
    s = sorted(utrs, key=lambda x: x.start)
    ret = []
    strand = s[0].strand
    if strand == "+":
        if end == "3UTR":
            ii = 0
            for i in range(1, len(s)):
                ret.append(s[i])
        if end == "5UTR":
            for i in range(0, len(s)-1):
                ret.append(s[i])
    else:
        if end == "3UTR":
            for i in range(0, len(s)-1):
                ret.append(s[i])
        if end == "5UTR":
            for i in range(1, len(s)):
                ret.append(s[i])        
    return ret
    
class UTR():
    def __init__(self, transcript_id, chrom, start, end, strand, exon_number):
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.exon_number = exon_number
    def __repr__(self):
        return "\t".join([self.chrom, str(self.start), str(self.end), self.transcript_id + "_" + str(self.exon_number), "0", self.strand])
        
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
        if transcript_id not in n_exons_3utr:
            n_exons_3utr[transcript_id] = exon_number
        else:
            if exon_number < n_exons_3utr[transcript_id]:
                n_exons_3utr[transcript_id] = exon_number
    if category == "5UTR":
        if transcript_id not in n_exons_5utr:
            n_exons_5utr[transcript_id] = exon_number
        else:
            if exon_number > n_exons_5utr[transcript_id]:
                n_exons_5utr[transcript_id] = exon_number

    if category == "3UTR":
        if transcript_id not in three_utr:
            three_utr[transcript_id] = []
        # print UTR(transcript_id, chrom, start, end, strand, exon_number)
        # three_utr[transcript_id].append(UTR(transcript_id, chrom, start, end, strand, exon_number))
        three_utr[transcript_id].append(UTR(transcript_id, chrom, start, end, strand, exon_number))
    # elif category == "5UTR":
    #     if transcript_id not in five_utr:    
    #         five_utr[transcript_id] = []
    #     five_utr[transcript_id].append(UTR(transcript_id, chrom, start, end, strand, exon_number))

for u in three_utr.keys():
    if len(three_utr[u]) > 1:
        l = lonely_utr(three_utr[u], "3UTR")
        for ll in l:
            print ll

