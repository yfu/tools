#!/usr/bin/env python

# Given two bed files (sorted by names instead of coordinate) with the first
# one being the alignments of the "before" parts
# and the second one being the alignments of the "align+after" parts,
# this script outputs a list of intron mapping reads
# It basically "joins" the alignments of "before" and "align+after" part

import sys

bef_fn = sys.argv[1]
bef = open(bef_fn, "r")
# align + after
alaf_fn = sys.argv[2]
alaf = open(alaf_fn, "r")

class Segment:
    def __init__(self, a):
        # input is an entry in a bed file
        self.chrom, self.start, self.end, self.name, self.score, self.strand = a
        self.start = int(self.start)
        self.end = int(self.end)
    def __str__(self):
        return "\t".join([self.chrom, str(self.start), str(self.end), self.strand])
    def output_all(self):
        return "\t".join([self.chrom, str(self.start), str(self.end), self.name, self.score, self.strand])

finished = False
def parse_line(fh):
    t = fh.readline()
    if t == '':
        global finished
        finished = True
    else:
        return Segment(t.split())
print parse_line(bef)

# print parse_line(bef)
a = parse_line(bef)
b = parse_line(alaf)

while True:
    print "a: " + str(a) + a.name
    print "b: " + str(b) + b.name
    if finished:
        print >> sys.stderr, "Finished!"
        sys.exit(0)
    if a.name == b.name:
        print "Find one concordant read: " + str(a) + \
          "\t" + str(b)
        a = parse_line(bef)
        b = parse_line(alaf)
    elif a.name > b.name:
        b = parse_line(alaf)
    else:
        a = parse_line(bef)
        
                


            
