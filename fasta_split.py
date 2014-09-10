#!/usr/bin/env python

# Split a fasta into multiple files. Each entry has its own file

from Bio import SeqIO
import sys

fh = open(sys.argv[1])
for record in SeqIO.parse(fh, "fasta"):
    out = open(str(record.id) + ".fa", "w")
    print >> out, ">" + str(record.id)
    print >> out, str(record.seq)
    out.close()

fh.close()
