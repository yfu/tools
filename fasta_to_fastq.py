#!/usr/bin/env python

import sys
from Bio import SeqIO
handle = open(sys.argv[1], "rU")
for record in SeqIO.parse(handle, "fasta") :
        n = str(record.id)
        s = str(record.seq)
        print "@" + n
        print s
        print "+"
        print "I" * len(s)
