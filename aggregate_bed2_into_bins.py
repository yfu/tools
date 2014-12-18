#!/usr/bin/env python

# Given a annotation file in bed format and a bed2 format, this script outputs the data for the meta analysis

import HTSeq
import gzip
import sys

annotation = sys.argv[1]
bed2 = sys.argv[2]

if annotation[-3:] == ".gz" or annotation[-5:] == ".gzip":
    f = gzip.open(annotation)
else:
    f = open(annotation)
    
for line in f.readlines():
    line = line.strip()
    col = line.split()
    
