#!/usr/bin/env python

"""
Join all the counts for transposons and output such a table:

Each row is a sample (barcode) and each column stands for a transposon
"""

import os

f = 0 # File descriptor number
for db in ['repBase', 'transposon']:
    counts = {}
    all_tranposons = set()
    for i in ['HEIL', 'JMAC', 'KNBD', 'OQSU', 'PRTV']:
        for j in ['id1', 'id2' , 'id3', 'id4', 'id5', 'id6']:
            exp = i + "_" + j
            counts[exp] = {}
            in_f = exp + '.' + db + '.ct'
            fh = open(in_f)
            for line in fh.readlines():
                e = line.split()
                all_tranposons.add(e[0])
                counts[exp][e[0]] = e[1]
                
    f += 1
    fd = os.fdopen(f, "w")

    # Print the header (the first row)
    for t in sorted(all_tranposons):
        fd.write("\t" + t)
    fd.write("\n")

    # Now print the data
    # The first table is output via STDOUT and the second table is output via STDERR
    for exp in sorted(counts.keys()):
        fd.write(exp)
        for transposon in sorted(all_tranposons):
            if transposon in counts[exp].keys():
                fd.write("\t" + counts[exp][transposon])
            else:
                fd.write("\t" + "0")
        fd.write("\n")
            
