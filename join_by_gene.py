#!/usr/bin/env python

"""
Join all the FPKMs for genes and output such a table:

Each row is a sample (barcode) and each column stands for a genes
"""

import os

f = 0 # File descriptor number
for db in ['cufflinks_flyBase_compatible_hits_norm', 'cufflinks_flyBase_total_hits_norm', 'cufflinks_iGenome_compatible_hits_norm', 'cufflinks_iGenome_total_hits_norm']:
    counts = {}
    all_genes = set()
    for i in ['HEIL', 'JMAC', 'KNBD', 'OQSU', 'PRTV']:
        for j in ['id1', 'id2' , 'id3', 'id4', 'id5', 'id6']:
            exp = i + "_" + j
            counts[exp] = {}
            in_f = exp + '_' + db + '.genes.fpkm_tracking'
            fh = open(in_f)
            fh.readline()
            for line in fh.readlines():
                e = line.split()
                # print e[4]
                all_genes.add(e[4])
                counts[exp][e[4]] = e[7]
                
    f += 1
    fd = os.fdopen(f, "w")

    # Print the header (the first row)
    for t in sorted(all_genes):
        fd.write("\t" + t)
    fd.write("\n")

    # Now print the data
    # The first table is output via STDOUT and the second table is output via STDERR
    for exp in sorted(counts.keys()):
        fd.write(exp)
        for gene in sorted(all_genes):
            if gene in counts[exp].keys():
                fd.write("\t" + counts[exp][gene])
            else:
                fd.write("\t" + "0")
        fd.write("\n")
            
