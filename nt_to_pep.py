#!/usr/bin/env python

from Bio import SeqIO

min_pro_len = 10

print "#name\tframe\tstrand\torf_pep"
with open("test.fa") as fh:
    for record in SeqIO.parse(fh, "fasta"):
        for strand, nuc in [("+", record.seq), ("-", record.seq.reverse_complement())]:
            for frame in range(3):
                length = 3 * ((len(record) - frame) // 3)
                print "%s\t%d\t%s\t%s" % (record.id, frame, strand, nuc[frame: frame+length].translate())
                # for pro in nuc[frame: frame+length].translate().split("*"):
                #     if len(pro) >= min_pro_len:
                #         print "%s\t%d\t%s\t%d\t%s" % (record.id, len(pro), strand, frame, pro)
        
