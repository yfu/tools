#!/usr/bin/env python

# Author: Yu Fu (yfu@yfu.me)

from Bio import SeqIO
from itertools import *
import sys

def find_n(seq, n_limit):
    ''' Given a DNA string and n_limit, find all the start and end positions
    for n_limit or more consecutive N's
    '''
    ret = []
    n_flag = False
    for i in range(len(seq)):
        if seq[i] == 'N' or seq[i] == 'n':
            if n_flag == False:
                start = i
                end = i + 1
                n_flag = True
            else:
                # end marks (one past) the end position
                end = i + 1
        else:
            if n_flag == True:
                if (end - start >= n_limit):
                    ret.append([start, end])
                    # print str(start) + ":" + str(end)
                n_flag = False
    return ret

def scaf_to_contig(seq, n_limit):
    ''' Given the name and seq of a scaffold, this function will break
    the scaffold into a list of contigs, if consecutive 25 (n_limit) nucleotides
    are N's. '''
    contigs = []
    coor = find_n(seq, n_limit)          # coor looks like this: [[start1, end1], [start2, end2],...]
    prev = 0
    for c in coor:
        contigs.append(seq[prev:c[0]])
        prev = c[1]
    if len(coor) != 0:
        contigs.append(seq[prev:])
    # print "contigs list" + str(contigs)
    return contigs
    

def main():
    # debug
    # print scaf_to_contig('ATGC' + 'N' * 20 + 'ATGC' + 'N'*25 + 'ATGC' + 'N' * 30 + 'ATGC' + 'N' * 5 + 'ATGC', 25)
    # sys.exit("debug")
    fh = open(sys.argv[1])
    n_limit = 25
    # The total length of all scaffolds with N's included
    total_scaf_len = 0
    total_contig_len = 0
    # A list containing the length of every scaffold
    scaf_len = []
    contig_len = []

    # The total length of all N's
    total_scaf_gap = 0
    total_contig_gap = 0

    for record in SeqIO.parse(fh, "fasta"):
        seq = str(record.seq)
        one_scaf_len = len(seq)
        total_scaf_len += one_scaf_len
        scaf_len.append(one_scaf_len)
        total_scaf_gap += seq.count('N') + seq.count('n')
        contigs = scaf_to_contig(seq, n_limit)
        for c in contigs:
            l = len(c)
            total_contig_len += l
            contig_len.append(l)
            total_contig_gap += seq.count('N') + seq.count('n')
        
    print "Total scaffold length" + "\t" + str(total_scaf_len)
    print "Total length of gaps" + "\t" + str(total_scaf_gap)

    avg_scaf_len = sum(scaf_len) / len(scaf_len)
    print "The average scaffold length" + "\t" + str(avg_scaf_len)
    scaf_len.sort(reverse=True)

    # Accumulted scaffold length, which is used to calculate N10, N20...
    scaf_len_accum = 0
    nxx_cutoff = map(lambda x: float(x) / 10 * total_scaf_len, range(1, 11)) # Include N100...
    print nxx_cutoff
    nxx_flag = 1

    for i in scaf_len:
        scaf_len_accum += i
        if scaf_len_accum >= nxx_cutoff[nxx_flag - 1]:
            nxx_flag += 1
            print "N" + str(nxx_flag - 1) + '0' + "\t" + str(i)
    
    print "Total contig length" + "\t"  + str(total_contig_len)
    print "Total lengths of gaps in contigs" + "\t"  + str(total_contig_gap)
if __name__ == '__main__':
    main()

