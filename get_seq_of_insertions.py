#!/usr/bin/env python

import fileinput, re
import sys

def pretty_print(d):
    # Pretty-print the content of the dictionary.
    for k in sorted(d.keys()):
        print k + "\t", 
        for i in d[k]:
            print i,
        print

pat = re.compile(r'\d+\w')
pat_num = re.compile(r'^\d+')
pat_char = re.compile(r'\w$')

# This is used to memorize the genomics location we are processing. If the location does not change,
# then we should not output the content. If the location changes, we know we have finished processing
# the current location and we go to the next location
prev_loc_first = 0
prev_loc_last = 0
for line in fileinput.input():
    line = line.strip()
    e = line.split()
    chrom = e[2]
    loc = int(e[3])
    seq = e[9]
    read_len = len(seq)

    # CIGAR string
    # print e[5]
    mat = re.findall(pat, e[5])
    if mat == None:
        continue
    gr = mat
    # print(gr)

    # Ignore those reads like this 49M, and 1S47M1S
    if len(gr) <= 1 or len(gr) >2:
        continue

    # Soft-clipping only occurs at the beginning or the end
    part_first = gr[0]
    # print part_first + "**"
    part_last = gr[-1]

    m = re.findall(pat_num, part_first)
    num = int(m[0])
    m = re.findall(pat_char, part_first)
    char = m[0]

    flag_soft_at_beg = False
    if char == 'S':
        print >>sys.stderr, "I found a soft-clipping at the beginning! \nCurrent CIGAR: " + char + ". And current number of the CIGAR: " + str(num)
        flag_soft_at_beg = True
        # This read is soft-clipped at the beginning
        my_key = chrom + "\t" + str(loc) + "\t" + "U"
        print my_key + "\t" + seq[0:num]

    # print part_last + "**"
    m = re.findall(pat_num, part_last)
    num = int(m[0])
    m = re.findall(pat_char, part_last)
    char = m[0]
    # print char + "**"

    if char == 'S':
        if flag_soft_at_beg == True:
            print >>sys.stderr, "Something might be wrong. This read has soft-clipping at both the beginning and the end"
        print >>sys.stderr, "I found a soft-clipping at the end! \nCurrent CIGAR: " + char + ". And current number of the CIGAR: " + str(num)
        # This read is soft-clipped at the end
        loc = loc + read_len - num
        my_key = chrom + "\t" +  str(loc) + "\t" + "D"
        print my_key + "\t" + seq[read_len-num:]


