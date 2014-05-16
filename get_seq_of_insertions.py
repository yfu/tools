#!/usr/bin/env python

import fileinput, re

# Upstream insertions: it contains the insertions that happen before a genomic location
# The genomic location is the key and the insertions are combined to form a list as the value
u_ins = {}
# Downstream insertions: it contains the insertions that happen after a genomic location
d_ins = {}
pat = re.compile(r'\d+\w')
pat_num = re.compile(r'^\d+')
pat_char = re.compile(r'\w$')
for line in fileinput.input():
    line = line.strip()
    e = line.split()
    chrom = e[2]
    loc = e[3]
    seq = e[9]
    read_len = len(seq)

    # CIGAR string
    # print e[5]
    mat = re.findall(pat, e[5])
    if mat == None:
        continue
    gr = mat
    print(gr)

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
        print "I found a soft-clipping at the beginning! \nCurrent CIGAR: " + char + ". And current number of the CIGAR: " + str(num)
        flag_soft_at_beg = True
        # This read is soft-clipped at the beginning
        my_key = chrom + " " + loc
        if my_key not in u_ins.keys():
            u_ins[my_key] = [ seq[0:num] ]
        else:
            u_ins[my_key].append( seq[0:num] )

    print part_last + "**"
    m = re.findall(pat_num, part_last)
    num = int(m[0])
    m = re.findall(pat_char, part_last)
    char = m[0]
    print char + "**"

    if char == 'S':
        if flag_soft_at_beg == True:
            print "Something might be wrong. This read has soft-clipping at both the beginning and the end"
        print "I found a soft-clipping at the end! \nCurrent CIGAR: " + char + ". And current number of the CIGAR: " + str(num)
        # This read is soft-clipped at the end
        loc = read_len - num
        my_key = chrom + " " +  str(loc)
        if my_key not in d_ins.keys():
            d_ins[my_key] = [ seq[read_len-num+1:] ]
        else:
            d_ins[my_key].append( seq[read_len-num+1:] )
print u_ins
print d_ins
