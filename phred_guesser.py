#!/usr/bin/env python

import sys, optparse


def calc_qual_range(l):
    """Calculate the range of quality scores for just one sequence
    """
    l_min = 255
    l_max = 0
    for c in l:
        if l_min > ord(c):
            l_min = ord(c)
        if l_max < ord(c):
            l_max = ord(c)
    return (l_min, l_max)


def main():
    counter = 0
    g_min = 255
    g_max = 0
    for line in sys.stdin:
        counter += 1
        if counter == 4:
            # Calculate the local min and max quality score
            l_min, l_max = calc_qual_range(line.strip())
            if g_min > l_min:
                g_min = l_min
            if g_max < l_max:
                g_max = l_max
            counter = 0
    print "The range of the quality scores is: [%d, %d]" % (g_min, g_max)
    print report(g_min, g_max)

def report(g_min, g_max):
    """Report phred33 or phred64
    """
    ranges = {
        'phred33': (33, 84), # Sanger
        # 'Solexa': (59, 104),
        'phred64': (64, 115), # Illumina-1.3
        # 'Illumina-1.5': (67, 104)
    }
    phred33_flag = False
    phred64_flag = False
    for k in ranges.keys():
        if g_min >= ranges[k][0] and g_max <= ranges[k][1]:
            phred33_flag = True
        if g_min >= ranges[k][0] and g_max <= ranges[k][1]:
            phred64_flag = True
    if phred33_flag and phred64_flag:
        print "Ambiguous"
    elif phred33_flag and not phred64_flag:
        print "phred33"
    elif phred64_flag and not phred33_flag:
        print "phred64"
    else:
        print "Something is wrong?"
    
if __name__ == "__main__":
    main()
