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

RANGES = {
    'Sanger': (33, 73),
    'Solexa': (59, 104),
    'Illumina-1.3': (64, 104),
    'Illumina-1.5': (67, 104)
}
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
        
if __name__ == "__main__":
    main()
