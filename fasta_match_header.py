#!/usr/bin/env python

# Print out fasta entries that match the pattern.
# seqtk subseq is only able to find exact matches.

import re, argparse, sys

parser = argparse.ArgumentParser(description='Print out fasta entries that match the pattern')

parser.add_argument('-p', '--pattern', help='The pattern', required=True)

args = parser.parse_args()
pat = args.pattern

flag = False
for line in sys.stdin:
    line = line.strip()
    if line[0] == '>':          # A fasta header
        name = line[1:]
        if re.search(pat, name) is None:
            flag = False
            next
        else:
            flag = True
            print line
    else:
        if flag == True:
            print line
