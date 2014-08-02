#!/usr/bin/env python

# Read in a fasta file and split it into multiple files
# Author: Yu Fu (yfu at yfu dot me)

import argparse

parser = argparse.ArgumentParser(description='Split a fasta file into multiple files')

parser.add_argument('-p', '--parts', help='The number of (almost) equal parts', required=True)
parser.add_argument('-f', '--file', help='The fasta file', required=True)
parser.add_argument('-v', '--verbose', help='Print lots of useless information', action="store_true")

args = parser.parse_args()
n = int(args.parts)
fn = args.file
verbose = args.verbose

total = 0
fh = open(fn, 'r')
for line in fh.readlines():
    line = line.strip()
    if line[0] == '>':
        total += 1

fh.close()

# Do notice that total might not be a multiple of parts, say 151 total line and 3 parts.

each = int(total / float(n))

fh = open(fn, 'r')
output = []

# Notice that inside the program, the index of files starts from 0
# and the filenames start from 1
for i in range(n):
    output.append( open(fn + '.' + str(i+1), "w") )

counter = -1;
for line in fh.readlines():
    if(line[0] == '>'):
        counter += 1
    line = line.strip()
    file_number = int(counter / each)
    # In order to put the last bit of the file into the last file...
    if( counter / each > n-1 ):
        file_number = n-1
    # print file_number, line
    if(verbose==True):
        print str(file_number) +"\t" + line
    print >>output[file_number], line

for i in range(n):
    output[i].close()
