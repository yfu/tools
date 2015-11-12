#!/usr/bin/env python

import gzip
import sys
import argparse

def get_header_with_umi(header, umi):
    """This function inserts a UMI after the '@' symbol, making the downstream analysis easier
"""
    return "@" + umi + "_" + header[1:]
def is_good_phred(phred, qc):
    """Check if all the qualities are above qc, that is, if one 
nt has a quality score below qc, this function returns False
"""
    ret = True
    for i in phred:
        if i < qc:
            ret = False
            break
    return ret

# print phred_checker("BBC", "B")
# print phred_checker("ABC", "B")

parser = argparse.ArgumentParser(description='A script to reformat r1 reads in a UMI fastq file so that the name of each record contains the UMI')
parser.add_argument('-l', '--left', help='the input fastq file for r1. If you want to pipe in the input, you can use /dev/stdin', required=True)
parser.add_argument('-r', '--right', help='the input fastq file for r2. If you want to pipe in the input, you can use /dev/stdin', required=True)
parser.add_argument('-L', '--left-out', help='the output fastq file for r1', required=True)
parser.add_argument('-R', '--right-out', help='the output fastq file for r2', required=True)
parser.add_argument('-q', '--quality', help='Quality (phred quality score) cutoff for UMI. Default is 20, that is UMI with qualities >= 20 will be kept. This program assumes the phred quality scores in the fastq file are using sanger format', required=False, type=int, default=20)
parser.add_argument('-D', '--debug', help='Turn on debugging mode', action="store_true")

# Sanger format can encode a Phred quality score from 0 to 93 using ASCII 33 to 126 (although in raw read data the Phred quality score rarely exceeds 60, higher scores are possible in assemblies or read maps). Also used in SAM format.[4] Coming to the end of February 2011, Illumina's newest version (1.8) of their pipeline CASAVA will directly produce fastq in Sanger format, according to the announcement on seqanswers.com forum.[5]

args = parser.parse_args()
DEBUG=args.debug
if DEBUG:
    print >>sys.stderr, "Debugging mode is on"
# Quality cutoff
qc = chr(args.quality + 33)
print >>sys.stderr, "Quality cutoff in ASCII:\t" + qc
fn1 = args.left
fn2 = args.right

out1 = open(args.left_out, "w")
out2 = open(args.right_out, "w")

c = 0

umi_len = 5
umi_locator = 'GGG'
umi_locator_len = len(umi_locator)

# Trim one nucleotide after the GGG
umi_downstream = 'T'
umi_downstream_len = len(umi_downstream)
# Those reads without or with GGG
n_without_locator = 0
n_with_locator = 0
# Those reads with N's before GGG
n_with_ambiguous_umi = 0
# Those with A/C/G after GGG
n_with_wrong_padding = 0
# Those having bad quality UMI w/ GGG and T and w/o N's
n_bad_quality_umi = 0
# The rest
n_good_reads = 0
f1 = open(fn1)
f2 = open(fn2)

while True:
    c += 1
    if c%4==1:
        r1_name = f1.readline().strip()
        r2_name = f2.readline().strip()
        if not r1_name:
            break
    elif c%4==2:
        r1_seq = f1.readline().strip()
        r2_seq = f2.readline().strip()
    elif c%4==3:
        r1_info = f1.readline().strip()
        r2_info = f2.readline().strip()
    else:
        r1_qual = f1.readline().strip()
        r2_qual = f2.readline().strip()
        if DEBUG == True:
            print '-' * 80            
            print "Original qual:\t" + r1_qual            
            print "Original read:\t" + r1_seq
            print "Locator:\t" + ' ' * 5 + r1_seq[umi_len : umi_len + umi_locator_len]
            print "UMI:\t\t" + r1_seq[0: umi_len]
            print "What's left:\t" + ' ' * 9 + r1_seq[umi_len + umi_locator_len + umi_downstream_len: ]
        if r1_seq[umi_len : umi_len + umi_locator_len] == umi_locator:
            n_with_locator += 1
            my_umi = r1_seq[0:umi_len]
            my_umi_qual = r1_qual[0:umi_len]
            if my_umi.find("N") == -1:
                if r1_seq[umi_len + umi_locator_len] == umi_downstream:
                    if is_good_phred(my_umi_qual, qc):
                        my_seq = r1_seq[umi_len + umi_locator_len + umi_downstream_len: ]
                        n_good_reads += 1
                        print >>out1, get_header_with_umi(r1_name, my_umi)
                        print >>out1, my_seq
                        print >>out1, r1_info
                        print >>out1, r1_qual[umi_len + umi_locator_len + umi_downstream_len: ]
                        print >>out2, get_header_with_umi(r2_name, my_umi)
                        print >>out2, r2_seq
                        print >>out2, r2_info
                        print >>out2, r2_qual
                    else:
                        n_bad_quality_umi += 1
                        if DEBUG:
                            print "*Found one with low quality UMI:\t" + r1_seq
                else:
                    n_with_wrong_padding += 1
                    if DEBUG:
                        print "*Found one with wrong padding nucleotide (A/C/G):\t" + r1_seq
            else:
                if DEBUG == True:
                    print "*Found one with N's in UMI:\t" + r1_seq
                n_with_ambiguous_umi += 1
        else:
            if DEBUG == True:
                print "*Found one without locator:\t" + r1_seq
            n_without_locator += 1
f1.close()
f2.close()
out1.close()
out2.close()
print >>sys.stderr, "-" * 80
print >>sys.stderr, "Total:\t" + str((c-1)/4)
print >>sys.stderr, "Reads w/o locator:\t" + str(n_without_locator)
print >>sys.stderr, "Reads w/ locator:\t" + str(n_with_locator)
print >>sys.stderr, "-" * 80
print >>sys.stderr, "Reads w/ N's in UMI:\t" + str(n_with_ambiguous_umi)
print >>sys.stderr, "Reads w/ wrong padding nt (A/C/G):\t" + str(n_with_wrong_padding)
print >>sys.stderr, "Reads w/ low-quality UMI:\t" + str(n_bad_quality_umi)
print >>sys.stderr, "-" * 80
print >>sys.stderr, "Reads w/ proper UMI:\t" + str(n_good_reads)
