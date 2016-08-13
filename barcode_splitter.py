#!/usr/bin/env python

#
# Split a big fastq file into several by barcode.
# TODO: parallelize this script, if possible.
# Author: Yu Fu
#

import sys, os, re
import optparse # For compatibility with Python 2.6 and lower
from optparse import OptionParser
import gzip
from struct import unpack


def is_gzipped(filename):
    # 1F 8B 08 00 / gz magic number
    magic = ('\x1f', '\x8b', '\x08')
    with open(filename, 'rb') as handle:
        s = unpack('ccc', handle.read(3))
        print s
    return s == magic

    
def calc_str_dist(barcode, s):
    """Calculate the distance between two strings of the same or different lengths:
    In either case, use the first few letters to determine the distance. For example,
    @FCD249BACXX:8:1101:1301:2095#ATCACGAT/1 belongs to barcode 'id1     ATCACG'
    """
    l = len(barcode)
    dist = 0
    for i in range(l):
        if barcode[i] != s[i]:
            dist += 1
    return dist


def main():
    usage = """A simple script to split a fastq file by barcode
Author: Yu Fu (yfu at yfu dot me)"""
    
    parser = OptionParser(usage=usage, version="%prog 0.1")
    parser.add_option("-f", "--fastq", action="store", type="string",
                      dest="fastq_fn", help="Input fastq file")
    
    parser.add_option("-b", "--barcode", action="store", type="string",
                      dest="barcode_fn", help="Barcode file")
    
    parser.add_option("-m", "--mismatches", action="store", type="int",
                      dest="mismatches", help="Number of mismatches allowed")
    
    parser.add_option("-p", "--prefix", action="store", type="string",
                      dest="fn_prefix", help="Optional filename prefix")
    parser.add_option("-v", "--verbose", action="store_true",
                      dest="verbose", help="Verbose mode")
    
    # parser.add_option("-n", "--name", action="store", type="string",
    #                   dest="name", help="The string that is added for each new fastq files")
    (options, args) = parser.parse_args()
    if options.fastq_fn == None:
        parser.print_help()
        sys.exit(1)

    fastq_fn = options.fastq_fn
    barcode_fn = options.barcode_fn
    mismatches = options.mismatches
    fn_prefix = options.fn_prefix
    if fn_prefix == None:
        fn_prefix = ""
    # print "*%s*" % (fn_prefix)
    verbose = options.verbose
    # name = options.name
    
    barcodes_names = {} # Seq as the key, name as the value
    barcodes_fh = {} # Seq as the key, file handle for output as value
    with open(barcode_fn) as barcode_in:
        for line in barcode_in.readlines():
            if not line:
                break;
            tmp = re.split(r"\s*", line)
            barcodes_names[tmp[1]] = tmp[0]
            barcodes_fh[tmp[1]] = ""
        
    p, s = os.path.split(fastq_fn)
    s = s.replace(".gz", "")
    s = s.replace(".fastq", "")
    s = s.replace(".fq", "")
    p = p + s
    print >>sys.stderr, "Output prefix: " + p
    barcodes_fh["unmatched"] = open(p + fn_prefix + "_unmatched.fq", "w")
    
    for k in barcodes_names.keys():
        # a.fastq would be split into a.id1.fastq, a.id2.fastq....
        barcodes_fh[k] = open(p  + fn_prefix + "_" + barcodes_names[k] + ".fq", "w")

    if is_gzipped(fastq_fn):
        print >>sys.stderr, "fastq file is gzipped"
        fastq_in = gzip.open(fastq_fn)
    else:
        print >>sys.stderr, "fastq file is not gzipped"
        fastq_in = open(fasdq_fn)

    while True:
        line1 = fastq_in.readline()
        rest_of_line = fastq_in.readline() + fastq_in.readline() \
          + fastq_in.readline()
        if not rest_of_line:
            break
        line1 = line1.strip()
        p1, p2 = re.split(r'[\s#]', line1)  # line1.split()
        seq_barcode = p2.split(":")[-1]

        if verbose == True:
            print >> sys.stderr, "Current sequence: " + line1
            print >> sys.stderr, "Current barcode: " + seq_barcode
        # print "The barcode from the current sequence is " + seq_barcode
        dist = 0
        matched = False
        for i in barcodes_names.keys():
            dist = calc_str_dist(i, seq_barcode)
            # print "Distance: " + str(dist) + "between " + i + "and " + seq_barcode
            # print "Trying to match: " + i
            if dist <= mismatches:
                if verbose == True:
                    print >> sys.stderr, "Current sequence matches barcode " + str(i)
                print >> barcodes_fh[i], line1
                print >> barcodes_fh[i], rest_of_line,
                # print "Successful match"
                matched = True
                # if the barcode from one seq already matches,
                # we don't need to try to match it with other barcodes
                break;
        if matched == False:
            if verbose == True:
                print >> sys.stderr, "Current sequence does not match any barcodes"
            print >> barcodes_fh["unmatched"], line1
            print >> barcodes_fh["unmatched"], rest_of_line,
        # print "-" * 72
    for k in barcodes_fh.keys():
        barcodes_fh[k].close()
                                
if __name__ == "__main__":
    main()

