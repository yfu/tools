#!/usr/bin/env python

# Usage: python mark_duplicates_umi.py -f test.sorted.bam >test.dup_marked.bam 2>test.dup_marked.log

DEBUG = False
# A FASTQ file is first processed by reformat_umi_fastq.py to put the
# barcode info in the fasta header. Then it is aligned by aligner such as
# STAR. The bam file can then be processed by this script.
# It will go through the bam file and search for r1 reads that align to the same
# position and 
import pysam as ps
from Bio.Seq import Seq
import sys
import argparse
import sys
from multiprocessing import Process, Queue, Pool

# Two cases:
# ------------------------->
#                                     <-------------------------
#             r1                                   r2
#
# |
# r1fwd
#######################################################################
# ------------------------->
#                                     <-------------------------
#             r2                                   r1
# |
# r2fwd
######################################################################
# For a sorted bam file, in case 1, the 5' end of r1 and template length can determine the template;
# in case 2, the 5' end of r2 and template length can determine the template
# 

parser = argparse.ArgumentParser(description='A pair of FASTQ files are first reformatted using reformat_umi_fastq.py and then is aligned to get the bam file. This script can parse the umi barcode in the name of each read to mark duplicates. You need to first sort the bam file by name and add the FM tag using -a option. Then run it without -a option to mark the duplicates. You need to tell which chromosome you want to process so that you can parallelize this script.')
parser.add_argument('-f', '--file', help='the input bam file', required=True)
parser.add_argument('-p', '--processes', help='number of processes', required=False, type=int, default=8)
# parser.add_argument('-c', '--chromosome', help='the chromosome that you want to process', required=True)
# parser.add_argument('-a', '--add-tag', help='add a FM (five prime end of the mate) tag as the preprocessing step', action="store_true")
parser.add_argument('-d', '--debug', help='turn on debug mode', action="store_true")
# parser.add_argument('-l', '--read-length', help='if read length is given, it can be used to more accurately mark duplicates', type=int, default=-1)
# parser.add_argument('-o', '--output', help='the output file', required=True)
args = parser.parse_args()
processes = args.processes
# True: output the sequences on the reference strand
# False: output the sequences in the orignal direction
global infile
infile = args.file
# bam = ps.AlignmentFile(infile, "rb")
# readlength = args.read_length
DEBUG = args.debug

def mark_duplicates(infile, chromosome):
    bam = ps.AlignmentFile(infile, "rb")
    out = ps.AlignmentFile(infile + "." + chromosome + ".bam", "wb", template=bam)
    print >>sys.stderr, "Processing chromosome: " + chromosome
    r1fwd = -1
    r2fwd = -1
    # 
    # dup_read_names = []
    # # Position + barcode as the key; count as the value
    # posbc = {}
    c = 0
    # Unique read IDs: read5 + barcode + template length
    r1fwd_ids = []
    r2fwd_ids = []
    # Stores the reads to be marked as duplicates (only the 2nd and later ones are marked; the first one is not)
    r1fwd_dup = []
    r2fwd_dup = []
    for read in bam.fetch(reference=chromosome):
        c += 1
        if c % 100000 == 0:
            print >>sys.stderr, "Processed " + str(c) + " entries..."
        read_n = read.query_name
        read_bc = read_n.split("_")[0]
        read_chr = read.reference_id
        read5 = read.reference_start
        if DEBUG:
            print >>sys.stderr, "read_info: 5'end: " + str(read5)
            print >>sys.stderr, "read_info: Reference ID: " + str(read.reference_id)
            print >>sys.stderr, "read_info: Barcode: " + read_bc
            print >>sys.stderr, "read_info: Read: " + str(read)        
        if not read.is_reverse and read.is_read1:
            read_id = str(read5) + read_bc + str(read.template_length)
            if DEBUG:
                print >>sys.stderr, read_id + "\t" + str(read)
            if read_id in r1fwd_ids:
                if DEBUG:
                    print >>sys.stderr, "Found a duplicate (r1fwd): " + str(read)
                r1fwd_dup.append(read_n)
            else:
                r1fwd_ids.append(read_id)
        elif not read.is_reverse and read.is_read2:
            # for a read that maps to the reverse strand, the read.template_length is still positive (not like those in the
            # bam file)
            # read_id = str(read5) + read_bc + str(-read.template_length)
            read_id = str(read5) + read_bc + str(read.template_length)
            if DEBUG:
                print >>sys.stderr, read_id + "\t" + str(read)
            if read_id in r2fwd_ids:
                if DEBUG:
                    print >>sys.stderr, read_n, "Found a duplicate (r2fwd)" + str(read)
                r2fwd_dup.append(read_n)
            else:
                r2fwd_ids.append(read_id)
    bam.close()
    # print >>sys.stderr, r1fwd_ids
    # print >>sys.stderr, r1fwd_dup
    # print >>sys.stderr, r2fwd_ids    
    
    # The 2nd pass:
    bam = ps.AlignmentFile(infile, "rb")
    for read in bam.fetch(reference=chromosome):
        read_n = read.query_name
        if read_n in r1fwd_dup or read_n in r2fwd_dup:
            # Add 0x400 flags for PCR duplicates
            read.flag = read.flag | 0x400
        out.write(read)
    bam.close()
    out.close()

def mark_duplicates_worker(chromosome):
    mark_duplicates(infile, chromosome)
    
if __name__ == "__main__":
    # mark_duplicates(bam, processes=processes)
    # mark_duplicates_helper(infile, processes=processes)
    bam_tmp = ps.AlignmentFile(infile, "rb")
    refs = bam_tmp.references
    f = lambda x: mark_duplicates(infile, x)
    p = Pool(processes)
    # Goddamned Pool does not accept lambda functions!
    p.map(mark_duplicates_worker, refs)
    
