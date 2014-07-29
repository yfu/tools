#!/usr/bin/env python

# This is just an ad-hoc script for calculating histograms for Ping-Pong signal.
# Consider generalizing it before you use it for other purposes
# This dictionary records NTMs when the reads are mapped to the reverse index
# Usage: samtools -F 0x10 xxx.bam | python pp_hist.py my.bed2 > hist.txt
# Author: Yu Fu (yfu at yfu dot me)
import sys

def main():
    # Read the file to get the ntm
    ntm = get_ntm(sys.argv[1])
    # Read the file again to get the histogram
    fin = open(sys.argv[1])
    pp_hist = {}
    for line in fin.readlines():
        line = line.strip()
        e = line.split()
        if e[1] != '16':
            # Just in case someone forgets to filter the sam file
            sys.stderr.write("Found " + e[1] + " at column 2.")
            sys.exit("The second column is expected to be 16. Please only send me sam files that have 16 (samtools view -f 0x10) as their second column. Bye!")
        # Example: scaffold28950_1156_1177_1_1_-
        [read_chrom, read_start, read_end, read_copy, read_ntm, read_strand] = e[0].split("_")
        read_start, read_end, read_copy, read_ntm = int(read_start), int(read_end), int(read_copy), int(read_ntm)
        # Example: orig_scaffold28950_1155_1175_+_ext_scaffold28950_1125_1178_+_1_3
        [_, orig_chrom, orig_start, orig_end, orig_strand, _, ext_chrom, ext_start, ext_end, ext_strand, orig_copy, orig_ntm] = e[2].split("_")
        orig_start, orig_end, ext_start, ext_end, orig_copy, orig_ntm = int(orig_start), int(orig_end), int(ext_start), int(ext_end), int(orig_copy), int(orig_ntm)
        if from_the_same_locus(read_chrom, read_start, read_end, read_strand, ext_chrom, ext_start, ext_end, ext_strand) == False:
            # print line
            # These two variables exist solely for convenience
            ext_copy = orig_copy
            ext_ntm = orig_ntm
            seq = e[9]
            # Use 0-based coordinates, as those in the bed format
            pos = int(e[3]) - 1 
            weight = (float(read_copy) / read_ntm) * (float(orig_copy) / orig_ntm) / ntm[e[9]]
            # print "\t" .join( [str(read_copy), str(read_ntm), str(orig_copy), str(orig_ntm)] )
            # print weight
            offset = calc_offset(orig_start, orig_end, orig_strand, ext_start, ext_end, ext_strand, read_start, read_end, read_strand, pos)
            # print offset
            if offset not in pp_hist.keys():
                pp_hist[offset] = weight
            else:
                pp_hist[offset] += weight
        else:
            continue
    for k in pp_hist.keys():
        print str(k) + "\t" + str(pp_hist[k])
    
def calc_offset(orig_start, orig_end, orig_strand, ext_start, ext_end, ext_strand, read_start, read_end, read_strand, pos):
    ''' Given the information for the reads and the reverse index, this function
    calculates the offset needed for the Ping-Pong signal.
    '''
    read_len = read_end - read_start
    if(ext_strand == '+'):
        # return 0
        return pos + read_len - (orig_start - ext_start)
    else:
        # Calculation for the '-' case is a bit tricky...
        # return 0
        # return (ext_end - orig_start) - (pos - 1)
        # return pos + 1 - (ext_end - orig_end)
        return pos + read_len - (ext_end - orig_end)
        
def from_the_same_locus(read_chrom, read_start, read_end, read_strand, ext_chrom, ext_start, ext_end, ext_strand):
    ''' In the upstream analysis, it is often found that the reads map to another reads from the same locus. For example, the first read comes 
    from chromosome 1: 1-20, and the second read comes from chromosome 1:2-21. When I extend the second read to make the reference, obvisouly, 
    I can still map read 1 to the reference made from read 2. This function exists to exclude these pairs by determining if they are from overlapping
    genomic regions.'''

    if(read_chrom == ext_chrom):        
        # As long as one end of the read falls into the range of the reference, report True.
        if( (read_start > ext_start and read_start < ext_end) or (read_end > ext_start and read_end < ext_end) ):
            return True
    # print "\t" . join([str(read_chrom), str(read_start), str(read_end), str(read_strand), str(ext_chrom), str(ext_start), str(ext_end), str(ext_strand)])
    return False

def get_ntm(fn):
    ''' Read a sam file and report NTMs
    '''
    ntm = {}
    fin = open(fn)
    counter = 0
    for line in fin.readlines():
        if (counter % 10000 == 0):
            print "Processed " + str(counter) + " alignments..."
        counter += 1
        line = line.strip()
        e = line.split()
        if e[1] != '16':
            # Just in case someone forgets to filter the sam file
            sys.stderr.write("Found " + e[1] + " at column 2.")
            sys.exit("The second column is expected to be 16. Please only send me sam files that have 16 (samtools view -f 0x10) as their second column. Bye!")
        # Example: scaffold28950_1156_1177_1_1_-
        [read_chrom, read_start, read_end, read_copy, read_ntm, read_strand] = e[0].split("_")
        read_start, read_end, read_copy, read_ntm = int(read_start), int(read_end), int(read_copy), int(read_ntm)
        # Example: orig_scaffold28950_1155_1175_+_ext_scaffold28950_1125_1178_+_1_3
        [_, orig_chrom, orig_start, orig_end, orig_strand, _, ext_chrom, ext_start, ext_end, ext_strand, orig_copy, orig_ntm] = e[2].split("_")
        orig_start, orig_end, ext_start, ext_end, orig_copy, orig_ntm = int(orig_start), int(orig_end), int(ext_start), int(ext_end), int(orig_copy), int(orig_ntm)
        if from_the_same_locus(read_chrom, read_start, read_end, read_strand, ext_chrom, ext_start, ext_end, ext_strand) == False:
            seq = e[9]
            if(seq not in ntm.keys()):
                ntm[seq] = 1
            else:
                ntm[seq] += 1
        else:
            continue
    return ntm

if __name__ == "__main__":
    main()
