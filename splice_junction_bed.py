#!/usr/bin/env python
# Note that some short exons will be omitted.
# Too short exons are not very common:
# fuy2@ummsres18:~/data/prepachytene/results/2015-01-14$ python get_splice_junction_bed.py ~/data/piPipes/common/mm9/UCSC.refSeq.Genes.bed12.gz 2>&1 > /dev/null | wc
# 31716  507456 2381304
# fuy2@ummsres18:~/data/prepachytene/results/2015-01-14$ python get_splice_junction_bed.py ~/data/piPipes/common/mm9/UCSC.refSeq.Genes.bed12.gz 2>/dev/null | wc
# 227241 2499651 13115805

# Usage: python splice_junction_bed.py ~/data/piPipes/common/mm9/UCSC.refSeq.Genes.bed12.gz 2>/dev/null  > splice_junctions.bed12
# Author: Yu Fu

import gzip, sys

# Use the left 35nt and right 35nt around the ss
s_range = 35

fh = gzip.open(sys.argv[1])
for line in fh:
    line = line.strip()
    col = line.split("\t")
    chrom, start, end, name, signal, strand, thick_start, thick_end, _, _, block_sizes_s, block_starts_s = col
    start = int(start)
    end = int(end)
    block_sizes_s = block_sizes_s.split(",")
    block_sizes = []
    for b in block_sizes_s:
        if b == "":             # Skip the last one, which is empty for every line
            continue
        block_sizes.append(int(b))
            
    block_starts_s = block_starts_s.split(",")
    block_starts = []
    for b in block_starts_s:
        if b == "":
            continue
        block_starts.append(int(b))
    
    for i in range(len(block_sizes)-1):
        if block_sizes[i] < 2 * s_range:
            print >>sys.stderr, "Found a exon that is shorter than 2 x " + str(s_range) + ". I will skip this one: " + name
            continue
        p = block_sizes[i] + block_starts[i]
        # These positions are relative to original chromStart
        first_half_end = p
        first_half_start = first_half_end - s_range
        second_half_start = block_starts[i+1]
        second_half_end = second_half_start + s_range
        # Notice that now chromStart has changed
        # Block count is always 2 since we are looking at the ss
        # Intron size + 2 * s_range
        second_half_start_new = (block_starts[i+1] - p) + s_range
        print "\t".join([chrom, str(start + first_half_start), str(start + second_half_end), name + "_sj_" + str(i), "0", strand, str(start + first_half_start), str(start + second_half_end), "0", "2", str(s_range) + "," + str(s_range) + "," , "0," + str(second_half_start_new) + ","])
        
        

