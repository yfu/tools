#!/usr/bin/env python

# Note that some short exons will be omitted.
# Too short exons are not very common:
# fuy2@ummsres18:~/data/prepachytene/results/2015-01-14$ python get_splice_junction_bed.py ~/data/piPipes/common/mm9/UCSC.refSeq.Genes.bed12.gz 2>&1 > /dev/null | wc
# 31716  507456 2381304
# fuy2@ummsres18:~/data/prepachytene/results/2015-01-14$ python get_splice_junction_bed.py ~/data/piPipes/common/mm9/UCSC.refSeq.Genes.bed12.gz 2>/dev/null | wc
# 227241 2499651 13115805

# Usage: python exon_intron_junction_bed.py ~/data/piPipes/common/mm9/UCSC.refSeq.Genes.bed12.gz > exon_intron_boudaries.bed6 
# To get fasta: cat splice_junctions.50nt.bed12 | awk 'BEGIN{OFS=FS="\t"} { $4=$4 "_" $1 "_" $2 "_" $3 "_" $6; print }' | bedtools getfasta -split -bed - -fi ~/data/piPipes/common/mm9/mm9.fa -fo - -fullHeader -s -name  | less

# Author: Yu Fu

import gzip, sys

# Use the left 35nt and right 35nt around the ss
s_range = 50

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
    
    for i in range(len(block_sizes)):
        # In this case, The exon is too short and the two reported ranges will overlap. To reduce multimappers, ignore such exons
        if block_sizes[i] < 2 * s_range:
            continue
        if i == 0:
            print "\t" . join([chrom, str(start + block_sizes[i] - s_range), str(start + block_sizes[i] + s_range), name + "_exon_intron_boundary_" + str(i) + "_r", "0", strand])
            # print "!"
        elif i == len(block_sizes)-1:
            print "\t" . join([chrom, str(start + block_starts[i] - s_range), str(start + block_starts[i] + s_range), name + "_exon_intron_boundary_" + str(i) + "_l", "0", strand])
            # print "?"
        else:
            print "\t" . join([chrom, str(start + block_starts[i] - s_range), str(start + block_starts[i] + s_range), name + "_exon_intron_boundary_" + str(i) + "_l", "0", strand])
            print "\t" . join([chrom, str(start + block_starts[i] + block_sizes[i] - s_range), str(start + block_starts[i] + block_sizes[i] + s_range), name + "_exon_intron_boundary_" + str(i) + "_r", "0", strand])

            # print "?!"
