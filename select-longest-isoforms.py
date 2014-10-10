#!/usr/bin python

# Select the longest isoform for each gene
# Accepts a bed file as input
# Example: gtfToGenePred ~/data/piPipes/common/dm3/dm3.genes.gtf stdout -genePredExt | awk '{ $1=$12; print }' | genePredToBed stdin stdout | python run3-select-longest-isoforms.py | sort -k1,1 -k2,2n > longest_isoform_per_gene.bed

import fileinput

strand_info = {}
start_coor = {}
end_coor = {}
chrom_info = {}

for line in fileinput.input():
    line = line.strip()
    col = line.split()
    chrom, start, end, name, signal, strand = col[0:6]
    start = int(start)
    end = int(end)

    chrom_info[name] = chrom
    strand_info[name] = strand
    if name not in start_coor:
        start_coor[name] = start
    elif start < start_coor[name]:
        start_coor[name] = start

    if name not in end_coor:
        end_coor[name] = end
    elif end > end_coor[name]:
        end_coor[name] = end

for i in chrom_info.keys():
    print chrom_info[i] + "\t" + str(start_coor[i]) + "\t" +  str(end_coor[i]) + "\t" + i + "\t" + "0" + "\t" +  strand_info[i]
