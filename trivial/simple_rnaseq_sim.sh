#!/usr/bin/env bash

# Notice that in the name column of the bed file, the coordinates are for read2, which is forward. Read1 is reverse


bedtools makewindows -b 42AB.bed -w 1000 -s 1 -i srcwinnum | awk '{ if($3-$2==1000) print }' | awk 'BEGIN{OFS="\t"} {$4 = $4 "_" $1 ":" $2 "-" $3 "/2"; print $0, ".", "+"} ' > read2.bed
bedtools makewindows -b 42AB.bed -w 1000 -s 1 -i srcwinnum | awk '{ if($3-$2==1000) print }' | awk 'BEGIN{OFS="\t"} {$4 = $4 "_" $1 ":" $2 "-" $3 "/1"; $2 += 400; $3 += 400; print $0, ".", "-"} ' > read1.bed

# Add more depth
# This should be able to provide 100 coverage for each base (except for the bases at the very beginning and the very end)
bedtools getfasta -fi ~/data/piPipes/common/dm3/dm3.fa -bed read1.bed -fo - -name -s > read1.fa
bedtools getfasta -fi ~/data/piPipes/common/dm3/dm3.fa -bed read2.bed -fo - -name -s > read2.fa
