#!/usr/bin/env bash


# zcat ~/data/piPipes/common/dm3/dm3.piRNAcluster.bed.gz | grep 42AB > 42AB.bed6
# bedtools maskfasta -fi ~/data/piPipes/common/dm3/dm3.fa -bed 42AB.bed6 -fo dm3.42AB_masked.fa
# bedtools getfasta -fi ~/data/piPipes/common/dm3/dm3.fa -bed 42AB.bed6 -fo - -name > 42AB.fasta

win_size=23

bedtools makewindows -b 42AB.bed6 -w ${win_size} -s 1 | awk -v l=${win_size} '{ if($3-$2==l) { print } }' > 42AB.${win_size}nt_win.bed
bedtools getfasta -fi ~/data/piPipes/common/dm3/dm3.fa -bed 42AB.${win_size}nt_win.bed -fo - > 42AB.${win_size}nt_win.fa

bowtie -f -v 0 -a --best --strata -p 32 -S ~/data/piPipes/common/dm3_chengjian/BowtieIndex/genome 42AB.${win_size}nt_win.fa | samtools view -b -F 0x4 - > 42AB.${win_size}nt_win.bam

bedtools bamtobed -i 42AB.${win_size}nt_win.bam > 42AB.${win_size}nt_win.bed
# cat  42AB.23nt_win.sam | awk '{ if($3!="GSV6") { print $1 } }' | sed -E 's/([0-9a-zA-Z]+):([0-9]+)-([0-9]+)/\1\t\2\t\3/' | sort -k1,1 -k2,2n | bedtools merge -i - > GSV6.non_unique.bed
cat 42AB.${win_size}nt_win.bed | sed -E 's/([0-9a-zA-Z]+):([0-9]+)-([0-9]+)/\1\t\2\t\3/' | awk '{ if($1==$4 && $2==$5 && $3 == $6) {} else {print} }' | cut -f4,5,6,7,8 | awk -v OFS="\t" '{ print $1, $2, $3, ".", 0, $5 }' | sort -k1,1 -k2,2n --parallel=8 | bedtools merge -i - > 42AB.${win_size}nt_win.non_uniq.bed

cat 42AB.${win_size}nt_win.non_uniq.bed | awk -v OFS="\t" '{ print $1, $2, $3, 1}' > 42AB.${win_size}nt_win.non_uniq.bedGraph
bedGraphToBigWig 42AB.${win_size}nt_win.non_uniq.bedGraph ~/data/piPipes/common/dm3/dm3.ChromInfo.txt 42AB.${win_size}nt_win.non_uniq.bigWig

## 23 and 29nt do not make much difference
bedtools complement -i 42AB.23nt_win.non_uniq.bedGraph -g ~/data/piPipes/common/dm3/dm3.ChromInfo.txt | awk -v OFS="\t" '{ print $0, $3-$2 }' | sort -k4,4nr > 42AB.23nt_win.uniq.bedGraph
