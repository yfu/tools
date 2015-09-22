#!/usr/bin/env bash

# Given a SAM file, output the # reads that overlap with transposon annotations
sam=$1

prefix=${sam%.sam}
prefix=${prefix%.bam}
cpu=8
genome=mm10

samtools view -b -F 0x100 ${sam} | bedtools bamtobed -bed12 -i - | bedtools bed12tobed6 -i /dev/stdin | awk -v OFS="\t" '{ l=length($4); r=substr($4, l, 1); if(r==1) { $6=$6=="+"?"-":"+"; print } else { print } }' | sort -k1,1 -k2,2n --parallel=${cpu} > ${prefix}.bed6

cmd=${prefix}.transposon_intersection.cmd
echo "" > ${cmd}
if [[ ${genome} == "mm10" ]]; then
    for i in DNA LINE LTR rRNA Satellite Simple_repeat SINE tRNA; do
	echo "bedtools intersect -sorted -s -a ${prefix}.bed6 -b ~/data/shared/mm10/UCSC.rmsk.${i}.sorted.bed -wa -wb -f 0.9 | sort -k1,1 -k2,2n | bedtools groupby -i - -g 7,8,9,10 -c 4 -o count_distinct " >> ${cmd}
    done
fi

cat ${cmd} | parallel -j ${cpu}
