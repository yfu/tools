#!/usr/bin/env bash

# zcat ~/data/piPipes/common/dm3/Zamore.transposon.group*.gz > Zamore.transposon.all_groups.bed
# for i in 0 1 2 3; do
#     zcat ~/data/piPipes/common/dm3/Zamore.transposon.group${i}.bed.gz | cut -f4 | sed -E 's/(.+)\.[0-9]+/\1/' | sort -u | awk -v number=${i} '{ print $1 "\t" "g" number }' > group${i}.name
# done
# cat group{0,1,2,3}.name > groups.name


# input=~/data/gfp_pirna/data/RSRZ.CC4.br1.tr1.ovary/genome_mapping/RSRZ.CC4.br1.tr1.ovary.x_rRNA.dm3_chengjian.sorted.unique.bed12
# nf=0.01
input=$1
nf=$2
prefix=$(basename ${input})
prefix=${prefix%.bed12}
transposon_bed=Zamore.transposon.all_groups.bed
sense_bed=${prefix}.vs_transposon.S.bed
antisense_bed=${prefix}.vs_transposon.AS.bed
# bedtools intersect -f 0.99 -a ${input} -b ${transposon_bed} -s -split -wa -wb > ${sense_bed}
# bedtools intersect -f 0.99 -a ${input} -b ${transposon_bed} -S -split -wa -wb > ${antisense_bed}
sense_count=${prefix}.vs_transposon.S.count
antisense_count=${prefix}.vs_transposon.AS.count
count=${prefix}.vs_transposon.count

# Remember to divide the numbers by 2 because the unit is pairs of reads
cat ${sense_bed} | awk '{ print $16 "\t" $4 }' | sort -u | awk '{ s[$1]+=1 } END{ for(i in s) { print i "\t" s[i] } }'     | sed -E 's/(.+)\.[0-9]+\t([0-9]+)/\1\t\2/' | awk '{ s[$1]+=$2 } END{ for(i in s) { print i "\t" s[i]/2 } }' > ${sense_count}
cat ${antisense_bed} | awk '{ print $16 "\t" $4 }' | sort -u | awk '{ s[$1]+=1 } END{ for(i in s) { print i "\t" s[i] } }' | sed -E 's/(.+)\.[0-9]+\t([0-9]+)/\1\t\2/' | awk '{ s[$1]+=$2 } END{ for(i in s) { print i "\t" s[i]/2 } }' > ${antisense_count}

awk -v OFS="\t" -v nf=${nf} '{ if(NR==FNR) {c[$1]=$2} else { if($1 in c) {my=c[$1]} else {my=0}; print $1, "sense", $2, my, my * nf }}' ${sense_count} groups.name > ${count}.tmp
awk -v OFS="\t" -v nf=${nf} '{ if(NR==FNR) {c[$1]=$2} else { if($1 in c) {my=c[$1]} else {my=0}; print $1, "antisense", $2, my, my * nf }}' ${antisense_count} groups.name >> ${count}.tmp

sort -k3,3 -k1,1 -k2,2 ${count}.tmp > ${count}
