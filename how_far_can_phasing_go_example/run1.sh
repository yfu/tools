#!/usr/bin/env bash

# cp ../updown2k/GSV6_*.best_hit.bed ./

for i in ../GSV6_upstream_in_42A18.fa  ../GSV6_upstream_in_zip.fa; do
    prefix=$(basename ${i%.fa})
    ext=$(cat ${prefix}.dm3.best_hit.bed | awk '{ print 2000 - ($3 - $2) }')
    bedtools slop -l ${ext} -i ${prefix}.dm3.best_hit.bed -s -r 0 -g ~/data/shared/dm3/dm3.ChromInfo.txt > ${prefix}.dm3.2k.bed
done

for i in ../GSV6_downstream_in_42A18.fa  ../GSV6_downstream_in_zip.fa; do
    prefix=$(basename ${i%.fa})
    ext=$(cat ${prefix}.dm3.best_hit.bed | awk '{ print 2000 - ($3 - $2) }')
    bedtools slop -r ${ext} -i ${prefix}.dm3.best_hit.bed -s -l 0 -g ~/data/shared/dm3/dm3.ChromInfo.txt > ${prefix}.dm3.2k.bed
done

