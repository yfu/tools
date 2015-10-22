#!/usr/bin/env bash

chrominfo=$1
watson=$2
crick=$3
output=$4

if [[ "$#" -ne 4 ]]; then
    echo "Usage: merge_bigwig_watson_crick.sh chrominfo.txt Watson.bigWig Crick.bigWig output.bigWig"
    exit 1
else
    bigWigToBedGraph ${crick} stdout | awk '{ $4=-$4; print }' > ${crick}.reverse.bedGraph.tmp
    bedGraphToBigWig ${crick}.reverse.bedGraph.tmp ~/data/shared/hg19/hg19.ChromInfo.txt ${crick}.reverse.bigWig.tmp
    bigWigMerge ${watson} ${crick}.reverse.bigWig.tmp ${output}.bedGraph.tmp
    bedGraphToBigWig ${output}.bedGraph.tmp ${chrominfo} ${output}
    rm ${crick}.reverse.bedGraph.tmp {crick}.reverse.bigWig.tmp ${output}.bedGraph.tmp
fi
