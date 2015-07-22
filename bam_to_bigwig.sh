#!/usr/bin/env bash

# Usage: bash -x run6-bam-to-bigwig.sh ~/data/shared/hg19/hg19.ChromInfo.txt Mattick.RNaseR.HeLa.BP.rep1.fq.gz.hg19.Aligned.out.sam 0.0001
# input is sorted sam

if [[ "$#" -lt 2 ]]; then
    echo "Usage: bam_to_bigwig.sh your.ChromInfo.txt your.sam normalization_factor"
    exit 2
fi

chrominfo=$1
input=$2
nf=$3
cpu=16
prefix=${input%.sam}
bam_sorted=${prefix}.bam
samtools view -b ${input} | samtools sort -l 0 -@ ${cpu} -O bam -T ${RANDOM} - > ${bam_sorted}
samtools index ${bam_sorted}

bed_sorted=${prefix}.bed
bedtools bamtobed -bed12 -i ${bam_sorted} > ${bed_sorted}

para_cmd=${RANDOM}${RANDOM}.parallel.commands
echo "bedtools genomecov -split -bg -i ${bed_sorted} -strand + -g ${chrominfo} | awk -v nf=$nf 'BEGIN {OFS=\"\t\"} {\$4=\$4*nf; print}' > ${prefix}.Watson.bedGraph && bedGraphToBigWig ${prefix}.Watson.bedGraph ${chrominfo} ${prefix}.Watson.bigWig" > ${para_cmd}
echo "bedtools genomecov -split -bg -i ${bed_sorted} -strand - -g ${chrominfo} | awk -v nf=$nf 'BEGIN {OFS=\"\t\"} {\$4=-\$4 * nf; print}' > ${prefix}.Crick.bedGraph && bedGraphToBigWig ${prefix}.Crick.bedGraph ${chrominfo} ${prefix}.Crick.bigWig" >> ${para_cmd}
cat ${para_cmd} | parallel --progress
