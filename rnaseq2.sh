#!/usr/bin/env bash

export genome_fa=/data/fuy2/shared/dm3.fa
# transcriptome_gtf=/data/fuy2/shared/iGenome.gtf
# TODO: Use an array to iterate through many gff files
# transcriptome_gtf="/data/fuy2/shared/dmel_r5.55/gff/dmel-all-r5.55.gff.gz"
transcriptome_gtf="/data/fuy2/shared/flybase_r5.45/dmel-all-r5.45.gff.FlyBase.gene.gtf"

library_type="fr-firststrand"

label_a=$1
label_b=$2

sample_a_bams=$3 # Separated by comma
sample_b_bams=$4

out_dir=$(pwd)

step=1
cpu=$5

echo "Use $cpu CPUs."

[ ! -f .${step}.cuffdiff ] && \
cuffdiff -p $cpu -o $out_dir -L $label_a,$label_b -u --compatible-hits-norm -b $genome_fa --library-type $library_type \
     --library-norm-method geometric \
     -v --no-update-check \
     $transcriptome_gtf \
     ${sample_a_bams} \
     ${sample_b_bams} \
     2> $out_dir/cuffdiff.log && \
touch .${step}.cuffdiff

step=$((STEP+1))


