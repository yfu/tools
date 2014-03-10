#!/usr/bin/env bash

export genome_fa=/data/fuy2/shared/dm3.fa
transcriptome_gtf=/data/fuy2/shared/iGenome.gtf
library_type="fr-firststrand"

label_a=$1
label_b=$2

sample_a_prefix=$3
sample_b_prefix=$4

out_dir=$(pwd)

step=1
cpu=20

[ ! -f .${step}.cuffdiff ]
cuffdiff -o $out_dir -L $label_a,$label_b -p $CPU -u --compatible-hits-norm -b $genome_fa --library-type $library_type \
     --library-norm-method geometric \
     --no-update-check \
     $transcriptome_gtf \
     $out_dir/${sample_a_prefix}*.sorted.bam \
     $out_dir/${sample_b_prefix}*.sorted.bam \
     2> $out_dir/cuffdiff.log && \
touch .${step}.cuffdiff

step=$((STEP+1))


