#!/usr/bin/env bash

# Given a bed2 file (excluding miRNA reads, include both unique and multiple mappers), output the following bed2 files (based on length and mappability):

# all mappers
# all mappers including only siRNAs
# all mappers including only piRNAs

# unique mappers
# unique mappers including only siRNAs
# unique mappers including only piRNAs

# Usage: bed2_length_separation.sh mm10 a/b/c/prefix.all.bed2 output_dir
assembly=$1
bed2=$2
output_dir=$3
source /data/fuy2/piPipes/common/${assembly}/variables

if [[ -z ${output_dir} ]]; then
    output_dir="."
else
    mkdir -p ${output_dir}
fi
# Just in case the bed2 file does not end up with .all.bed2
prefix=${bed2%.bed2}
prefix=${prefix%.all}
prefix=$(basename $prefix)

# all mappers
cp ${bed2} $output_dir/

echo ${output_dir}/${prefix}.all.siRNA.bed2
awk -v b=${siRNA_bot} -v t=${siRNA_top} '{ l=$3-$2; if( l>=b && l<=t ) { print } }' ${bed2} > ${output_dir}/${prefix}.all.siRNA.bed2
awk -v b=${siRNA_bot} -v t=${siRNA_top} '{ l=$3-$2; if( l>=b && l<=t && $5==1) { print } }' ${bed2} > ${output_dir}/${prefix}.unique.siRNA.bed2
awk -v b=${piRNA_bot} -v t=${piRNA_top} '{ l=$3-$2; if( l>=b && l<=t ) { print } }' ${bed2} > ${output_dir}/${prefix}.all.piRNA.bed2
awk -v b=${piRNA_bot} -v t=${piRNA_top} '{ l=$3-$2; if( l>=b && l<=t && $5==1) { print } }' ${bed2} > ${output_dir}/${prefix}.unique.piRNA.bed2




