#!/usr/bin/env bash

# This script will extract miRNA sequences from mature.fa and hairpin.fa by miRBase.
## seqtk does not output spaces in the header line of fasta. Replace those spaces with comma and replace them back in the end
input=$1
tmp_fa=${input}.tmp
tmp_namelist=${tmp_fa}.namelist

tr ' ' ',' < ${input} > ${tmp_fa}
grep -E 'melanogaster|mori' ${tmp_fa} | sed -E 's/^>//' > ${tmp_namelist}

output=${input%.fa}
output=${output%.fasta}.dme_bmo.fa

seqtk subseq ${tmp_fa} ${tmp_namelist} | tr ',' ' ' > $output
# rm ${tmp_fa} ${tmp_namelist}
