#!/usr/bin/env bash

# Given a fasta file from Trinity, an evalue and a sequence, this script will output the sequences in the Trinity assembly that are similar to the input sequence
# You probably need to trim the sequence headers in the input fasta file
msg="Usage: trinity_find_similar_seq.sh Trinity.fasta 1e-30 myseq.fasta"

if [ $# -lt 3 ]; then
    echo $msg
    exit 1
fi
transcriptome=$1
evalue=$2
seq=$3
# Change this variable if you want to use other types of BLAST
blast_cmd=tblastx
n_threads=8

prefix=${transcriptome%.fasta}
prefix=${prefix%.fa}

output_prefix=${seq%.fasta}
output_prefix=${seq%.fa}.${blast_cmd}

if [ ! -f Trinity.short_header.fasta.nhr ]; then
    makeblastdb -dbtype nucl -in ${transcriptome}
fi

output_blast=${output_prefix}.txt
output_blast_table=${output_prefix}.table.txt

tblastx -db ${transcriptome} -query ${seq} -num_threads ${n_threads} > ${output_blast}
tblastx -db ${transcriptome} -query ${seq} -num_threads ${n_threads} -outfmt 7 > ${output_blast_table}

namelist=${output_prefix}.RNA1.e_${evalue}.namelist
grep -v "#" ${output_blast_table} | awk -v evalue=${evalue} '{ if($11 < evalue) {print} }' | cut -f2 | sed -E 's/_seq[0-9]+//' | sort | uniq > ${namelist}

while read gene_pat; do
    echo ${gene_pat}
    cat ${transcriptome} | fasta_match_header.py -p ${gene_pat} > ${output_prefix}.${gene_pat}.fasta
done < ${namelist}
