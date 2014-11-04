
scafSeq=$1

prefix=${scafSeq%.scafSeq}
prefix=${prefix%.fa}
tmp_file=${prefix}.assemblathon.fa
scaf_len_file=${prefix}.assemblathon.scafLengths

# The file contains only lengths of sigleton contigs
contig_only_len_file=${prefix}.assemblathon.contigLengths

echo "Outputing scaffold lengths into ${scaf_len_file}" 1>&2
awk '{ if($0 ~ /^>/) { if($0~/^>scaff/){print; flag=1} else {flag=0}} else { if(flag==1){print} } }' ${scafSeq} | tr ' ' '_' | fasta_lengths.pl > ${scaf_len_file}

echo "Outputing length of singleton contigs into ${contig_only_len_file}" 1>&2
awk '{ if($0 ~ /^>/) { if($0~/^>C/){print; flag=1} else {flag=0}} else { if(flag==1){print} } }' ${scafSeq} | tr ' ' '_' | fasta_lengths.pl > ${contig_only_len_file}

echo "assemblathon_stats.pl is evaluating ${scafSeq}" 1>&2
awk '{ if($0 ~ /^>/) { if($0~/>scaff/){print; flag=1} else {flag=0}} else { if(flag==1){print} } }' ${scafSeq} > ${tmp_file}

output_summary=${prefix}.assemblathon.summary
# output_csv=${prefix}.assemblathon.csv
# output_graph=${prefix}.assemblathon.graph
assemblathon_stats.pl -n 100 -graph -csv ${tmp_file} &> ${output_summary}

rm ${tmp_file}
