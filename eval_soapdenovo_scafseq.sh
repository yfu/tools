
scafSeq=$1

prefix=${scafSeq%.scafSeq}
prefix=${prefix%.fa}
tmp_file=${prefix}.assemblathon.fa

awk '{ if($0 ~ /^>/) { if($0~/>scaff/){print; flag=1} else {flag=0}} else { if(flag==1){print} } }' ${scafSeq} > ${tmp_file}

output_summary=${prefix}.assemblathon.summary
# output_csv=${prefix}.assemblathon.csv
# output_graph=${prefix}.assemblathon.graph
assemblathon_stats.pl -n 100 -graph -csv ${tmp_file} &> ${output_summary}
rm ${tmp_file}
