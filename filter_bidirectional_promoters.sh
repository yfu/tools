gtf=$1
cat ${gtf} | grep start_codon | awk 'BEGIN{OFS="\t"} {if($7=="+")print }' > start_codons_forward_strand.gtf
cat ${gtf} | grep start_codon | awk 'BEGIN{OFS="\t"} {if($7=="-")print }' > start_codons_reverse_strand.gtf

bedtools closest -a start_codons_forward_strand.gtf -b start_codons_reverse_strand.gtf -S -D a > tmp.gtf
cat tmp.gtf | awk 'BEGIN{FS=OFS="\t"} {d=$19; if(d<0&&d>-10000){print} }' | cut -f 9,18 | sed -E 's/gene_id \"([^"]+)\".+gene_id \"([^"]+)\".+/\1\t\2/' | awk '{ print $1; print $2 }' | sort | uniq > genes_with_bidirectional_promoters.txt
gawk '{ if(FNR==NR) {d[$1]=1} else{ match($0, /gene_id "([^"]+)".+/, m); g=m[1]; if(g in d) {} else {print} } } ' genes_with_bidirectional_promoters.txt mm9.genes.gtf > mm9.no_genes_with_bidirectional_promoters.gtf

gtfToGenePred mm9.no_genes_with_bidirectional_promoters.gtf stdout | genePredToBed stdin stdout > mm9.no_genes_with_bidirectional_promoters.gene.bed
