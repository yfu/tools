#!/usr/bin/env bash


# Given a bed12 file for just one gene, this script outputs the gene plot (for all_mappers, uniq_mappers)
# 
# Usage: bed12_to_depth.sh one_gene.bed12 gene.name gene.length n_all_mappers n_uniq_mapper
# chrom_info=$1
bed12=$1
gene_n=$2
gene_l=$3

prefix=${bed12%.bed12}
prefix=${prefix%.bed6}
prefix=${prefix%.bed}

# Numbers of all mappers and of unique mappers for normalization purposes
n_all_mapper=$4
n_uniq_mapper=$5
nf_all=$( awk -v n_all_mapper=$n_all_mapper 'BEGIN{ print 1e6/n_all_mapper }' )
nf_uniq=$( awk -v n_uniq_mapper=${n_uniq_mapper} 'BEGIN{ print 1e6/n_uniq_mapper }' )
echo "Multiplier for all mappers: $nf_all"
echo "Multiplier for unique mappers: $nf_uniq"

prefix_all=${prefix}.all_mapper
bed6=${prefix_all}.bed6
prefix_uniq=${prefix}.uniq_mapper
bed6_uniq=${prefix_uniq}.bed6

bedtools bed12tobed6 -i ${bed12} | grep -v 'track ' > ${bed6}
awk '{ if($5=="1") {print} }' ${bed6} > ${bed6_uniq}

watson_bg_all=${prefix_all}.Watson.bedGraph
watson_depth_all=${prefix_all}.norm.Watson.depth
crick_bg_all=${prefix_all}.Crick.bedGraph
crick_depth_all=${prefix_all}.norm.Crick.depth

watson_bg_uniq=${prefix_uniq}.Watson.bedGraph
watson_depth_uniq=${prefix_uniq}.norm.Watson.depth
crick_bg_uniq=${prefix_uniq}.Crick.bedGraph
crick_depth_uniq=${prefix_uniq}.norm.Crick.depth

chrom_info=${gene_n}.${gene_l}.ChromInfo.txt
echo -e "${gene_n}\t${gene_l}" > ${chrom_info}

bed6_to_bedgraph.py ${chrom_info} ${bed6} ${prefix_all}
bed6_to_bedgraph.py ${chrom_info} ${bed6_uniq} ${prefix_uniq}

feature_l=$gene_l
bedgraph_to_depth.py ${watson_bg_all} 0 ${feature_l} | awk -v nf_all=${nf_all} 'BEGIN{ OFS="\t" } {$3=$3 * nf_all; print}' > ${watson_depth_all}
bedgraph_to_depth.py ${crick_bg_all} 0 ${feature_l}  | awk -v nf_all=${nf_all} 'BEGIN{ OFS="\t" } { $3=-$3 * nf_all; print }' > ${crick_depth_all}

bedgraph_to_depth.py ${watson_bg_uniq} 0 ${feature_l} | awk -v nf_uniq=${nf_uniq} 'BEGIN{OFS="\t"} { $3=$3 * nf_uniq; print }' > ${watson_depth_uniq}
bedgraph_to_depth.py ${crick_bg_uniq} 0 ${feature_l}  | awk -v nf_uniq=${nf_uniq} 'BEGIN{OFS="\t"} { $3=-$3 * nf_uniq; print }' > ${crick_depth_uniq}

combined_depth_all=${prefix_all}.norm.combined.depth
combined_depth_uniq=${prefix_uniq}.norm.combined.depth
paste ${watson_depth_all} ${crick_depth_all} | cut -f1,2,3,6 > ${combined_depth_all}
paste ${watson_depth_uniq} ${crick_depth_uniq} | cut -f1,2,3,6 > ${combined_depth_uniq}
plot_gene_signal.R ${combined_depth_all}  ${prefix_all}   ${prefix_all}.pdf
plot_gene_signal.R ${combined_depth_uniq} ${prefix_uniq}  ${prefix_uniq}.pdf
		     
