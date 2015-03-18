#!/usr/bin/env bash

# Given a sorted BAM file, this script will output a bed12 file, and process the bed12 file to get a bedgraph file (with multimappers considered)
# Usage: bash run1-bam-to-bed12.sh ../2015-03-07/RSRZ.CC1.br1.tr1.ovary/genome_mapping/RSRZ.CC1.br1.tr1.ovary.x_rRNA.dm3_chengjian.sorted.bam test
END_TO_REVERSE_STRAND=1

chrom_info=extra.ChromInfo.txt
for i in RSRZ.CC1.br1.tr1.ovary RSRZ.GS1CT1.br1.tr1.ovary; do
    input=../2015-03-07/${i}/genome_mapping/${i}.x_rRNA.dm3_chengjian.sorted.bam
    output_prefix=${i}.constructs_only.all_mappers
    bed12=${output_prefix}.bed12
    bed6=${output_prefix}.bed6
    output_uniq_prefix=${i}.constructs_only.uniq_mappers
    bed6_uniq=${output_uniq_prefix}.bed6
    
    # bedtools bamtobed -bed12 -tag NH -i ${input} | awk -v strand=$END_TO_REVERSE_STRAND 'BEGIN{FS=OFS="\t"}{e=substr($4,length($4)); if (e==strand) $6=($6=="+"?"-":"+"); print $0; }' | grep -E 'GSV6|nos-gal4-vp16|UASp-EGFP' > ${bed12}
    for j in GSV6 nos-gal4-vp16 UASp-EGFP; do
	gl=$(grep ${j} ${chrom_info} | cut -f2)	
	n_mapper=$(sed -n '/genomie_mapper_reads/ s/genomie_mapper_reads:\t// p' ../2015-03-07/${i}/${i}.basic_stats)
	n_uniq_mapper=$(sed -n '/genomie_unique_mapper_reads/ s/genomie_unique_mapper_reads:\t// p' ../2015-03-07/${i}/${i}.basic_stats)
	bed12_gene=${i}.${j}.bed12
	cat ${bed12} | awk -v g=$j 'BEGIN{OFS="\t"} { if($1==g) {print} }' > ${bed12_gene}
	echo "bed12_to_gene_plots.sh ${bed12_gene} ${j} ${gl} ${n_mapper} ${n_uniq_mapper}"
    done
    # bedtools bed12tobed6 -i ${bed12} | grep -v 'track ' > ${bed6}
    # awk '{ if($5=="1") {print} }' ${bed6} > ${bed6_uniq}
    
    # bed6_to_bedgraph.py extra.ChromInfo.txt ${bed6} ${output_prefix}
    # bed6_to_bedgraph.py extra.ChromInfo.txt ${bed6_uniq} ${output_uniq_prefix}
    # echo ${n_mapper}
    # echo ${n_uniq_mapper}
    
    # nf=$( awk -v n_mapper=$n_mapper 'BEGIN{ print 1e6/n_mapper }' )
    # nf_uniq=$( awk -v n_uniq_mapper=$n_uniq_mapper 'BEGIN{ print 1e6/n_uniq_mapper }' )

    # watson_bg=${output_prefix}.Watson.bedGraph
    # crick_bg=${output_prefix}.Crick.bedGraph
    # for j in GSV6 nos-gal4-vp UASp-EGFP; do
    # 	feature_l=$(grep GSV6 *.txt | cut -f2)
    # 	watson_bg_one=${watson_bg%.Watson.bedGraph}.${j}.Watson.bedGraph
    # 	crick_bg_one=${crick_bg%.Crick.bedGraph}.${j}.Crick.bedGraph	
    # 	grep ${j} ${watson_bg} > ${watson_bg_one}
    # 	grep ${j} ${crick_bg}  > ${crick_bg_one}
    # 	feature_l=$(grep ${j} $chrom_info | cut -f2)
    # 	watson_depth=${watson_bg%.Watson.bedGraph}.${j}.Watson.depth
    # 	crick_depth=${crick_bg%.Crick.bedGraph}.${j}.Crick.depth
    # 	bedgraph_to_depth.py ${watson_bg_one} 0 ${feature_l} > ${watson_depth}
    # 	bedgraph_to_depth.py ${crick_bg_one} 0 ${feature_l} | awk 'BEGIN{OFS="\t"} {$3=-$3; print}' > ${crick_depth}
    # 	combined_depth=${watson_bg%.Watson.bedGraph}.${j}.combined.depth
    # 	paste ${watson_depth} ${crick_depth} | cut -f1,2,3,6 > ${combined_depth}
    # 	plot_gene_signal.R ${combined_depth} ${j}.all_mappers ${combined_depth}.pdf
    # done
done
