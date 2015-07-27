#!/usr/bin/env bash

for D in SRS.GS1CT1.br1.tr1.ox.ovary SRS.GS1CT1.br2.tr1.ox.ovary SRS.GS1CT1.br3.tr1_2.ox.ovary SRS.GS1CT1.br3.tr1.ox.ovary SRS.GS1CT1.br3.tr2.ox.ovary; do
# for D in SRS.GS1CT1.br2.tr1.ox.ovary; do
    # for D in SRS.GS1CT1.br2.tr1.ox.ovary SRS.GS1CT1.br3.tr1_2.ox.ovary SRS.GS1CT1.br3.tr1.ox.ovary SRS.GS1CT1.br3.tr2.ox.ovary; do
    # D=SRS.GS1CT1.br1.tr1.ox.ovary
    PREFIX=$D
    GCR_D=${D}/gene+cluster+repBase
    GCR_ALLMAP_BED2=${GCR_D}/${PREFIX}.gene+cluster+repBase.bed2
    PREFIX_COUNT=${GCR_D}/${PREFIX}.gene+cluster+repBase.bed2
    TABLE=${D}/${PREFIX}.table

    # PREFIX_MAPPER=SRS.GS1CT1.br1.tr1.ox.ovary.x_rRNA.x_hairpin.dm3_chengjianv0
    # D=../../SRS.GS1CT1.br1.tr1.ox.ovary
    PREFIX_MAPPER=${D}/genome_mapping/${PREFIX}.x_rRNA.x_hairpin.dm3_chengjianv0
    ALL_MAPPER=${PREFIX_MAPPER}.all.bed2
    SIRNA_MAPPER=${PREFIX_MAPPER}.all.siRNA.bed2
    
    MIRNA_MAPPER="TODO"
    
    total_reads=$(awk '{ s+=$4/$5 } END{ print s }' ${ALL_MAPPER})
    echo -e "Total\t${total_reads}" > ${TABLE}
    uniq_reads=$(awk '{ if($5==1) { s+=$4 } } END{ print s }' ${ALL_MAPPER})
    echo -e "Total unique reads\t${uniq_reads}" >> ${TABLE}
    
    # Useful parameters, such as lengths of piRNAs and siRNAs
    # TODO: maybe use *.all.siRNA.bed2?
    source ~/repo/tools/piPipes_aux/dm3.parameters.sh
    CISNAT=/data/fuy2/piPipes/common/dm3/cisNATs.bed.gz
    cisnat_sirna_reads=$(cat ${SIRNA_MAPPER} | bedtools intersect -u -a - -b ${CISNAT} | awk '{ s+=$4/$5 } END{ print s }')
    echo -e "cis-NAT\t${cisnat_sirna_reads}" >> ${TABLE}

    # 42AB unique mappers
    # TODO: consider making the piPipes to generate this file for me...
    CLUSTER_42AB=/data/fuy2/piPipes/common/dm3/Brennecke.piRNAcluster.42AB.bed6.gz
    cluster_42ab_reads=$(cat ${ALL_MAPPER} | awk '{ if($5==1) { print } }' | bedtools intersect -u -a - -b ${CLUSTER_42AB} | awk '{ s+=$4/$5 } END{ print s }')
    echo -e "42AB unique mappers\t${cluster_42ab_reads}" >> ${TABLE}
    
    # flam unique mappers
    # TODO: consider making the piPipes to generate this file for me...
    CLUSTER_FLAM=/data/fuy2/piPipes/common/dm3/Brennecke.piRNAcluster.flam.bed6.gz
    cluster_flam_reads=$(cat ${ALL_MAPPER} | awk '{ if($5==1) { print } }' | bedtools intersect -u -a - -b ${CLUSTER_FLAM} | awk '{ s+=$4/$5 } END{ print s }')
    echo -e "flam unique mappers\t${cluster_flam_reads}" >> ${TABLE}

    # non-transposon mappers
    GCR_D=${D}/gene+cluster+repBase
    GCR_BED2=${GCR_D}/${PREFIX}.gene+cluster+repBase.bed2
    transposon_sirna_reads=$(cat ${GCR_BED2} | awk -v top=$siRNA_top -v bot=$siRNA_bot '{ l=$3-$2; if( l>=bot && l <= top && $1 ~ /FBgn/) { s+=$4/$5 } } END{ print s }')
    sirna_reads=$( awk '{ s+=$4/$5 } END{ print s }' ${SIRNA_MAPPER} )
    non_transposon_sirna_reads=$( cat ${GCR_BED2} | awk -v sirna_reads=${sirna_reads} -v transposon_sirna_reads=${transposon_sirna_reads} 'BEGIN{ print sirna_reads-transposon_sirna_reads }')
    echo -e "non-transposon siRNA reads\t${non_transposon_sirna_reads}" >> ${TABLE}
done
