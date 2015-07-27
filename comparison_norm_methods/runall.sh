#!/usr/bin/env bash

# Given a piPipes output dir as the input
source /home/fuy2/repo/tools/piPipes_aux/global.sh
D=$1
D=~/data/gfp_pirna/test/test_data/SRS.CC1.Ago3IP.br1.tr1.ox.ovary
# D=SRS.GS1CT1.br1.tr1.ox.ovary
PREFIX=${D##*/}
OUTPUT_D=${PREFIX}/gene+cluster+repBase
mkdir -p ${OUTPUT_D}

INPUT=${D}/input_read_files/${PREFIX}.x_rRNA.insert
source ~/data/piPipes/common/dm3/variables

CPU=8
GENOME=genes+cluster+repBase
GCR_IDX=~/data/piPipes/common/dm3_chengjian/BowtieIndex/gene+cluster+repBase

GCR_ALLMAP_LOG=${OUTPUT_D}/${PREFIX}.gene+cluster+repBase.log
GCR_ALLMAP_BED2=${OUTPUT_D}/${PREFIX}.gene+cluster+repBase.bed2
## echo "!!!" && echo ${GCR_ALLMAP_BED2}

bowtie -r -v $genome_MM -a --best --strata -p $CPU \
       -S \
       ${GCR_IDX} \
       ${INPUT} \
       2> $GCR_ALLMAP_LOG | \
	samtools view -uS -F0x4 - 2>/dev/null | \
	${BEDTOOLS} bamtobed -i - > ${OUTPUT_D}/${PREFIX}.${GENOME}v${genome_MM}a.insert.bed && \
	${INSERT_BED_TO_BED2} $INPUT ${OUTPUT_D}/${PREFIX}.${GENOME}v${genome_MM}a.insert.bed > ${GCR_ALLMAP_BED2} && \
	rm -rf ${INSERT%.insert}.${GENOME}v${genome_MM}a.insert.bed

GCR_D=${PREFIX}/gene+cluster+repBase
GCR_ALLMAP_BED2=${GCR_D}/${PREFIX}.gene+cluster+repBase.bed2
PREFIX_COUNT=${GCR_D}/${PREFIX}.gene+cluster+repBase.bed2

count_cmd=${PREFIX}/count.${RANDOM}${RANDOM}.commands
echo "awk '{ if(\$1 ~ /cluster|42AB|flam/) {c[\$1]+=\$4/\$5} } END{ for(i in c) { print i \"\t\" c[i] } }' ${GCR_ALLMAP_BED2} | sort -k1,1 > ${PREFIX_COUNT}.cluster.count" > ${count_cmd}
echo "awk '{ if(\$1 ~ /^FBgn/) {c[\$1]+=\$4/\$5} } END{ for(i in c) { print i \"\t\" c[i] } }' ${GCR_ALLMAP_BED2} | sort -k1,1 > ${PREFIX_COUNT}.repBase.count" >> ${count_cmd}
echo "awk '{ if(\$1~/^NM_|^NR_/) {c[\$1]+=\$4/\$5} } END{ for(i in c) { print i \"\t\" c[i] } }' ${GCR_ALLMAP_BED2} | sort -k1,1 > ${PREFIX_COUNT}.gene.count"  >> ${count_cmd}
echo "awk '{if(\$1~/nos-gal4-vp16|GSV6|UASp-EGFP/) {c[\$1]+=\$4/\$5} } END{ for(i in c) { print i \"\t\" c[i] } }' ${GCR_ALLMAP_BED2} | sort -k1,1 > ${PREFIX_COUNT}.construct.count" >> ${count_cmd}
cat ${count_cmd} | parallel --progress 
function get_normfactor() {
    D=$1
    PREFIX=${D##*/}
    # GCR_D=${PREFIX}/gene+cluster+repBase
    GCR_ALLMAP_BED2=${GCR_D}/${PREFIX}.gene+cluster+repBase.bed2
    PREFIX_COUNT=${GCR_D}/${PREFIX}.gene+cluster+repBase.bed2
    TABLE=${PREFIX}/${PREFIX}.table

    # PREFIX_MAPPER=${D}/genome_mapping/${PREFIX}.x_rRNA.x_hairpin.dm3_chengjianv0
    ALL_MAPPER=${D}/genome_mapping/*.x_rRNA.x_hairpin.*.all.bed2
    SIRNA_MAPPER=${D}/genome_mapping/*.x_rRNA.x_hairpin.*.all.siRNA.bed2
    
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
    CLUSTER_42AB=/data/fuy2/piPipes/common/dm3/Brennecke.piRNAcluster.42AB.bed6.gz
    cluster_42ab_reads=$(cat ${ALL_MAPPER} | awk '{ if($5==1) { print } }' | bedtools intersect -u -a - -b ${CLUSTER_42AB} | awk '{ s+=$4/$5 } END{ print s }')
    echo -e "42AB unique mappers\t${cluster_42ab_reads}" >> ${TABLE}    
    # flam unique mappers
    CLUSTER_FLAM=/data/fuy2/piPipes/common/dm3/Brennecke.piRNAcluster.flam.bed6.gz
    cluster_flam_reads=$(cat ${ALL_MAPPER} | awk '{ if($5==1) { print } }' | bedtools intersect -u -a - -b ${CLUSTER_FLAM} | awk '{ s+=$4/$5 } END{ print s }')
    echo -e "flam unique mappers\t${cluster_flam_reads}" >> ${TABLE}
    # non-transposon mappers
    GCR_D=${PREFIX}/gene+cluster+repBase
    GCR_BED2=${GCR_D}/${PREFIX}.gene+cluster+repBase.bed2
    transposon_sirna_reads=$(cat ${GCR_BED2} | awk -v top=$siRNA_top -v bot=$siRNA_bot '{ l=$3-$2; if( l>=bot && l <= top && $1 ~ /FBgn/) { s+=$4/$5 } } END{ printf s }')
    echo "**" ${transposon_sirna_reads} "**"
    sirna_reads=$( awk '{ s+=$4/$5 } END{ print s }' ${SIRNA_MAPPER} )
    non_transposon_sirna_reads=$( cat ${GCR_BED2} | awk -v sirna_reads=${sirna_reads} -v transposon_sirna_reads=${transposon_sirna_reads} 'BEGIN{ print sirna_reads-transposon_sirna_reads }')
    echo -e "non-transposon siRNA reads\t${non_transposon_sirna_reads}" >> ${TABLE}
}

export -f get_normfactor
get_normfactor $D

