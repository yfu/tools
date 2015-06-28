#!/usr/bin/env bash

source /home/fuy2/repo/tools/piPipes_aux/global.sh

for D in SRS.GS1CT1.br1.tr1.ox.ovary SRS.GS1CT1.br2.tr1.ox.ovary SRS.GS1CT1.br3.tr1_2.ox.ovary SRS.GS1CT1.br3.tr1.ox.ovary SRS.GS1CT1.br3.tr2.ox.ovary; do
# for D in SRS.GS1CT1.br2.tr1.ox.ovary SRS.GS1CT1.br3.tr1_2.ox.ovary SRS.GS1CT1.br3.tr1.ox.ovary SRS.GS1CT1.br3.tr2.ox.ovary; do    
    # D=SRS.GS1CT1.br1.tr1.ox.ovary
    PREFIX=$D
    OUTPUT_D=${D}/gene+cluster+repBase
    mkdir -p ${OUTPUT_D}
    
    INPUT=${D}/input_read_files/${D}.x_rRNA.insert
    source ~/data/piPipes/common/dm3/variables
    
    CPU=8
    GENOME=genes+cluster+repBase
    GCR_IDX=~/data/piPipes/common/dm3_chengjian/BowtieIndex/gene+cluster+repBase

    GCR_ALLMAP_LOG=${OUTPUT_D}/${PREFIX}.gene+cluster+repBase.log
    GCR_ALLMAP_BED2=${OUTPUT_D}/${PREFIX}.gene+cluster+repBase.bed2
    echo "!!!" && echo ${GCR_ALLMAP_BED2}
    
    bowtie -r -v $genome_MM -a --best --strata -p $CPU \
       -S \
       ${GCR_IDX} \
       ${INPUT} \
       2> $GCR_ALLMAP_LOG | \
	samtools view -uS -F0x4 - 2>/dev/null | \
	${BEDTOOLS} bamtobed -i - > ${OUTPUT_D}/${PREFIX}.${GENOME}v${genome_MM}a.insert.bed && \
	${INSERT_BED_TO_BED2} $INPUT ${OUTPUT_D}/${PREFIX}.${GENOME}v${genome_MM}a.insert.bed > ${GCR_ALLMAP_BED2} && \
	rm -rf ${INSERT%.insert}.${GENOME}v${genome_MM}a.insert.bed
done
