#!/usr/bin/env bash
# This is for the unstranded GTEx data
source /home/fuy2/repo/tools/piPipes_aux/global.sh
source /home/fuy2/repo/tools/piPipes_aux/functions.sh

while getopts "hd:i:c:o:g:B:xvLD" OPTION; do
    case $OPTION in
	h)usage && exit 0 ;;
	d)PR=`readlink -f ${OPTARG}`;; # the directory of pipeline results
	o)OUTDIR=`readlink -f ${OPTARG}` ;;
	c)CPU=$OPTARG ;;
    esac
done
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8

GENOME=hg19
COMMON_FOLDER=~/data/piPipes/common/${GENOME}
CHROM=${COMMON_FOLDER}/${GENOME}.ChromInfo.txt
[ ! -z ${OUTDIR} ] || OUTDIR=$PWD # if -o is not specified, use current directory

xrRNA_LEFT_FQ=$(readlink -e ${PR}/input_read_files/*.x_rRNA.1.fq)
xrRNA_RIGHT_FQ=$(readlink -e ${PR}/input_read_files/*.x_rRNA.2.fq)
PREFIX=$(basename $PR)
PREFIX=${PREFIX%.x_rRNA.1.fq}

CUFFLINKS_DIR=${PR}/cufflinks_output
GENOMIC_MAPPING_DIR=${PR}/genome_mapping
MapMass=$(grep 'Normalized Map Mass' $CUFFLINKS_DIR/${PREFIX}.cufflinks.log | cut -d' ' -f4)
export NormScale=$(echo $MapMass | awk '{printf "%f",1000000.0/$1}')
echo $NormScale

mkdir -p ${OUTDIR} && cd ${OUTDIR}
BW_OUTDIR=$OUTDIR

echo2 "Making bigWig from sorted bam"
echo "${BEDTOOLS} bamtobed -bed12 -tag NH -i ${GENOMIC_MAPPING_DIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.bam > ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.bed12 && \
    $bedtools_piPipes genomecov -scale $NormScale -split -bg -i ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.bed12 -g $CHROM > ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.ns.bedGraph && \
    bedGraphToBigWig ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.ns.bedGraph $CHROM ${BW_OUTDIR}/${PREFIX}.x_rRNA.${GENOME}.sorted.unique.ns.bigWig"

