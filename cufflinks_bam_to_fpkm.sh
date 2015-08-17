#!/usr/bin/env bash

REPO=/home/fuy2/repo/tools
source ${REPO}/piPipes_aux/functions.sh
COUNTER=${REPO}/htseq_bam_to_count.sh

while getopts "hd:f:i:c:o:g:B:xvLD" OPTION; do
    case $OPTION in
	f)BAM=${OPTARG};; # the directory of pipeline results
	c)CPU=$OPTARG ;;
	g)GENOME=$OPTARG ;;
    esac
done
PREFIX=${BAM%.bam}
PREFIX=${PREFIX%.sorted}
PREFIX=$(basename ${PREFIX})

if [[ -z "$CPU" ]]; then
    CPU=8
fi

if [[ "$GENOME" == "mm10" ]]; then
    TRANSCRIPTOME_GTF=/home/fuy2/data/shared/mm10/gencode.vM4.corrected_chrom_names.gtf
    GENOME_FA=/home/fuy2/data/shared/mm10/mm10.fa
elif [[ "$GENOME" == "mm10_small" ]]; then
    TRANSCRIPTOME_GTF=/data/fuy2/shared/mm10_small/mm10_small.gtf
    GENOME_FA=/data/fuy2/shared/mm10_small/mm10_small.fa
else
    echo2 "Unrecognized genome!"
    exit 1
fi

if [[ -z "$BAM" ]]; then
   echo2 "Please provide a bam file!"
   exit 1
fi
LIBRARY_TYPE=fr-firststrand
NO_LEN_CORRECTION=""
STATUS_F=.${PREFIX}.cufflinks.quantification_by_cuff_compatible_hits_norm
CUFFLINKS_DIR=${PREFIX}_cufflinks
mkdir -p ${CUFFLINKS_DIR}
[ ! -f ${STATUS_F} ] && \
    cufflinks \
	-o $CUFFLINKS_DIR \
	-p $CPU \
	$NO_LEN_CORRECTION \
	-G ${TRANSCRIPTOME_GTF} \
	-b $GENOME_FA \
	-u \
	--library-type $LIBRARY_TYPE \
	--compatible-hits-norm \
	--no-update-check \
	${BAM} \
	2> $CUFFLINKS_DIR/${PREFIX}.cufflinks.log && \
    touch .${JOBUID}.status.${STEP}.quantification_by_cuff_compatible_hits_norm
MapMass=`grep 'Normalized Map Mass' $CUFFLINKS_DIR/${PREFIX}.cufflinks.log | cut -d' ' -f4`
export NormScale=`echo $MapMass | awk '{printf "%f",1000000.0/$1}'`
echo -ne "$NormScale" > .${PREFIX}.cufflinks_nf
