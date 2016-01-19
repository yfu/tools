#!/usr/bin/env bash

# Given a bam file and a GTF file, this script sorts the bam by name then run it using htseq
if [[ "$#" != 2 ]]; then
    echo "Usage: htseq_bam_to_count.sh your.gtf your.unsorted.bam"
    exit 0
fi

FEATURE=gene_id
HTSEQ_MODE=intersection-strict
GTF=$1
BAM=$(readlink -e $2)
# PREFIX=$(basename ${BAM%.bam})
# PREFIX=${BAM%.bam}
PREFIX=$(basename ${BAM%.bam})
# Sorted by name required by htseq-count
BAM_SORTED=${PREFIX}.sorted_by_name.bam
STATUS_SORT_F=${PREFIX}.sorted_by_name.ok
LOG=${PREFIX}.${FEATURE}.${HTSEQ_MODE}.htseq.log
OUTPUT=${PREFIX}.${FEATURE}.${HTSEQ_MODE}.htseq
CPU=8

if [ -f ${STATUS_SORT_F} ]; then
    echo "Sorted bam file already exists."
else
    echo "Sorting the bam file by name..."
    samtools sort -n -m 2G -o ${BAM_SORTED} -O bam -T ${RANDOM} -@ ${CPU} ${BAM}
    touch ${STATUS_SORT_F}
fi

STATUS_HTSEQ_F=${PREFIX}.${FEATURE}.${HTSEQ_MODE}.ok
if [ -f ${STATUS_HTSEQ_F} ]; then
    echo "htseq-count has already been run."
else
    echo "Running htseq-count..."
    htseq-count -f bam -r name -s reverse -t exon -i ${FEATURE} -m ${HTSEQ_MODE} ${BAM_SORTED} $GTF 2>${LOG} | sort -k1,1 >${OUTPUT}
    touch ${STATUS_HTSEQ_F}
fi
