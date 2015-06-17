#!/usr/bin/env bash

# Given a bam file and a GTF file, this script sorts the bam by name then run it using htseq
FEATURE=gene_id
BAM=$1
PREFIX=$(basename ${BAM%.bam})
# Sorted by name required by htseq-count
BAM_SORTED=${PREFIX}.sorted.bam
STATUS_SORT_F=${PREFIX}.sorted_by_name.ok
GTF=$2
LOG=${PREFIX}.${FEATURE}.${HTSEQ_MODE}.htseq.log
HTSEQ_MODE=union
OUTPUT=${PREFIX}.${FEATURE}.${HTSEQ_MODE}.htseq

if [ -f ${STATUS_SORT_F} ]; then
    echo "Sorted bam file already exists."
else
    echo "Sorting the bam file by name..."
    samtools sort -n -m 2G -o ${BAM_SORTED} -O bam -T ${RANDOM} -@ 8 ${BAM}
    touch ${STATUS_SORT_F}
fi

STATUS_HTSEQ_F=${PREFIX}.${FEATURE}.${HTSEQ_MODE}.ok
if [ -f ${STATUS_HTSEQ_F} ]; then
    echo "htseq-count has already been run."
else
    echo "Running htseq-count..."
    htseq-count -f bam -r name -s reverse -t exon -i gene_id -m ${HTSEQ_MODE} ${BAM_SORTED} $GTF >${OUTPUT} 2>${LOG}
    touch ${STATUS_HTSEQ_F}
fi
