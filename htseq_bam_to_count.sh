#!/usr/bin/env bash

# Given a bam file and a GTF file, this script sorts the bam by name then run it using htseq
if [[ "$#" != 3 ]]; then
    echo "Usage: for dUTP RNA-seq libraries: htseq_bam_to_count.sh your.gtf your.unsorted.bam"
    echo "Usage: for other cases, specify the orientation you want: htseq_bam_to_count.sh your.gtf your.unsorted.bam [yes|reverse]"    
    exit 0
fi

FEATURE=gene_id
HTSEQ_MODE=intersection-strict
GTF=$1
BAM=$2

if [[ -z "$3" ]]; then
    htseq_strand_option=reverse
else
    htseq_strand_option=$3
fi

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
    htseq-count -f bam -r name -s ${htseq_strand_option} -t exon -i ${FEATURE} -m ${HTSEQ_MODE} ${BAM_SORTED} $GTF 2>${LOG} | sort -k1,1 >${OUTPUT}
    touch ${STATUS_HTSEQ_F}
fi
