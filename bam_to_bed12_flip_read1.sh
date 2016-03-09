#!/usr/bin/env bash

bam=$1
END_TO_REVERSE_STRAND=1		# For dUTP RNA-seq library
bedtools bamtobed -bed12 -tag NH -i ${bam} | awk -v strand=$END_TO_REVERSE_STRAND 'BEGIN{FS=OFS="\t"}{e=substr($4,length($4)); if (e==strand) $6=($6=="+"?"-":"+"); if ($5==1) print $0; }'
