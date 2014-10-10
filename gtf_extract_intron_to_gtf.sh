#!/usr/bin/env bash

# Given a GTF file from stdin, for example dm3.genes.gtf, extract and output the introns in GTF format

gtfToGenePred stdin stdout -genePredExt | cut -f1-12 | awk '{if($8>=2) print}' | awk 'BEGIN{OFS="\t"} { A=""; B=""; split($9, a, ","); split($10, b, ","); for(i=1;i!=length(a)-1;i++) {A=A""b[i]","; B=B""a[i+1]","} $9=A; $10=B; $8=$8-1; $4=$6=b[1]; $5=$7=a[2]; print }' | genePredToGtf file stdin stdout | awk 'BEGIN{FS="\t"; OFS="\t"} { if($3="exon") {$3="intron"; print} }' | sed 's/exon_number/intron_number/' | sed 's/exon_id/intron_id/'
