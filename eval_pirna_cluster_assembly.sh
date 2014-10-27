#!/usr/bin/env bash

assembly_fa=$1
bedtools getfasta -fi ~/data/piPipes/common/dm3/dm3.fa -bed ~/data/piPipes/common/dm3/Brennecke.piRNAcluster.bed6.gz -fo - -s | fold -w72 > dmel.pirna.clusters.fa

makeblastdb -dbtype nucl -in dmel.pirna.clusters.fa 
blastn -db dmel.pirna.clusters.fa -outfmt 6 -num_threads 32 -query $assembly_fa -evalue 1e-20 > blast_transcriptome_against_clusters.blastn.table.txt

cat blast_transcriptome_against_clusters.blastn.table.txt | grep 'chr2R:2144348' | awk '{ if($3>95) print }' | cut -f2,9,10 | awk 'BEGIN{OFS="\t"} { if($2>$3) {a=$2; $2=$3; $3=a}; $1="42AB"; print }' | sort -k2,2n | bedtools merge -i - > 42AB.transcripts.bed

cat 42AB.transcripts.bed | awk '{OFS="\t"} {$1="chr2R"; $2=$2 + 2144348; $3=$3+2144348; print}' > 42AB.transcripts.genome_coor.bed 
