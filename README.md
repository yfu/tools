tools
=====

## pirna_cluster_caller.py

A simple script that calls piRNA clusters given a bed2 file and a file containing the length of each chromosome.

Usage: pirna_cluster_caller.py -f FILE -g GENOME [-u]

## columnize_fasta.sh

Given a fasta file, this script outputs two columns, which is suitable for line-oriented programs.

## splice_junction_bed.py

Given a bed12 format file, this script will output the bed12 describing left and right flanking regions of the splice site.

## pps_plot.R
Given the scores from pps_simple.py, this script will plot them.

## bed12_to_bedgraph.sh

Given a bed12 file for just one gene, this script outputs the gene plot (for all_mappers, uniq_mappers)

Usage: bed12_to_depth.sh one_gene.bed12 gene.name gene.length n_all_mappers n_uniq_mapper                           