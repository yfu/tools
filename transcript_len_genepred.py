#!/usr/bin/env python

# Output transcript lengths together with gene names given a file in extended genePred format
# gtfToGenePred -genePredExt -geneNameAsName2 ~/data/piPipes/common/mm9/mm9.genes.gtf stdout

# Usage: cat mm9.genes.genePred | transcript_len_genepred.py > mm9.transcript.length
# Next step: awk 'BEGIN{OFS=FS="\t"} { c[$1]+=1; t[$1]+=$3 } END{ for(i in c) {print i, t[i]/c[i]} }' mm9.transcript.length > mm9.gene.avg_len
import sys
import fileinput

transcripts_len = {}
transcripts_to_genes = {}

for line in fileinput.input():
    line = line.strip()
    col = line.split()
    transcript_id = col[0]
    exon_starts = col[8].split(",")
    exon_ends = col[9].split(",")
    gene_name = col[11]
    transcripts_to_genes[transcript_id] = gene_name
    # print transcripts_len.keys()
    if transcript_id in transcripts_len:
        print >> sys.stderr, "Duplicate transcript id found in genePred file: " + line
    for i in range(len(exon_starts)):
        if exon_starts[i] == "" or exon_ends[i] == "":
            continue
        if transcript_id not in transcripts_len:
            transcripts_len[transcript_id] = int(exon_ends[i]) - int(exon_starts[i])
        else:
            transcripts_len[transcript_id] += int(exon_ends[i]) - int(exon_starts[i])

for t in transcripts_len:
    print transcripts_to_genes[t] + "\t" + t + "\t" + str(transcripts_len[t])
