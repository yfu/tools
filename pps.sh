#!/bin/bash -x

# pps: a Ping-Pong score calculator
# Accept bed2 as input
# Yu Fu (yfu at yfu dot me)

ext5=30;
ext3=30
# awk -v ext=$ext 'BEGIN{OFS="\t" } { $2-=ext; if($2<0) {$2=0}; print  }' r1.no_adaptor_hi5_unox.1mm.bed2 | head

# TODO: replace bedtools slop with with my own implementation, which is able to tell:
# which reads are extended to contain N's and 
# which reads reached the borders of scaffolds.

# For simplicity, the intermediate files have fixed names (like STAR output).
input=$1

# If the genome length file does not exist, I'll create one.
if [ ! -f genome.length ]; then 
    python get_fasta_seq_length.py genome.fa > genome.length
fi

bedtools slop -s -l ${ext5} -r ${ext3} -i ${input} -g genome.length > ext.tmp.bed2
# Hack the name column of the bed file in order to keep passing the information to the sam file to make my life easier.
# The name field format: (original chrom)_(original start)_(orginal end)_(original strand)_(ext chrom)_(ext start)_(ext end)_(ext strand)_(copy #)_(ntm)
paste ext.tmp.bed2 ${input} | awk '{ print $1 "\t" $2 "\t" $3 "\t" "orig_" $8 "_" $9 "_" $10 "_" $13 "_ext_" $1 "_" $2 "_" $3 "_" $6 "_" $4 "_" $5 "\t" 0 "\t" $6 "\t" $7 }' > ext.bed2 # && rm ext.tmp.bed2

bedtools getfasta -fi genome.fa -bed ext.bed2 -s -name -fo ext.fa
if [ ! -f ext.idx.1.ebwt ]; then
    bowtie-build ext.fa ext.idx 1>/dev/null
fi

cat ${input} | awk 'BEGIN{u="_"} { print ">" $1 u $2 u $3 u $4 u $5 u $6 "\t"; print $7}' > orig_reads

# TODO: remove -m 100 to see if this affects the final result.
bowtie -p 30 -f -v 0 -a -m 100 --best --strata -S ext.idx orig_reads | samtools view -Sb - -F0x4 > mapping_result.bam
samtools sort -@ 20 -m 2G mapping_result.bam  mapping_result.sorted && rm mapping_result.bam
samtools view -b -h -f 0x10 mapping_result.sorted.bam > mapping_result.reverse_strand.sorted.bam

# Get NTM for each reads
# samtools view mapping_result.reverse_strand.sorted.bam | awk '{ print $10 "\t" $3 }' | sort -k1 | uniq | cut -f1 | uniq -c | sed -r 's/\s+([0-9]+)\s+([ATCGN]+)/\2\t\1/' > ntm_reads_mapped_to_reverse_index.txt
samtools view mapping_result.reverse_strand.sorted.bam > mapping_result.reverse_strand.noheader.sorted.sam

# samtools view mapping_result.reverse_strand.sorted.bam | python pp_hist.py ntm_reads_mapped_to_reverse_index.txt
