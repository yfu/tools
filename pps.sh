#!/bin/bash -x

# pps: a Ping-Pong score calculator
# Accept bed2 as input
# Yu Fu (yfu at yfu dot me)

# Usage: put a genome.fa containing the genome sequence and a bed2 file containing the alignment. 
# A Bed2 file is just like a normal bed file, but it has the # copies as the 4th column and the NTM as the 5th column
# TODO: include a script (bed_2_bed2.py) that, given the bed file (converted from sam file), outputs a bed2 file.

ext5=30;
ext3=30
# This applied to a single instance of bowtie only. Always set this to one.
cpu=1
# $thread is used to specify the number of threads used when this script splits the fasta file and call bowtie
thread=8

# awk -v ext=$ext 'BEGIN{OFS="\t" } { $2-=ext; if($2<0) {$2=0}; print  }' r1.no_adaptor_hi5_unox.1mm.bed2 | head

# TODO: replace bedtools slop with with my own implementation, which is able to tell:
# which reads are extended to contain N's and 
# which reads reached the borders of scaffolds.

# For simplicity, the intermediate files have fixed names (like STAR output).
input=$1

# If the genome length file does not exist, I'll create one.
if [ ! -f genome.length ]; then 
    get_fasta_seq_length.py genome.fa > genome.length
fi

bedtools slop -s -l ${ext5} -r ${ext3} -i ${input} -g genome.length > ext.tmp.bed2
# Hack the name column of the bed file in order to keep passing the information to the sam file to make my life easier.
# The name field format: (original chrom)_(original start)_(orginal end)_(original strand)_(ext chrom)_(ext start)_(ext end)_(ext strand)_(copy #)_(ntm)
paste ext.tmp.bed2 ${input} | awk '{ print $1 "\t" $2 "\t" $3 "\t" "orig_" $8 "_" $9 "_" $10 "_" $13 "_ext_" $1 "_" $2 "_" $3 "_" $6 "_" $4 "_" $5 "\t" 0 "\t" $6 "\t" $7 }' > ext.bed2 && rm ext.tmp.bed2

bedtools getfasta -fi genome.fa -bed ext.bed2 -s -name -fo ext.fa

# Split the fasta file into multiple files so that building index is faster and and bowtie will not complain about the size of the fasta file
fasta_splitter.py -p $thread -f ext.fa

# Instead of the loop, I am using GNU Parallel to parallelize the bowtie building process
# for i in $(seq 1 $thread); do
#     idx_name=ext.idx.p${i}
#     if [ ! -f ${idx_name}.1.ebwt ]; then
# 	bowtie-build ext.fa.${i} ext.idx.p${i} 1>ext.idx.p${i}.buildlog
#     fi
# done

echo "Building bowtie index(es)"
parallel -j ${thread} 'bowtie-build ext.fa.{} ext.idx.p{} 1>ext.idx.p{}.buildlog' ::: $(seq 1 $thread)

cat ${input} | awk 'BEGIN{u="_"} { print ">" $1 u $2 u $3 u $4 u $5 u $6 "\t"; print $7}' > orig_reads

# TODO: remove -m 100 to see if this affects the final result.
map_sort_convert() {
    i=$1
    bowtie -p 1 -f -v 0 -a -m 100 --best --strata -S ext.idx.p${i} orig_reads | samtools view -u -Sb - -F0x4 > mapping_result.p${i}.bam
    samtools sort -@ 1 -m 5G -l 0 mapping_result.p${i}.bam mapping_result.p${i}.sorted && rm mapping_result.p${i}.bam
    samtools view -h -f 0x10 mapping_result.p${i}.sorted.bam > mapping_result.reverse_strand.p${i}.sorted.sam
    rm mapping_result.p${i}.sorted.bam
}
export -f map_sort_convert
parallel -j ${thread} 'map_sort_convert {}' ::: $(seq 1 $thread)

cat $(ls mapping_result.reverse_strand.p*.sorted.sam) | grep -v '@' > mapping_result.reverse_strand.noheader.sorted.sam
# samtools view -S -h mapping_result.reverse_strand.sorted.sam > mapping_result.reverse_strand.noheader.sorted.sam
pp_hist.py mapping_result.reverse_strand.noheader.sorted.sam >pp.log 2>pp.err

# Get NTM for each reads
# samtools view mapping_result.reverse_strand.sorted.bam | awk '{ print $10 "\t" $3 }' | sort -k1 | uniq | cut -f1 | uniq -c | sed -r 's/\s+([0-9]+)\s+([ATCGN]+)/\2\t\1/' > ntm_reads_mapped_to_reverse_index.txt

grep -E "[0-9]+\s[0-9]" pp.log > pp.dataframe
pingpong_zscore.R ${input} pp.dataframe 
# samtools view mapping_result.reverse_strand.sorted.bam | python pp_hist.py ntm_reads_mapped_to_reverse_index.txt
