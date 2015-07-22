#!/usr/bin/env bash

###############################################################################
# Part 1: mapping reads to the genome and ignore them for downstream analysis #
###############################################################################
# fq=Mattick.RNaseR.HeLa.BP.rep1.fq.gz
FQ=$1
ln -s $FQ ./
PREFIX=${FQ%.fq.gz}
PREFIX=${FQ%.fq}
PREFIX=$(basename $PREFIX)
# zcat ${fq} | wc -l | awk '{ print $0/4 }' > ${PREFIX}.total_reads

# head -n1000000 Mattick.CSQ.K562.BP.rep1.fq > Mattick.CSQ.K562.BP.rep1.1000000.fq
# The same with bp.sh

GENOME=hg19
CPU=16
GTF=/data/fuy2/shared/${GENOME}/${GENOME}.gencode.gtf
GENOME_FA=/data/fuy2/shared/${GENOME}/${GENOME}.fa
PIPELINE_DIR=~/repo/tools/bp

OUTPUT_DIR=.
COMMON_DIR=${PIPELINE_DIR}/common/${GENOME}

# FQ=Mattick.CSQ.K562.BP.rep1.1000000.fq

if [[ ${FQ##*.} == gz ]]; then
    star_read_file_command='--readFilesCommand gunzip -c '
else
    star_read_file_command=''
fi

OUTPUT_PREFIX=${PREFIX}.${GENOME}.
# Set --genomeLoad to NoSharedMemory if you want to run it in the cluster
STAR \
    --runMode alignReads \
    --genomeDir ${COMMON_DIR}/STARIndex \
    --readFilesIn ${FQ}\
    ${star_read_file_command} \
    --runThreadN $CPU \
    --outFilterScoreMin 0 \
    --outFilterScoreMinOverLread 0.72 \
    --outFilterMatchNmin 0 \
    --outFilterMatchNminOverLread 0.72 \
    --outFilterMultimapScoreRange 1 \
    --outFilterMultimapNmax -1 \
    --outFilterMismatchNmax 10 \
    --outFilterMismatchNoverLmax 0.05 \
    --alignIntronMax 0 \
    --alignIntronMin 21 \
    --outFilterIntronMotifs None \
    --genomeLoad NoSharedMemory \
    --outFileNamePrefix $OUTPUT_PREFIX \
    --outSAMunmapped None \
    --outReadsUnmapped Fastx \
    --outSJfilterReads Unique \
    --seedSearchStartLmax 20 \
    --seedSearchStartLmaxOverLread 1.0 \
    --chimSegmentMin 0 2>&1 1> ${OUTPUT_PREFIX}STAR.log

# Make bigWig files from the regular RNA-seq reads
# Normalization factor = 1
bam_to_bigwig.sh ~/data/shared/hg19/hg19.ChromInfo.txt test.100000.hg19.Aligned.out.sam 1

###########################################################################
# Part 2: mapping reads to 5' introns to extract useful part of the reads #
###########################################################################

input=${PREFIX}.${GENOME}.Unmapped.out.mate1
prefix=${input%.hg19.Unmapped.out.mate1}
cpu=16

# Keep the original reads or the reverse complements
# seqtk seq -r ${input} > ${prefix}.hg19.Unmapped.out.mate1.rc
# input=${prefix}.hg19.Unmapped.out.mate1.rc
# prefix=${prefix}.rc

# First 23nt of introns
index_5pintron=~/repo/tools/bp/common/hg19/hg19.gencode.intron.uniq.23nt
# Whole introns
index_intron=~/repo/tools/bp/common/hg19/hg19.gencode.intron.uniq
# prefix=Mattick.CSQ.K562.BP.rep1.1000000
# input=${prefix}.hg19.Unmapped.out.mate1
map_to_5_intron=${prefix}.hg19.map_to_5intron.bam
map_to_5_intron_s=${prefix}.hg19.map_to_5intron.s.sam
map_to_5_intron_as=${prefix}.hg19.map_to_5intron.as.sam
# map_to_5_intron_fa=${prefix}.lariat_and_5intron.s.fa

log=${prefix}.hg19.map_to_5intron.log
# Require that MAPQ >= 10
bowtie2 --local --score-min L,45,0 -D 20 -R 2 -N 0 -L 20 -i L,1,0 --phred33 -p ${cpu} -x ${index_5pintron} -U ${input} 2> ${log} | samtools view -q 10 -bS -F 0x4 - | samtools sort -@ ${cpu} -T ${RANDOM} -O bam > ${map_to_5_intron}
samtools index ${map_to_5_intron}
samtools view -f 16 ${map_to_5_intron} > ${map_to_5_intron_as}
samtools view -F 16 ${map_to_5_intron} > ${map_to_5_intron_s}
# CONSIDER how to use those mapping to the antisense strand of intons (notice that all the intron indexes are in the sense strand)
# gawk -v fpl=23 -v tl=20 '{ pat="(^[0-9]+)S" fpl "M"; match($6, pat, ary); a=ary[1]; if(a>=tl) { print ">" $3 "|" substr($10, a+1, fpl); print substr($10, 1, a); }}' ${map_to_5_intron_s} > ${map_to_5_intron_fa}
echo "Mapping to 23nt introns:" > ${log}

# Output the softclipped part on the left end of reads, aligned part, and softclipped part on the right end of reads (always on the reference strand)
~/repo/tools/parse_bam/get_before_alignment_after.py -s -f ${map_to_5_intron} 2>> ${log} | awk -v p=${prefix} 'BEGIN{bef=p".before.fa"; ali=p".align.fa"; aft=p".after.fa";} { if($2!="") {print ">" $1 "\n" $2> bef}; if($3!="") {print ">" $1 "\n" $3 > ali}; if($4!="") {print ">" $1 "\n" $4 > aft}  }'

# Use those reads that map to the antisense strand of introns
# fpl: five prime length (23) requiring part of the read to match the 5' end of an intron
# tl: After trimming, there should be at least tl (20) nt left on the reads
# gawk -v fpl=23 -v tl=20 '{ pat="(^[0-9]+)S" fpl "M"; match($6, pat, ary); a=ary[1]; if(a>=tl) { print ">" $3 "|" substr($10, a, fpl); print substr($10, 1, a); }}' ${map_to_5_intron_c} >> ${map_to_5_intron_fa}
map_to_5_intron_after=${prefix}.after.fa
map_to_intron=${prefix}.lariat_and_5intron.map_to_intron.bam
map_to_intron_log=${prefix}.lariat_and_5intron.map_to_intron.log
bowtie2 -x ${index_intron} -f -U ${map_to_5_intron_after} -p ${cpu} 2> ${map_to_intron_log} | samtools view -bS -F 0x4 -q 10 - > ${map_to_intron}

# Get the weblogo
p=${prefix}.lariat_and_5intron.map_to_intron
bedtools bamtobed -i ${p}.bam | bedtools getfasta -s -fi ~/repo/tools/bp/common/hg19/hg19.gencode.intron.uniq.fa -bed - -fo - | awk '{ if(NR%2==1) {print $0} else { print substr($0, length($0)-10+1) } }' | weblogo -D fasta -F pdf > ${p}.last10.logo.pdf

