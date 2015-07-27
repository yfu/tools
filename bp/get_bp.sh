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
bam_to_bigwig.sh ~/data/shared/hg19/hg19.ChromInfo.txt ${OUTPUT_PREFIX}Aligned.out.sam 1

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
echo "Mapping to 5' introns:" > ${log}

# Output the softclipped part on the left end of reads, aligned part, and softclipped part on the right end of reads (always on the reference strand)
segments_tab=${prefix}.before.align.after.table
~/repo/tools/parse_bam/get_before_alignment_after.py -s -f ${map_to_5_intron} 2>> ${log} > ${segments_tab}
cat ${segments_tab} | awk -v p=${prefix} 'BEGIN{bef=p".before.fa"; ali=p".align.fa"; aft=p".after.fa";} { if($2!="") {print ">" $1 "\n" $2> bef}; if($3!="") {print ">" $1 "\n" $3 > ali}; if($4!="") {print ">" $1 "\n" $4 > aft}  }'
before=${prefix}.before.fa 	# This is already generated by the previous line
align=${prefix}.align.fa 	# This is already generated by the previous line
align_after=${prefix}.align+after.fa
cat ${segments_tab} |awk '{ print ">" $1; print $3 $4 }' > ${align_after}
# Use those reads that map to the antisense strand of introns
# fpl: five prime length (23) requiring part of the read to match the 5' end of an intron
# tl: After trimming, there should be at least tl (20) nt left on the reads
# gawk -v fpl=23 -v tl=20 '{ pat="(^[0-9]+)S" fpl "M"; match($6, pat, ary); a=ary[1]; if(a>=tl) { print ">" $3 "|" substr($10, a, fpl); print substr($10, 1, a); }}' ${map_to_5_intron_c} >> ${map_to_5_intron_fa}
# Map the "before" part of the read to the introns

genome_index=/data/fuy2/shared/hg19/bowtie2index/hg19
map_to_genome_before=${prefix}.before.map_to_genome.bam
map_to_genome_before_log=${prefix}.before.map_to_genome.log
bowtie2 -x ${genome_index} -f -U ${before} -p ${cpu} 2> ${map_to_genome_before_log} | samtools view -bS -F 0x4 -q 10 - | samtools sort -@ ${CPU} -O bam -T ${RANDOM}${RANDOM} - > ${map_to_genome_before}
map_to_genome_before_bed=${prefix}.before.map_to_genome.bed
bedtools bamtobed -i ${map_to_genome_before} | sort -k4,4 >${map_to_genome_before_bed}

map_to_genome_align=${prefix}.align.map_to_genome.bam
map_to_genome_align_log=${prefix}.align.map_to_genome.log
bowtie2 -x ${genome_index} -f -U ${align} -p ${cpu} 2> ${map_to_genome_align_log} | samtools view -bS -F 0x4 -q 10 - | samtools sort -@ ${CPU} -O bam -T ${RANDOM}${RANDOM} - > ${map_to_genome_align}
map_to_genome_align_bed=${prefix}.align.map_to_genome.bed
bedtools bamtobed -i ${map_to_genome_align} | sort -k4,4 >${map_to_genome_align_bed}
map_to_genome_before_align_bed=${prefix}.before+after.map_to_genome.bed
join -1 4 -2 4 ${map_to_genome_before_bed} ${map_to_genome_align_bed} > ${map_to_genome_before_align_bed}

# Results after the first filter: strandedness and chromosome
f1=${prefix}.before+align.map_to_genome.f1.bed
cat ${map_to_genome_before_align_bed} | awk '{ if($6==$11 && $2==$7) { print } }' | awk 'BEGIN{ OFS="\t" } { c1=$2; s1=$3; e1=$4; st1=$6; c2=$7; s2=$8; e2=$9; st2=$11; s=s1<s2?s1:s2; e=e1<e2?e2:e1; print c1 "\t" s "\t" e "\t" ".\t0\t" st1 "\t" $0}' | sort -k1,1 -k2,2n > ${f1}

# 2nd filter: the before and align part has to be in the same intron
# It also requires the 'before' and 'align' segments should not be right next to each other
f2=${prefix}.before+align.map_to_genome.f2.bed
bedtools intersect -f 1.00 -a  ${f1} -b ~/repo/tools/bp/common/hg19/hg19.gencode.intron.uniq.bed -wa -u | awk '{ if($6=="+" && $10!=$14) {print} if($6=="-" &&$9!=$15) {print} }' | sort -k1,1 -k2,2n > ${f2}

f2_before=${prefix}.before+align.map_to_genome.f2.before.bed
f2_align=${prefix}.before+align.map_to_genome.f2.align.bed
cat ${f2} | awk -v l=100 '{st=$6; s=$9; e=$10; if(st=="+") { s=e-l } else { e=s+l }; print $1 "\t" s "\t" e "\t" $7 "\t" 0 "\t" st }' > ${f2_before}
cat ${f2} | awk -v l=100 '{st=$6; s=$14; e=$15; if(st=="+") { e=s+l } else { s=e-l }; print $1 "\t" s "\t" e "\t" $7 "\t" 0 "\t" st }' >  ${f2_align}

f2_before_tab=${prefix}.before+align.map_to_genome.f2.before.tab
f2_align_tab=${prefix}.before+align.map_to_genome.f2.align.tab
genome_fa=~/data/shared/hg19/hg19.fa
bedtools getfasta -s -tab -fi ${genome_fa} -bed ${f2_before} -fo - > ${f2_before_tab}
bedtools getfasta -s -tab -fi ${genome_fa} -bed ${f2_align}  -fo - > ${f2_align_tab}

# putative lariat index
put_lar_fa=${prefix}.put_lar.fa
put_lar_idx=${prefix}.put_lar
paste ${f2_before_tab} ${f2_align_tab} | awk 'BEGIN{i=1} { print ">" $1 "|" $3 "|" i++; print $2 $4; }' > ${put_lar_fa}
bowtie2-build -o 1 ${put_lar_fa} ${put_lar_idx}

map_to_put_lar=${prefix}.map_to_put_lar.bam
bowtie2 --local -x ${put_lar_idx} -U ${FQ} | samtools view -b -F 0x4 - | samtools sort -@ ${CPU} -T ${RANDOM}${RANDOM} -O bam - > ${map_to_put_lar}


# By testing, the 3' ends of reads are likely to be junk and I cannot align the full
# length align+after parts, and thus these are deprecated.
# map_to_genome_align_after=${prefix}.align+after.map_to_genome.bam
# map_to_genome_align_after_log=${prefix}.align+after.map_to_genome.log
# bowtie2 -x ${genome_index} -f -U ${align_after} -p ${cpu} 2> ${map_to_genome_align_after_log} | samtools view -bS -F 0x4 -q 10 - | samtools sort -@ ${CPU} -O bam -T - > ${map_to_genome_align_after}
# map_to_genome_align_after_bed=${prefix}.align+after.map_to_genome.bed
# bedtools bamtobed -i ${map_to_genome_align_after} | sort -k4,4 >${map_to_genome_align_after_bed}



# Get the weblogo
# p=${prefix}.lariat_and_5intron.map_to_intron
# bedtools bamtobed -i ${p}.bam | bedtools getfasta -s -fi ~/repo/tools/bp/common/hg19/hg19.gencode.intron.uniq.fa -bed - -fo - | awk '{ if(NR%2==1) {print $0} else { print substr($0, length($0)-10+1) } }' | weblogo -D fasta -F pdf > ${p}.last10.logo.pdf

# Extract the concordant reads and get bigWig files
# concord=${prefix}.lariat_and_5intron.map_to_intron.concord.fastq
# samtools view -b -F0x10 ${map_to_intron} | bedtools bamtofastq -i - -fq /dev/stdout > ${concord}
# genome_index=~/data/shared/hg19/bowtie2index/hg19
# concord_to_genome=${prefix}.lariat_and_5intron.map_to_intron.concord.map_to_genome.bam
# concord_to_genome_log=${prefix}.lariat_and_5intron.map_to_intron.concord.map_to_genome.log
# bowtie2 -p ${CPU} -x ${genome_index} -U ${concord} 2> ${concord_to_genome_log} | samtools view -b -F 0x4 - > ${concord_to_genome} 
# bam_to_bigwig.sh ~/data/shared/hg19/hg19.ChromInfo.txt ${concord_to_genome} 1

# Make lariat index of (100+100)nt: upstream of branchpoint and 5'ss

# TODO
# samtools view -F0x10 Mattick.RNaseR.HeLa.BP.rep1.1000000.fq.gz.lariat_and_5intron.map_to_intron.bam | awk '{ split($1, a, "|"); match(a[1], /(.+):[0-9]+-[0-9]+\([+-]\)/, b); q_chr=b[1]; match($3, /(.+):[0-9]+-[0-9]+\([+-]\)/, c); ref_chr=c[1]; if(q_chr==ref_chr) { print } }' | less
# STATS
# Number of reads (before and align are both mappable)
# cat Mattick.RNaseR.HeLa.BP.rep1.fq.gz.before+align.map_to_genome.bed | wc -l
# Number of concordant reads (before and align map to the strand strand on the same chromosome)
# cat Mattick.RNaseR.HeLa.BP.rep1.fq.gz.before+align.map_to_genome.bed | awk '{ if($6==$11 && $2==$7) { print } }' | wc -l
# wc -l ${f2}
