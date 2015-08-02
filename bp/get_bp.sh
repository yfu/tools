#!/usr/bin/env bash

binpath=$HOME/repo/tools/bp
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
statusf=.${PREFIX}.regular_mapping.ok
if [[ -f ${statusf} ]]; then
    echo "Regular RNA-seq mapping has been done. If you want to redo it, please delete ${statusf}"
else
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
	--chimSegmentMin 0 2>&1 1> ${OUTPUT_PREFIX}STAR.log && touch ${statusf}
fi
# Make bigWig files from the regular RNA-seq reads
# Normalization factor = 1
bam_to_bigwig.sh ~/data/shared/hg19/hg19.ChromInfo.txt ${OUTPUT_PREFIX}Aligned.out.sam 1

###########################################################################
# Part 2: mapping reads to 5' introns to extract useful part of the reads #
###########################################################################

input=${PREFIX}.${GENOME}.Unmapped.out.mate1
prefix=${input%.hg19.Unmapped.out.mate1}
cpu=${CPU}

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
statusf=.${prefix}.5intron_mapping.ok
if [[ -f ${statusf} ]]; then
    echo "5' intron mapping has been done. If you want to redo it, please delete ${statusf}"
else
    bowtie2 --local --score-min L,45,0 -D 20 -R 2 -N 0 -L 20 -i L,1,0 --phred33 -p ${cpu} -x ${index_5pintron} -U ${input} 2> ${log} | samtools view -q 10 -bS -F 0x4 - | samtools sort -@ ${cpu} -T ${RANDOM} -O bam > ${map_to_5_intron} && touch ${statusf}
fi
    
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
bedtools getfasta -s -tab -fi ${genome_fa} -bed ${f2_before} -fo - | awk '{ print substr($2, length($2)) "|" $1 "\t" $2 }' > ${f2_before_tab}
bedtools getfasta -s -tab -fi ${genome_fa} -bed ${f2_align}  -fo - > ${f2_align_tab}

# putative lariat index
put_lar_fa=${prefix}.put_lar.fa
put_lar_pdf=${prefix}.put_lar.pdf
put_lar_idx=${prefix}.put_lar
put_lar_idx_log=${prefix}.put_lar.idx.log
# Only keep the uniq lariat sequence
paste ${f2_before_tab} ${f2_align_tab} | sort -k1,1 | uniq | awk 'BEGIN{i=1} { print ">" $1 "|" $3 "|" i++; print $2 $4; }' > ${put_lar_fa}
samtools faidx ${put_lar_fa}
cat ${put_lar_fa} | weblogo -F pdf > ${put_lar_pdf}
bowtie2-build -o 1 ${put_lar_fa} ${put_lar_idx} &> ${put_lar_idx_log}

map_to_put_lar=${prefix}.map_to_put_lar.bam
map_to_put_lar_log=${prefix}.map_to_put_lar.log
# Minimum of 20 nt match
# Note that multi-mappers should be kept because lariats may have overlaps
bowtie2 -p ${CPU} --local --score-min L,40,0 -D 20 -R 2 -N 0 -L 20 -i L,1,0 -x ${put_lar_idx} -U ${input} 2> ${map_to_put_lar_log} | samtools view -b -F 0x4 - | samtools sort -@ ${CPU} -T ${RANDOM}${RANDOM} -O bam - > ${map_to_put_lar}
samtools index ${map_to_put_lar}

# Require alignment to cover the up-/down-stream of branchpoint
# baa: before, align and after
# map_to_put_lar_baa=${prefix}.map_to_put_lar.baa.tab
# map_to_put_lar_baa_log=${prefix}.map_to_put_lar.baa.log
# python ~/repo/tools/parse_bam/get_before_alignment_after.py -s -f ${map_to_put_lar} > ${map_to_put_lar_baa} 2>${map_to_put_lar_baa_log}

# Visualize the alignment by
# samtools view Mattick.RNaseR.HeLa.BP.rep1.fq.gz.map_to_put_lar.csrbp.bam | cut -f3 | awk '{ d[$1]+=1 } END{ for(i in d) { print i "\t" d[i] } }' | sort -k2,2nr | head
# samtools tview -p 'chr12:56213921-56214021(+)|chr12:56213277-56213377(+)|1733' Mattick.RNaseR.HeLa.BP.rep1.fq.gz.map_to_put_lar.csrbp.bam Mattick.RNaseR.HeLa.BP.rep1.fq.gz.put_lar.fa
crsbp=${prefix}.map_to_put_lar.crsbp.bam
python ~/repo/tools/parse_bam/filter_ol_bp.py -a 100 -l 5 -r 5 -f ${map_to_put_lar} > ${crsbp}
samtools index ${crsbp}
crsbp_list=${prefix}.map_to_put_lar.crsbp.list
samtools view ${crsbp} | awk '{ d[$3]+=1 } END{ for(i in d) { print i "\t" d[i] } }' | sort -k2,2nr > ${crsbp_list}

# crsbp_bcf=${prefix}.map_to_put_lar.crsbp.bcf
crsbp_bcf_log=${prefix}.map_to_put_lar.crsbp.pileup.log
crsbp_vcf=${prefix}.map_to_put_lar.crsbp.vcf
samtools mpileup -d10000000 -u -g -f ${put_lar_fa} ${crsbp} 2> ${crsbp_bcf_log} | bcftools call -A -m - > ${crsbp_vcf}
freq_prefix=${prefix}.map_to_put_lar.vcf
freq_log=${freq_prefix}.vcftools.log
vcftools --vcf ${crsbp_vcf} --counts --out ${freq_prefix} &> ${freq_log}
freq=${freq_prefix}.frq

sub=${prefix}.map_to_put_lar.sub
python ~/repo/tools/parse_bam/bp_err.py -f ${map_to_put_lar} -r ${put_lar_fa} > ${sub}

# NOTICE that only lariat supporting reads are considered
err_rate=${prefix}.map_to_put_lar.error_rate
grep -v '^#' ${sub} | awk '{ pos=$2; ref=toupper($3); d["A"]=$4; d["C"]=$5; d["G"]=$6; d["T"]=$7; d["N"]=$8; total[pos]+=$4+$5+$6+$7; t=$4+$5+$6+$7; correct=d[ref]; err=t-correct; pos_err[pos]+=err} END{ for(i in pos_err) { print i "\t" total[i] "\t" pos_err[i] } }' > ${err_rate}

final_lar_support=${prefix}.final.lar.support
final_lar=${prefix}.final.lar.bed
final_lar_bp=${prefix}.final.lar.bp.bed
samtools view ${map_to_put_lar} | awk '{ d[$3]+=1 } END{ for(i in d) { print i "\t" d[i] } }' | sort -k2,2nr > ${final_lar_support}
cat ${final_lar_support} | awk '{ split($0, a, "|"); match(a[2], /(.+):([0-9]+)-([0-9]+)\(([+-])\)/, m); match(a[3], /(.+):([0-9]+)-([0-9]+)\(([+-])\)/, n); print m[1] "\t"  m[2] "\t" m[3] "\t" ++i "\t" $2 "\t" m[4] "\t" n[1] "\t" n[2] "\t" n[3] "\t" n[4] }' > ${final_lar}
cat ${final_lar} | awk 'BEGIN{OFS="\t"} { if($6=="+") { $2=$3-1; } else { $3=$2+1 }; print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 }' > ${final_lar_bp}

anno_mattick=~/data/bp/shared/Mattick.bp.hg19.bed
inters=${prefix}.final.lar.bp.closest.bed
bedtools closest -t first -s -a ${final_lar_bp} -b ${anno_mattick} | sort -k11nr > ${inters}
stats=${prefix}.final.lar.stat
${binpath}/get_bp_stats.sh ${prefix} > ${stats}

# bcftools call -v -m Mattick.RNaseR.HeLa.BP.rep1.fq.gz.map_to_put_lar.crsbp.raw.bcf | grep -v '^##' | head

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

# cat Mattick.RNaseR.HeLa.BP.rep1.fq.gz.map_to_put_lar.before.align.after.tab | awk '{ if(length($3)>40 && length($2)<80 && length($4)<80) {print $2}}' | awk '{ match($0, /.+(AG|GT)$/, a); if(length(a[1])==0) { print } }' | sort^C uniq | awk '{ if(length($0)>20) { print substr($0, length($0)-10, 10) }}' | awk '{print ">" i++; print $0}' | weblogo -F pdf > test.pdf
