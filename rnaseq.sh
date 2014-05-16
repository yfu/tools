#!/usr/bin/env bash

# Usage: rnaseq.sh -l 10131619.r1_id1.fq -r 10131619.r2_id1.fq 

cpu=4

while getopts "hfl:r:o:s:c:FR" OPTION
do
     case $OPTION in
          l)
               left=`readlink -f $OPTARG`
          ;;
          r)
               right=`readlink -f $OPTARG`
          ;;
          o)
               out_dir=$OPTARG
          ;;
          # F)
          #      strand=$((STRAND+2))
          # ;;
          # R)
          #      strand=$((STRAND+1))
          # ;;
     esac
done


step=1
common="$HOME/data/PE_Pipeline/common_files/dm3"
phred_option="--phred33"
out_dir=$(pwd)
strand=1

read_length=`head -2 $left | awk '{getline; printf "%d", length($1)}'`
# overhang
star_sjdboverhang=$((read_length-1))

echo "########################################################################"
out_dir1=${out_dir}/x_rRNA
mkdir -p $out_dir1
[ -f ${out_dir}/.${step}.x_rRNA ] && echo "x_rRNA is done already! I will skip this step"
[ ! -f ${out_dir}/.${step}.x_rRNA ] &&
bowtie2 -x ${common}/rRNAIndex/rRNA -1 $left -2 $right -q $phred_option --very-fast -k 1 \
    --no-mixed --no-discordant --un-conc ${out_dir1}/x_rRNA.fq \
    -p $cpu -S /dev/null 2> ${out_dir1}/x_rRNA.log && \
touch ${out_dir}/.${step}.x_rRNA && \
echo "x_rRNA is done."
step=$((step+1))

echo "########################################################################"
# Genome mapping

star_genome_index=${common}/STARgenomeIndex/${star_sjdboverhang}

out_dir2=${out_dir}/genomeMapping
mkdir -p ${out_dir2}
cd ${out_dir2} && \
left=${out_dir1}/x_rRNA.1.fq && \
right=${out_dir1}/x_rRNA.2.fq
[ -f ${out_dir}/.${step}.genome_mapping ] && echo "Genome mapping is done already! I will skip this step."
[ ! -f ${out_dir}/.${step}.genome_mapping ] && \
STAR \
     --runMode alignReads \
     --genomeDir $star_genome_index \
     --readFilesIn ${left} ${right} \
     --runThreadN $cpu \
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
     --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
     --genomeLoad NoSharedMemory \
     --outFileNamePrefix x_rRNA. \
     --outSAMunmapped None \
     --outReadsUnmapped Fastx \
     --outSJfilterReads Unique \
     --seedSearchStartLmax 20 \
     --seedSearchStartLmaxOverLread 1.0 \
     --chimSegmentMin 0 && \
    touch ${out_dir}/.${step}.genome_mapping && \
    echo "Genome mapping is done."

step=$((step+1))

echo "########################################################################"
echo "Getting Statistics..."

# getting statistics
input_reads=`grep 'Number of input reads' x_rRNA.Log.final.out | awk '{print $NF}'`
uniq_reads=`grep 'Uniquely mapped reads number' x_rRNA.Log.final.out | awk '{print $NF}'`
multi_reads=`grep 'Number of reads mapped to multiple loci' x_rRNA.Log.final.out | awk '{print $NF}'`
all_mapped_reads=$((uniq_reads+multi_reads))
unmapped_reads=$((input_reads-uniq_reads-multi_reads))

cat <<EOF
Input reads: $input_reads
Unique reads: $uniq_reads
Multi reads: $multi_reads
All mapped reads: $all_mapped_reads
Unmapped reads: $unmapped_reads
EOF

echo "########################################################################"
# Genome bam processing
cd $out_dir

[ -f ${out_dir}/.${step}.genome_bam_processing ] && echo "Genome bam processing is done already. I will skip this step."

# echo ${out_dir2}/x_rRNA.Aligned.out.sam

[ ! -f ${out_dir}/.${step}.genome_bam_processing ] && \
samtools view -bS ${out_dir2}/x_rRNA.Aligned.out.sam | tee ${out_dir2}/x_rRNA.Aligned.out.bam | bedtools bamtobed -bedpe -mate1 -i - > ${out_dir2}/x_rRNA.bedpe && \

samtools sort -@ $cpu ${out_dir2}/x_rRNA.Aligned.out.bam ${out_dir2}/x_rRNA.sorted && \
    samtools index ${out_dir2}/x_rRNA.sorted.bam && \
    awk 'BEGIN{FS=OFS="\t"}{if (ARGIND==1) {++ct[$7]} else {$8=1.0/ct[$7]; print $0 > "/dev/stdout"; if (ct[$7]==1) print $0 > "/dev/stderr"}}' ${out_dir2}/x_rRNA.bedpe ${out_dir2}/x_rRNA.bedpe \
    1> ${out_dir2}/x_rRNA.all.bedpe \
    2> ${out_dir2}/x_rRNA.unique.bedpe && \
    rm -rf ${out_dir2}/x_rRNA.bedpe && \
    touch ${out_dir}/.${step}.genome_bam_processing
step=$((step+1))

echo "########################################################################"

# Prepare for cufflinks
# declare -a cufflinks_g=("iGenome")
declare -a cufflinks_g=("flyBase" "iGenome")
iGenome=/data/fuy2/shared/iGenome.gtf
# /home/hanb/nearline/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf
flyBase=/data/fuy2/PE_Pipeline/common_files/dm3/dmel-all-no-analysis-r5.54.gene+CDS.gtf
# /home/hanb/nearline/PE_Pipeline/common_files/dm3/dmel-all-no-analysis-r5.54.gene+CDS.gtf
library_type="fr-firststrand"
genome_fa=/data/fuy2/shared/dm3.fa

# Cufflinks
for t in "${cufflinks_g[@]}"; do \
     out_dir_sub=${out_dir}/cufflinks_${t}_compatible_hits_norm && mkdir -p $out_dir_sub
     [ ! -f ${out_dir}/.${step}.quantification_by_cuff_${t}_compatible_hits_norm ] && \
          cufflinks -v \
          -o $out_dir_sub \
          -p $cpu \
          -G ${!t} \
          -b $genome_fa \
          -u \
          --library-type $library_type \
          --compatible-hits-norm \
          --no-update-check \
          ${out_dir2}/x_rRNA.sorted.bam \
          2> $out_dir_sub/cufflinks.log && \
          touch ${out_dir}/.${step}.quantification_by_cuff_${t}_compatible_hits_norm
     step=$((step+1))

     out_dir_sub=${out_dir}/cufflinks_${t}_total_hits_norm
     mkdir -p $out_dir_sub
     [ ! -f ${out_dir}/.${step}.quantification_by_cuff_${t}_total_hits_norm ] && \
          cufflinks -v \
          -o $out_dir_sub \
          -p $cpu \
          -G ${!t} \
          -b $genome_fa \
          -u \
          --library-type $library_type \
          --total-hits-norm \
          --no-update-check \
          ${out_dir2}/x_rRNA.sorted.bam \
          2> $out_dir_sub/cufflinks.log && \
          touch ${out_dir}/.${step}.quantification_by_cuff_${t}_total_hits_norm
     step=$((step+1))
done

echo "########################################################################"


#use htseq-count
out_dir3=${out_dir}/htseqCount
mkdir -p ${out_dir3}
[ -f ${out_dir}/.${step}.overlap_by_htseq-count ] && echo "htseq-count is done already. I will skip this step."
[ ! -f ${out_dir}/.${step}.overlap_by_htseq-count ] && \
htseqCountPE.sh \
     ${out_dir2}/x_rRNA.Aligned.out.sam \
     $strand \
     fly \
     $out_dir3 \
     $cpu && \
     touch ${out_dir}/.${step}.overlap_by_htseq-count
step=$((step+1))
echo "########################################################################"
