# This script intersect a bed file with multiple features, regardless of strands

cp ~/data/shared/mm9ToMm10/piRNA.cluster.bed12 ./
cp ~/data/shared/mm9ToMm10/piRNA.cluster.bed6 ./

for i in prepachytene hybrid pachytene; do
    cp ~/data/shared/mm9ToMm10/piRNA.cluster.${i}.bed12 ./
    cp ~/data/shared/mm9ToMm10/piRNA.cluster.${i}.bed6 ./
done

rmsk_dir=/home/fuy2/data/shared/mm10
rmsk=${rmsk_dir}/UCSC.rmsk.sorted.bed

rmsk_dna=${rmsk_dir}/UCSC.rmsk.DNA.sorted.bed
rmsk_line=${rmsk_dir}/UCSC.rmsk.LINE.sorted.bed
rmsk_ltr=${rmsk_dir}/UCSC.rmsk.LTR.sorted.bed
rmsk_rrna=${rmsk_dir}/UCSC.rmsk.rRNA.sorted.bed
rmsk_satellite=${rmsk_dir}/UCSC.rmsk.Satellite.sorted.bed
rmsk_simplerepeats=${rmsk_dir}/UCSC.rmsk.Simple_repeat.sorted.bed
rmsk_sine=${rmsk_dir}/UCSC.rmsk.SINE.sorted.bed
rmsk_trna=${rmsk_dir}/UCSC.rmsk.tRNA.sorted.bed

cluster=piRNA.cluster.bed12
cluster_exon=piRNA.cluster.bed6
prepachytene=piRNA.cluster.prepachytene.bed12
prepachytene_exon=piRNA.cluster.prepachytene.bed6
hybrid=piRNA.cluster.hybrid.bed12
hybrid_exon=piRNA.cluster.hybrid.bed6
pachytene=piRNA.cluster.pachytene.bed12
pachytene_exon=piRNA.cluster.pachytene.bed6


# 1st step: intersect with rmsk
## dpp=$1
## mkdir -p ${dpp}
## Real data: ~/data/prepachytene/results/2014-11-21/SRA_mm10g/${dpp}/genome_mapping/Zamore.SRA.total.${dpp}.testis.trimmed.x_rRNA.x_hairpin.mm10gv1.all.piRNA.bed2
orig=$1
prefix=$2
input=${prefix}.sorted.bed2
sort --parallel=4 -k1,1 -k2,2n ${orig} > ${input}

rmsk_bed2=${prefix}.rmsk.bed2
x_rmsk_bed2=${prefix}.x_rmsk.bed2
bedtools intersect -sorted -a ${input} -b ${rmsk} -u > ${rmsk_bed2}
bedtools intersect -sorted -a ${input} -b ${rmsk} -v > ${x_rmsk_bed2}

for i in dna line ltr rrna satellite simplerepeats sine trna; do
    # different repeats usually do not overlap, simply do the following
    eval "rep=\${rmsk_${i}}"
    bedtools intersect -a ${rmsk_bed2} -b $rep -u -sorted > ${prefix}.rmsk_${i}.bed2
done

for i in protein_coding lincRNA pseudogene; do
    # GTF and bed
    awk -v i=${i} '{ if($3=="gene" && $14=="\""i "\";") {print} }' ~/data/shared/mm10/gencode.vM4.corrected_chrom_names.gtf > gene.${i}.gtf
    cat gene.${i}.gtf | tr -d '\";' | awk 'BEGIN{OFS="\t"} { print $1, $4-1, $5, $10, $18, $7 }' > gene.${i}.bed
    # Exon bed
    awk -v i=${i} '{ if($14=="\"" i "\";") {print} }' ~/data/shared/mm10/gencode.vM4.corrected_chrom_names.gtf | gtfToGenePred -ignoreGroupsWithoutExons stdin stdout | genePredToBed stdin stdout | sort -k1,1 -k2,2n > transcript.${i}.bed
    bedtools bed12tobed6 -i transcript.${i}.bed | sort -k1,1 -k2,2n > exon.${i}.bed
    # upstream
    bedtools flank -s -l 5000 -r 0 -i gene.${i}.bed -g /home/fuy2/data/piPipes/common/mm10g/mm10g.ChromInfo.txt | sort -k1,1 -k2,2n > gene.${i}.upstream5k.bed
    bedtools flank -s -r 5000 -l 0 -i gene.${i}.bed -g /home/fuy2/data/piPipes/common/mm10g/mm10g.ChromInfo.txt | sort -k1,1 -k2,2n > gene.${i}.downstream5k.bed

    # Set up the variables: gene_{protein_coding,lincRNA,pseudogene} and gene_{protein_coding,lincRNA,pseudogene}_exon
    eval ${i}=gene.${i}.bed
    eval ${i}_exon=exon.${i}.bed
    eval ${i}_upstream5k=gene.${i}.upstream5k.bed
    eval ${i}_downstream5k=gene.${i}.downstream5k.bed    
done

echo "Test: ${protein_coding_exon}" >&2

# Assign the reads to different categories one by one, because 
targets=(cluster lincRNA pseudogene protein_coding protein_coding_upstream5k protein_coding_downstream5k lincRNA_upstream5k lincRNA_downstream5k pseudogene_upstream5k pseudogene_downstream5k)
cur_bed=${x_rmsk_bed2}
for i in ${targets[@]}; do
    # Different categories overlap, so only keep those reads (alignments) that cannot be assigned to the current category for the next iteration
    eval t=\$$i
    echo "Processing target: $t" >&2
    # echo "bedtools intersect -a ${x_rmsk_bed2} -b $t -u -s -sorted > ${prefix}.${i}.s.bed2"
    bedtools intersect -a ${cur_bed} -b $t -u -sorted > ${prefix}.${i}.bed2
    # echo "bedtools intersect -a ${x_rmsk_bed2} -b $t -u -S -sorted > ${prefix}.${i}.as.bed2"
    bedtools intersect -a ${cur_bed} -b $t -v -sorted > ${cur_bed%.bed2}.x_${i}.bed2
    cur_bed=${cur_bed%.bed2}.x_${i}.bed2
    echo "Now using $cur_bed"
done

echo "What are left: ${cur_bed}" >&2
mv ${cur_bed} ${prefix}.other.bed2
