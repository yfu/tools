#!/usr/bin/env bash

# Count the number of reads falling into each type of feature. Special attention should be paid to those reads that can map to multiple features
# e.g. one read can map to two loci
# chr1	1	25	1	2	ATCGATCGATCGATCGATACGT
# When intersect this bed2 file with the following bed file for 5UTR and 3UTR
# chr1	1	100	5UTR
# chr1	1	120	3UTR

bed2=$1

base=/data/fuy2/piPipes/common/mm9

piRNA_cluster=${base}/mm9.piRNAcluster.bed.gz
Genes=${base}/UCSC.refSeq.Genes.bed12.gz
RepeatMasker=${base}/UCSC.RepeatMask.bed

piRNA_cluster_prepachytene=${base}/piRNA.cluster.prepachytene.bed12.gz
piRNA_cluster_hybrid=${base}/piRNA.cluster.hybrid.bed12.gz
piRNA_cluster_pachytene=${base}/piRNA.cluster.pachytene.bed12.gz

piRNA_cluster_prepachytene_exons=${base}/piRNA.cluster.prepachytene.bed6.gz
piRNA_cluster_hybrid_exons=${base}/piRNA.cluster.hybrid.bed6.gz
piRNA_cluster_pachytene_exons=${base}/piRNA.cluster.pachytene.bed6.gz

Genes_Exons=${base}/UCSC.refSeq.Exon.bed.gz
Genes_Exons_5UTR=${base}/UCSC.refSeq.5UTR.bed.gz
Genes_Exons_3UTR=${base}/UCSC.refSeq.3UTR.bed.gz
Genes_Intron=${base}/UCSC.refSeq.Intron.bed.gz

RepeatMasker_DNA=${base}/UCSC.RepeatMask.bed.DNA.gz
RepeatMasker_LINE=${base}/UCSC.RepeatMask.bed.LINE.gz
RepeatMasker_LTR=${base}/UCSC.RepeatMask.bed.LTR.gz
RepeatMasker_Satellite=${base}/UCSC.RepeatMask.bed.Satellite.gz
RepeatMasker_SimpleRepeat=${base}/UCSC.RepeatMask.bed.Simple_repeat.gz
RepeatMasker_SINE=${base}/UCSC.RepeatMask.bed.SINE.gz
RepeatMasker_rtRNA=${base}/UCSC.RepeatMask.bed.rtRNA.gz

bn=$(basename $bed2)
prefix=${bn%.bed2}
sort -k1,1 -k2,2n $bed2 > ${prefix}.sorted.bed2
bed2=${bn%.bed}.sorted.bed2

count() {
    targets=($@)
    exclusive_targets=${bn}.exclusive_targets
    
    if [ -f $exclusive_targets ]; then
	rm $exclusive_targets
    fi
    touch $exclusive_targets
    
    for t in "${targets[@]}"; do
	target=${!t}
	if [[ ${target##*.} == gz ]] || [[ ${target##*.gzip} == gzip ]]; then
	    zcat $target | awk -v t=$t '{ print $1 "\t" $2 "\t" $3 "\t" t }' >> $exclusive_targets
	else
	    cat $target | awk -v t=$t '{ print $1 "\t" $2 "\t" $3 "\t" t }' >> $exclusive_targets
	fi
    done
    
    sort -k1,1 -k2,2n $exclusive_targets > ${exclusive_targets}.tmp && mv ${exclusive_targets}.tmp ${exclusive_targets}
    
    bed2_intersect=$bn.vs.$exclusive_targets
    
    bedtools intersect -sorted -wa -wb -a $bed2 -b $exclusive_targets > $bed2_intersect
    
    ## Notice that one read might map to multiple isoforms...
    awk '{ if(NR==FNR) { ntm[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7] +=1 } else { print $0 "\t" ntm[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7] } }' $bed2_intersect $bed2_intersect | awk '{ cate[$(NF-1)] += $4/$5/$NF } END{ for (c in cate) { print c "\t" cate[c] } }'
    ## Also output read counts for the unannotated
    bedtools intersect -v -a $bed2 -b $exclusive_targets | awk '{ a += $4/$5 } END{ print "Other\t" a }'
}


targets=(piRNA_cluster Genes RepeatMasker)
count "${targets[@]}" > ${prefix}.piRNA_cluster.Genes.RepeatMasker.count

targets=(piRNA_cluster_prepachytene piRNA_cluster_hybrid piRNA_cluster_pachytene)
count "${targets[@]}" > ${prefix}.piRNA_cluster_class.count

targets=(piRNA_cluster_prepachytene_exons piRNA_cluster_hybrid_exons piRNA_cluster_pachytene_exons)
count "${targets[@]}" > ${prefix}.piRNA_cluster_exon_level.count

targets=(Genes_Exons Genes_Intron)
count "${targets[@]}" > ${prefix}.Genes_level.count

targets=(Genes_Exons_5UTR Genes_Exons_3UTR Genes_Intron)
count "${targets[@]}" > ${prefix}.Genes_utr_level.count

targets=(RepeatMasker_DNA RepeatMasker_LINE RepeatMasker_LTR RepeatMasker_Satellite RepeatMasker_SimpleRepeat RepeatMasker_SINE RepeatMasker_rtRNA)
count "${targets[@]}" > ${prefix}.RepeatMasker_leve.count
