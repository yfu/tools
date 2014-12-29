
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

Genes_Upstream10k=UCSC.refSeq.Genes.upstream10k.bed12
Genes_Downstream10k=UCSC.refSeq.Genes.downstream10k.bed12

## targets=(Genes piRNA_cluster piRNA_cluster_prepachytene piRNA_cluster_hybrid piRNA_cluster_pachytene Genes_Exons Genes_Exons_5UTR Genes_Exons_3UTR Genes_Intron)

targets=(Genes_Upstream10k Genes_Downstream10k)

## targets=(Genes_Exons Genes_Exons_5UTR Genes_Exons_3UTR Genes_Intron)
## Zamore.SRA.total.00dpp.testis.trimmed.x_rRNA.x_hairpin.mm9v1.unique.bed2
for i in 00 02 04 07 10 12 14 17 20 42; do
## for i in 02; do
    for j in "${targets[@]}"; do
	t=${!j}
	output=${i}dpp.${j}.aggregation.piRNA.100bin.signal
	log=${i}dpp.${j}.aggregation.piRNA.100bin.signal.log
	bed2=~/data/prepachytene/results/2014-11-21/SRA/${i}dpp/genome_mapping/Zamore.SRA.total.${i}dpp.testis.trimmed.x_rRNA.x_hairpin.mm9v1.all.piRNA.bed2
	## echo "bash -x aggregate_bed2_into_bins.sh ${t} ~/data/prepachytene/results/2014-11-21/SRA/${i}dpp/genome_mapping/Zamore.SRA.total.${i}dpp.testis.trimmed.x_rRNA.x_hairpin.mm9v1.unique.bed2 > ${output}"
	echo "aggregate_bed2_into_bins.py ${t} ${bed2} >${output} 2>${log}"
	output=${i}dpp.${j}.aggregation.siRNA.100bin.signal
	log=${i}dpp.${j}.aggregation.siRNA.100bin.signal.log
	bed2=~/data/prepachytene/results/2014-11-21/SRA/${i}dpp/genome_mapping/Zamore.SRA.total.${i}dpp.testis.trimmed.x_rRNA.x_hairpin.mm9v1.all.siRNA.bed2	
	echo "aggregate_bed2_into_bins.py ${t} ${bed2} >${output} 2>${log}"
    done
done


