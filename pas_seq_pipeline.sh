#!/usr/bin/env bash

# Start from sorted bam. Will output two sets of bigWig files (for unique mappers and all mappers, respectively).
# Get the norm_factor this way: cat /data/fuy2/prepachytene/results/2014-11-28/polya_removed_mm10g/Zamore.PAS.wildtype.testis.${i}.rep1.trimmed.polya_removed.fq.mm10g.Log.final.out | grep "Uniquely mapped reads number" - | grep -E -o '[0-9]+'
chrom_info=~/data/piPipes/common/mm10g/mm10g.ChromInfo.txt
cpu=8
# norm_factor=0.001

# bam=Zamore.PAS.wildtype.testis.07dpp.rep1.trimmed.polya_removed.fq.mm10g.Aligned.out.sorted.bam
bam=$1
# Zamore.PAS.wildtype.testis.07dpp.rep1.trimmed.polya_removed.fq.mm10g.Log.final.out
log=${bam%.Aligned.out.sorted.bam}.Log.final.out
uniq_reads=$(cat ${log} | grep "Uniquely mapped reads number" - | grep -E -o '[0-9]+')
echo "# reads mapping to unique loci: $uniq_reads" 
norm_factor_uniq=$(awk -v uniq_reads=$uniq_reads 'BEGIN{ print 1e6/uniq_reads }')
echo "norm_factor_uniq: $norm_factor_uniq"

# samtools view -h ${bam} | head -n 100000 | samtools view -Sb - > test.sorted.bam
# bam=test.sorted.bam

prefix=${bam%.sorted.bam}
# bed2_3prime only contains one nucleotide at the 3' end
bed2=${prefix}.bed2
bed2_3prime=${prefix}.3prime.bed2
samtools view -h ${bam} | awk '{ if($0 ~/^@/) {print} else { cigar=$6; if(cigar ~ /S$/) {} else{print} } }' | samtools view -Sb - | bedtools bamtobed -bed12 -i - -tag NH | awk 'BEGIN{ FS=OFS="\t" } { $4=1; print $0 }' | sort -k1,1 -k2,2n --parallel=8 | cut -f1-7 > ${bed2}
awk 'BEGIN{FS=OFS="\t"} { if($6=="+") {$2=$3-1} else { $3=$2+1 }; print $0 }' ${bed2} | sort -k1,1 -k2,2n --parallel=$cpu > ${bed2_3prime}

# For unique mappers
watson_uniq_bedgraph=${prefix}.Watson.uniq.bedGraph
watson_uniq_bigwig=${prefix}.Watson.uniq.bigWig
crick_uniq_bedgraph=${prefix}.Crick.uniq.bedGraph
crick_uniq_bigwig=${prefix}.Crick.uniq.bigWig
awk '{ if($5==1) {print} }' ${bed2_3prime} | bedtools genomecov -i test.3prime.bed2 -g ~/data/piPipes/common/mm10g/mm10g.ChromInfo.txt -bg -scale $norm_factor_uniq -strand + > $watson_uniq_bedgraph
awk '{ if($5==1) {print} }' ${bed2_3prime} | bedtools genomecov -i test.3prime.bed2 -g ~/data/piPipes/common/mm10g/mm10g.ChromInfo.txt -bg -scale $norm_factor_uniq -strand - | awk 'BEGIN{FS=OFS="\t"} { $4=-$4; print $0 }' > $crick_uniq_bedgraph

bedGraphToBigWig ${watson_uniq_bedgraph} ${chrom_info} ${watson_uniq_bigwig}
bedGraphToBigWig ${crick_uniq_bedgraph} ${chrom_info} ${crick_uniq_bigwig}

# For all mappers
watson_all_bedgraph=${prefix}.Watson.all.bedGraph
watson_all_bigwig=${prefix}.Watson.all.bigWig
crick_all_bedgraph=${prefix}.Crick.all.bedGraph
crick_all_bigwig=${prefix}.Crick.all.bigWig

# Number of reads mapped to multiple loci |       10122442
multi_reads=$(cat ${log} | grep "Number of reads mapped to multiple loci" - | grep -E -o '[0-9]+')
echo "# reads mapping to multiple loci: $multi_reads" 
norm_factor_all=$(awk -v multi_reads=$multi_reads -v uniq_reads=$uniq_reads 'BEGIN{ print 1e6/(multi_reads + uniq_reads) }')
echo "norm_factor_all: $norm_factor_all"

watson_all_bedgraph_tmp=${watson_all_bedgraph}.tmp
crick_all_bedgraph_tmp=${crick_all_bedgraph}.tmp
bed2_to_bedgraph.py ${bed2_3prime} $watson_all_bedgraph_tmp $crick_all_bedgraph_tmp
awk -v norm_factor_all=$norm_factor_all 'BEGIN{ FS=OFS="\t" } { $4=$4*norm_factor_all; print }' $watson_all_bedgraph_tmp > $watson_all_bedgraph
awk -v norm_factor_all=$norm_factor_all 'BEGIN{ FS=OFS="\t" } { $4=$4*norm_factor_all; print }' $crick_all_bedgraph_tmp > $crick_all_bedgraph

bedGraphToBigWig ${watson_all_bedgraph} ${chrom_info} ${watson_all_bigwig}
bedGraphToBigWig ${crick_all_bedgraph} ${chrom_info} ${crick_all_bigwig}

