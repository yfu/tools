
funky() {
    # Given the category for the bed2 file and the bed2 file, this function prints the read count
    category=$1
    bed2_category=$2
    echo -ne "$category\t"
    cat ${bed2_category} | awk '{ copy[$7]=$4; ntm[$7]=$5 }  END{ for (r in copy) { s += copy[r]/ntm[r]}; print s }'
}


# I put the commands in the function to parallelize it
one_dpp() {
    base=/home/fuy2/data/piPipes/common/mm9
    ii=${i}dpp
    stat_file=${ii}.stat
    echo "Small RNA reads mapping stats" > $stat_file
    
    bed2=~/data/prepachytene/results/2014-11-21/SRA/${ii}/genome_mapping/Zamore.SRA.total.${ii}.testis.trimmed.x_rRNA.x_hairpin.mm9v1.all.bed2
    # bed2=00dpp.bed2
    ########################################################################
    # Intersect bed2 with piRNA clusters
    ########################################################################
    bed2_intersect_with_pirna_cluster=${ii}.intersect_with_piRNA_cluster.bed2
    bedtools intersect -b ${base}/mm9.piRNAcluster.bed.gz -a $bed2 -u > ${bed2_intersect_with_pirna_cluster}
    funky piRNAcluster ${bed2_intersect_with_pirna_cluster} >> $stat_file
    
    # echo -ne "piRNA_cluster\t" >> $stat_file
    # cat ${bed2_intersect_with_pirna_cluster} | awk '{ copy[$7]=$4; ntm[$7]=$5 }  END{ for (r in copy) { s += copy[r]/ntm[r]}; print s }' >> $stat_file
    for j in prepachytene hybrid pachytene; do
	bedtools intersect -b ${base}/piRNA.cluster.${j}.bed12.gz -a ${bed2_intersect_with_pirna_cluster} -u > ${ii}.intersect_with_piRNA_cluster.${j}.bed2
	bedtools intersect -b ${base}/piRNA.cluster.${j}.bed6.gz -a ${bed2_intersect_with_pirna_cluster} -u > ${ii}.intersect_with_piRNA_cluster.${j}.Exon.bed2
	funky "piRNAcluster.${j}" ${ii}.intersect_with_piRNA_cluster.${j}.bed2 >> $stat_file
	funky "piRNAcluster.${j}.Exon" ${ii}.intersect_with_piRNA_cluster.${j}.Exon.bed2 >> $stat_file
    done

    ########################################################################
    # Intersect bed2 with Genes
    ########################################################################
    bed2_intersect_with_gene=${ii}.intersect_with_Genes.bed2
    bedtools intersect -b ${base}/UCSC.refSeq.Genes.bed12.gz -a $bed2 -u > ${bed2_intersect_with_gene}
    funky "Genes" ${bed2_intersect_with_gene} >> $stat_file
    for j in Exon Intron; do
	bedtools intersect -b ${base}/UCSC.refSeq.${j}.bed.gz -a $bed2_intersect_with_gene -u > ${ii}.intersect_with_Genes.${j}.bed2
	funky "Genes.${j}" ${ii}.intersect_with_Genes.${j}.bed2 >> $stat_file
    done
    for j in 5UTR 3UTR; do
	bedtools intersect -b ${base}/UCSC.refSeq.${j}.bed.gz -a ${ii}.intersect_with_Genes.Exon.bed2 -u > ${ii}.intersect_with_Genes.Exon.${j}.bed2
	funky "Genes.Exon.${j}" ${ii}.intersect_with_Genes.Exon.${j}.bed2 >> $stat_file
    done

    ########################################################################
    # Intersect with RepeatMasker
    ########################################################################
    bed2_intersect_with_repeatmasker=${ii}.intersect_with_repeatmasker.bed2
    bedtools intersect -b ${base}/UCSC.RepeatMask.bed -a $bed2 -u > ${bed2_intersect_with_repeatmasker}
    funky RepeatMask ${bed2_intersect_with_repeatmasker} >> $stat_file
    # Features: ls -1 UCSC.RepeatMask.*.gz | grep -E -o '[^.]+.gz' | tr -d '.gz'
    for j in DNA LINE LTR rtRNA Satellite Simple_repeat SINE; do
	bedtools intersect -b ${base}/UCSC.RepeatMask.bed.${j}.gz -a ${bed2_intersect_with_repeatmasker} -u > ${ii}.intersect_with_repeatmasker.${j}.bed2
	funky "RepeatMask.$j" ${ii}.intersect_with_repeatmasker.${j}.bed2 >> $stat_file
    done
}


for i in 00 02 04 07 10 12 14 17 20 42; do
    one_dpp ${i} &
done
