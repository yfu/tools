# wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz

# Lift over requires 99%+ identity
# Looks like every lift over is perfect
min_macth=0.99
for i in prepachytene hybrid pachytene; do
    tmp=piRNA.cluster.${i}.bed12.tmp
    sorted_bed12=piRNA.cluster.${i}.bed12
    liftOver -minMatch=0.99 ~/data/piPipes/common/mm9/piRNA.cluster.${i}.bed12.gz mm9ToMm10.over.chain.gz ${tmp} piRNA.cluster.${i}.bed12.unMapped
    sort -k1,1 -k2,2n $tmp > $sorted_bed12 && rm ${tmp}

    sorted_bed6=piRNA.cluster.${i}.bed6
    tmp=piRNA.cluster.${i}.bed6.tmp
    bedtools bed12tobed6 -i ${sorted_bed12} > ${tmp}
    sort -k1,1 -k2,2n ${tmp} > ${sorted_bed6}
    sorted_bed6_uniq_name=piRNA.cluster.${i}.uniq_name.bed6
    awk 'BEGIN{OFS="\t"; i=1} { $4=$4 "." i++; print }' ${sorted_bed6} > ${sorted_bed6_uniq_name}

    bigbed=piRNA.cluster.${i}.bigBed
    bedToBigBed ${sorted_bed12} ~/data/shared/mm10/mm10.chrom.sizes ${bigbed}
    
done

cat piRNA.cluster.{prepachytene,hybrid,pachytene}.bed6 | sort -k1,1 -k2,2n  > piRNA.cluster.bed6
cat piRNA.cluster.{prepachytene,hybrid,pachytene}.uniq_name.bed6 | sort -k1,1 -k2,2n > piRNA.cluster.uniq_name.bed6
cat piRNA.cluster.{prepachytene,hybrid,pachytene}.bed12 | sort -k1,1 -k2,2n > piRNA.cluster.bed12
