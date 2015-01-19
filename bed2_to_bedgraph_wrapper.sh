for i in 00 02 04 07 10 12 14 17 20 42; do
# for i in 02; do
    n=$(grep 'genome mapping reads (-rRNA; +miRNA_hairpin)' ${i}dpp.basic_stats | grep -E -o '[0-9]+')
    scale=$(awk -v n=$n 'BEGIN{ print 1e6 / n }')
    echo $scale
    
    awk '{ l=$3-$2; if(l>=21 && l<=23 && $5==1) {print} }' ../2014-11-21/SRA/${i}dpp/genome_mapping/Zamore.SRA.total.${i}dpp.testis.trimmed.x_rRNA.x_hairpin.mm9v1.all.bed2 | sort -k1,1 -k2,2n > ${i}dpp.uniq.siRNA.bed2 && \
    	bed2_to_bedgraph.py ${i}dpp.uniq.siRNA.bed2 ${i}dpp.uniq.Watson.bedgraph.unscaled ${i}dpp.uniq.Crick.bedgraph.unscaled && \
    	awk -v scale=$scale 'BEGIN{FS=OFS="\t"} { $4=$4*scale; print }' ${i}dpp.uniq.Watson.bedgraph.unscaled > ${i}dpp.uniq.Watson.bedgraph && \
    	awk -v scale=$scale 'BEGIN{FS=OFS="\t"} { $4=$4*scale; print }' ${i}dpp.uniq.Crick.bedgraph.unscaled > ${i}dpp.uniq.Crick.bedgraph && \
    	bedGraphToBigWig ${i}dpp.uniq.Watson.bedgraph mm9.chrom.sizes ${i}dpp.siRNA.uniq.Watson.bigWig && \
    	bedGraphToBigWig ${i}dpp.uniq.Crick.bedgraph mm9.chrom.sizes ${i}dpp.siRNA.uniq.Crick.bigWig &
    
    awk '{ l=$3-$2; if(l>=21 && l<=23) {print} }' ../2014-11-21/SRA/${i}dpp/genome_mapping/Zamore.SRA.total.${i}dpp.testis.trimmed.x_rRNA.x_hairpin.mm9v1.all.bed2 | sort -k1,1 -k2,2n > ${i}dpp.all.siRNA.bed2 && \
    	bed2_to_bedgraph.py ${i}dpp.all.siRNA.bed2 ${i}dpp.all.Watson.bedgraph.unscaled ${i}dpp.all.Crick.bedgraph.unscaled && \
    	awk -v scale=$scale 'BEGIN{FS=OFS="\t"} { $4=$4*scale; print }' ${i}dpp.all.Watson.bedgraph.unscaled > ${i}dpp.all.Watson.bedgraph && \
    	awk -v scale=$scale 'BEGIN{FS=OFS="\t"} { $4=$4*scale; print }' ${i}dpp.all.Crick.bedgraph.unscaled > ${i}dpp.all.Crick.bedgraph && \
    	bedGraphToBigWig ${i}dpp.all.Watson.bedgraph mm9.chrom.sizes ${i}dpp.siRNA.all.Watson.bigWig && \
    	bedGraphToBigWig ${i}dpp.all.Crick.bedgraph mm9.chrom.sizes ${i}dpp.siRNA.all.Crick.bigWig &

done
