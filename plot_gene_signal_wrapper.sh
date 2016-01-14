#!/usr/bin/env bash


# Test data:
# /data/fuy2/gfp_pirna/data/SRS.42A18AT1.br2.tr1.ox.ovary/bigWig_normalized_by_unique/SRS.42A18AT1.br2.tr1.ox.ovary.x_rRNA.x_hairpin.dm3_chengjianv0.all.piRNA.sorted.Watson.bigWig
# /data/fuy2/gfp_pirna/data/SRS.42A18AT1.br2.tr1.ox.ovary/bigWig_normalized_by_unique/SRS.42A18AT1.br2.tr1.ox.ovary.x_rRNA.x_hairpin.dm3_chengjianv0.all.piRNA.sorted.Crick.bigWig

watson=$1
crick=$2
chrominfo=~/data/piPipes/common/dm3_chengjian/dm3_chengjian.ChromInfo.txt

bigWigToBedGraph ${watson} stdout | awk '{ c=$1; if(c=="GSV6" || c=="nos-gal4-vp16" || c=="UASp-EGFP") {print} }' > ${watson}.bedGraph
bigWigToBedGraph ${crick} stdout | awk '{ c=$1; if(c=="GSV6" || c=="nos-gal4-vp16" || c=="UASp-EGFP") {print} }' > ${crick}.bedGraph

prefix=${watson%.Watson.bigWig}

for j in nos-gal4-vp16 GSV6 UASp-EGFP; do
    cons_l=$(grep ${j} ${chrominfo} | cut -f2)
    depth=${prefix}.depth
    cat ${watson}.bedGraph | awk -v c=${j} '{ if($1==c) {print} }' > tmp.watson
    cat ${crick}.bedGraph | awk -v c=${j} '{ if($1==c) {print} }' > tmp.crick
    
    bedgraph_to_depth.py tmp.watson 0 ${cons_l} > ${watson}.depth
    bedgraph_to_depth.py tmp.crick 0 ${cons_l} > ${crick}.depth
    paste ${watson}.depth ${crick}.depth | cut -f1,2,3,6 > ${prefix}.depth
    plot_gene_signal.R ${depth} "the name of the plot" ${prefix}.${j}.pdf
done
