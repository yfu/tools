#!/usr/bin/env bash

cpu=4
# Given a bam file, this script outputs the counts of tags in genes, transposons and clusters in dm3
bam=$1
prefix=${bam%.bam}
prefix=$(basename ${prefix})

# bedtools bamtobed -i ${bam} > ${prefix}.bed
# samtools view ${bam} | wc -l > ${prefix}.input.total_reads
nf=$(cat ${prefix}.input.total_reads | awk '{ print 1e6/$1 }')

declare -A t
t["gene"]=/data/fuy2/piPipes/common/dm3/UCSC.refSeq.Genes.bed12.gz
t["transposon"]=/data/fuy2/piPipes/common/dm3/Zamore.transposon.bed.gz
t["cluster"]=/data/fuy2/piPipes/common/dm3/Brennecke.piRNAcluster.bed6.gz
t["threeprimeutr"]=/data/fuy2/piPipes/common/dm3/UCSC.refSeq.3UTR.bed.gz

commands=${prefix}.parallel
if [ -f ${commands} ]; then
    rm ${commands}
fi
for i in gene transposon cluster threeprimeutr; do
    echo "cat ${prefix}.bed | bedtools intersect -a - -b ${t[${i}]} -wa -wb | bedtools groupby -i - -g 1,2,3,4,5,6 -c 10 -o collapse | awk '{ split(\$7, a, \",\"); l=length(a); for(i in a) { c[ a[i] ] += 1/l; } } END{ for(i in c) { print i \"\t\" c[i] } }' > ${prefix}.${i}.raw_count" >> ${commands}
done

cat ${commands} | parallel --progress -j ${cpu}

cat ${prefix}.transposon.raw_count | awk '{ match($1, /^(.+)\.[0-9]+$/, a); if(a[1]=="") { print } c[ a[1] ] += $2 } END{ for(i in c) { print i "\t" c[i] } }' > ${prefix}.transposon_family.raw_count

cat ${prefix}.{gene,transposon_family,cluster,threeprimeutr}.raw_count | awk -v nf=$nf 'BEGIN{ c["OTHERS"]=0; }{ if(FNR==NR) { g[$1]=$2; c[$2]=0 } else { if($1 in g) { c[g[$1]]+=$2; } else { c[$1]+=$2 } } } END{ for(i in c) { print i "\t" c[i] "\t" c[i] * nf } }' dm3.txid_gid - > ${prefix}.count

