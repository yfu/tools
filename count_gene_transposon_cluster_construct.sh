#!/usr/bin/env bash

cpu=4
# Given a bam file, this script outputs the counts of tags in genes, transposons and clusters in dm3
bed=$1
nf=$2
prefix=${bed%.bed}
prefix=${prefix%.bed12}
prefix=$(basename ${prefix})

# bedtools bamtobed -i ${bam} > ${prefix}.bed
# samtools view ${bam} | wc -l > ${prefix}.input.total_reads
# nf=$(cat ${prefix}.input.total_reads | awk '{ print 1e6/$1 }')

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
    echo "cat ${bed} | bedtools intersect -split -s -a - -b ${t[${i}]} -wa -wb | bedtools groupby -i - -g 1,2,3,4,5,6 -c 16 -o collapse | awk '{ split(\$7, a, \",\"); l=length(a); for(i in a) { c[ a[i] ] += 1/l; } } END{ for(i in c) { print i \"\t\" c[i] } }' > ${prefix}.${i}.S.raw_count" >> ${commands}
    echo "cat ${bed} | bedtools intersect -split -S -a - -b ${t[${i}]} -wa -wb | bedtools groupby -i - -g 1,2,3,4,5,6 -c 16 -o collapse | awk '{ split(\$7, a, \",\"); l=length(a); for(i in a) { c[ a[i] ] += 1/l; } } END{ for(i in c) { print i \"\t\" c[i] } }' > ${prefix}.${i}.AS.raw_count" >> ${commands}
done

echo "cat ${bed} | awk '{ if(\$6==\"+\") {print} }' | awk '{ if(\$1==\"GSV6\" || \$1==\"nos-gal4-vp16\" || \$1==\"UASp-EGFP\") {s[\$1]+=1} } END{ for(i in s) { print i \"\t\" s[i] } }' > ${prefix}.construct.Watson.raw_count && awk -v nf=$nf -v OFS=\"\t\" '{ print \$1, \$2, \$2*nf }' ${prefix}.construct.Watson.raw_count > ${prefix}.construct.Watson.count" >> ${commands}
echo "cat ${bed} | awk '{ if(\$6==\"-\") {print} }' | awk '{ if(\$1==\"GSV6\" || \$1==\"nos-gal4-vp16\" || \$1==\"UASp-EGFP\") {s[\$1]+=1} } END{ for(i in s) { print i \"\t\" s[i] } }' > ${prefix}.construct.Crick.raw_count && awk -v nf=$nf -v OFS=\"\t\" '{ print \$1, \$2, \$2*nf }' ${prefix}.construct.Crick.raw_count > ${prefix}.construct.Crick.count" >> ${commands}

# cat ${commands} | parallel --progress -j ${cpu}

cat ${prefix}.transposon.S.raw_count | awk '{ match($1, /^(.+)\.[0-9]+$/, a); if(a[1]=="") { print } c[ a[1] ] += $2 } END{ for(i in c) { print i "\t" c[i] } }' > ${prefix}.transposon_family.S.raw_count
cat ${prefix}.transposon.AS.raw_count | awk '{ match($1, /^(.+)\.[0-9]+$/, a); if(a[1]=="") { print } c[ a[1] ] += $2 } END{ for(i in c) { print i "\t" c[i] } }' > ${prefix}.transposon_family.AS.raw_count



## Consolidate the results
cat ${prefix}.${i}.S.raw_count | awk -v nf=$nf -v type=gene -v strand=S 'BEGIN{ c["OTHERS"]=0; }{ if(FNR==NR) { g[$1]=$2; c[$2]=0 } else { if($1 in g) { c[g[$1]]+=$2; } else { c[$1]+=$2 } } } END{ for(i in c) { print i "\t" type "\t" strand "\t" c[i] "\t" c[i] * nf } }' dm3.txid_gid - > ${prefix}.count
cat ${prefix}.${i}.AS.raw_count | awk -v nf=$nf -v type=gene -v strand=AS 'BEGIN{ c["OTHERS"]=0; }{ if(FNR==NR) { g[$1]=$2; c[$2]=0 } else { if($1 in g) { c[g[$1]]+=$2; } else { c[$1]+=$2 } } } END{ for(i in c) { print i "\t" type "\t" strand "\t" c[i] "\t" c[i] * nf } }' dm3.txid_gid - >> ${prefix}.count

for i in transposon_family cluster threeprimeutr; do
    for j in S AS; do
	cat ${prefix}.${i}.${j}.raw_count | awk -v nf=$nf -v type=${i} -v strand=${j} '{ if($1 in g) { c[g[$1]]+=$2; } else { c[$1]+=$2 } } END{ for(i in c) { print i "\t" type "\t" strand "\t" c[i] "\t" c[i] * nf } }' -
    done 
done >> ${prefix}.count


awk -v nf=$nf '{ print $1 "\t" "construct" "\t" "Watson" "\t" $2 "\t" $2 * nf }' ${prefix}.construct.Watson.raw_count >> ${prefix}.count
awk -v nf=$nf '{ print $1 "\t" "construct" "\t" "Crick" "\t" $2 "\t" $2 * nf }' ${prefix}.construct.Crick.raw_count >> ${prefix}.count
