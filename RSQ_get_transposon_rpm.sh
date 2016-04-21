#!/usr/bin/env bash

# Given a (soft-linked) SAM file, output the # reads that overlap with transposon annotations
sam=$1

prefix=${sam%.sam}
prefix=${prefix%.bam}
cpu=8
genome=mm10
## convert bam (sam) to bed6 and put r1 to the sense strand
samtools view -b -F 0x100 ${sam} | bedtools bamtobed -bed12 -i - | bedtools bed12tobed6 -i /dev/stdin | awk -v OFS="\t" '{ l=length($4); r=substr($4, l, 1); if(r==1) { $6 = $6=="+"?"-":"+"; print } else { print } }' | sort -k1,1 -k2,2n --parallel=${cpu} > ${prefix}.bed6

cmd=${prefix}.transposon_intersection.cmd
cmd2=${prefix}.transposon_intersection.cmd2
cmd3=${prefix}.transposon_intersection.cmd3
cmd4=${prefix}.transposon_intersection.cmd4

echo "" > ${cmd}
echo "" > ${cmd2}
echo "" > ${cmd3}
echo "" > ${cmd4}

if [[ ${genome} == "mm10" ]]; then
    for i in DNA LINE LTR rRNA Satellite Simple_repeat SINE tRNA; do
	# echo "bedtools intersect -sorted -s -a ${prefix}.bed6 -b ~/data/shared/mm10/UCSC.rmsk.${i}.sorted.bed -wa -wb -f 0.9 | bedtools groupby -i - -g 10 -c 4 -o count_distinct > ${prefix}.${i}.count" >> ${cmd}
	echo "bedtools intersect -sorted -s -a ${prefix}.bed6 -b ~/data/shared/mm10/UCSC.rmsk.${i}.sorted.bed -wa -wb -f 0.9 | sort -k1,1 -k2,2n | bedtools groupby -i - -g 7,8,9,10 -c 4 -o count_distinct > ${prefix}.${i}.count " >> ${cmd}
	echo "cat ${prefix}.${i}.count | awk '{ a[\$4]+=\$5 } END{ for(i in a) { print i \"\t\" a[i] } }' > ${prefix}.${i}.count.by_transposon" >> ${cmd2}
	# echo "cat ${prefix}.${i}.count | awk '{ a[\$4]+=\$5 } END{ for(i in a) { print i \"\t\" a[i] } }' > ${prefix}.${i}.count.by_transposon" >> ${cmd2}
	echo "awk '{ s+=\$2 } END{ print s }' ${prefix}.${i}.count.by_transposon > ${prefix}.${i}.count.by_transposon_family" >> ${cmd3}
    done

    if [ -f ${prefix}.all.count.by_transposon_family ]; then
	rm ${prefix}.all.count.by_transposon_family
	touch ${prefix}.all.count.by_transposon_family
    fi
    
    for i in DNA LINE LTR rRNA Satellite Simple_repeat SINE tRNA; do
	echo "awk -v i=$i '{ print i \"\t\" \$1 }' ${prefix}.${i}.count.by_transposon_family >> ${prefix}.all.count.by_transposon_family" >> ${cmd4}
    done
fi

cat ${cmd} | parallel -j ${cpu}
cat ${cmd2} | parallel -j ${cpu}
cat ${cmd3} | parallel -j ${cpu}
cat ${cmd4} | parallel -j ${cpu}

# Consolidate the numbers
