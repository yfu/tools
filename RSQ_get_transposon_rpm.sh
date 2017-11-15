#!/usr/bin/env bash


BEDTOOLS=/usr/local/bin/bedtools

if [[ "$#" < 1 ]]; then
    echo -e "Usage: \nfor dUTP RNA-seq where read1 is antisense\n\tRSQ_get_transposon_rpm.sh your.bam"
    echo -e "\tRSQ_get_transposon_rpm.sh your.bam reverse"
    echo -e "For Degrodome seq where read2 is antisense"
    echo -e "\tRSQ_get_transposon_rpm.sh your.bam yes"
    exit 0
fi
	  
# Given a (soft-linked) SAM file, output the # reads that overlap with transposon annotations
sam=$1
# The second option specifies which strand to use
# For dUTP RNA-seq, just skip this option
# For Degrodome-seq, use yes
if [[ -z "$2" ]]; then
    strand_to_reverse=1
elif [[ $2 == "reverse" ]]; then
    strand_to_reverse=1
elif [[ $2 == "yes" ]]; then
    strand_to_reverse=2
fi

prefix=${sam%.sam}
prefix=${prefix%.bam}

cpu=8
genome=mm10
## convert bam (sam) to bed6 and put r1 to the sense strand
## For each read, I only keep the best alignment result
status_f=${prefix}.convert_to_bed6.ok
samtools view -b -F 0x100 ${sam} | ${BEDTOOLS} bamtobed -bed12 -i - | ${BEDTOOLS} bed12tobed6 -i /dev/stdin | awk -v strand_to_reverse=${strand_to_reverse} -v OFS="\t" '{ l=length($4); r=substr($4, l, 1); if(r==strand_to_reverse) { $6 = $6=="+"?"-":"+"; print } else { print } }' | sort -k1,1 -k2,2n --parallel=${cpu} > ${prefix}.bed6

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
	echo "${BEDTOOLS} intersect -sorted -s -a ${prefix}.bed6 -b ~/data/shared/mm10/UCSC.rmsk.${i}.sorted.bed -wa -wb -f 0.9 | sort -k1,1 -k2,2n | ${BEDTOOLS} groupby -i - -g 7,8,9,10 -c 4 -o count_distinct > ${prefix}.${i}.count " >> ${cmd}
	echo "cat ${prefix}.${i}.count | awk '{ a[\$4]+=\$5 } END{ for(i in a) { print i \"\t\" a[i] } }' > ${prefix}.${i}.count.by_transposon" >> ${cmd2}
	# echo "cat ${prefix}.${i}.count | awk '{ a[\$4]+=\$5 } END{ for(i in a) { print i \"\t\" a[i] } }' > ${prefix}.${i}.count.by_transposon" >> ${cmd2}
	echo "awk '{ s+=\$2 } END{ print s }' ${prefix}.${i}.count.by_transposon > ${prefix}.${i}.count.by_transposon_class" >> ${cmd3}
    done

    if [ -f ${prefix}.all.count.by_transposon_class ]; then
	rm ${prefix}.all.count.by_transposon_class
	touch ${prefix}.all.count.by_transposon_class
    fi
    
    for i in DNA LINE LTR rRNA Satellite Simple_repeat SINE tRNA; do
	echo "awk -v i=$i '{ print i \"\t\" \$1 }' ${prefix}.${i}.count.by_transposon_class >> ${prefix}.all.count.by_transposon_class" >> ${cmd4}
    done
fi

cat ${cmd} | parallel -j ${cpu}
cat ${cmd2} | parallel -j ${cpu}
cat ${cmd3} | parallel -j ${cpu}
cat ${cmd4} | parallel -j ${cpu}

for i in DNA LINE LTR rRNA Satellite Simple_repeat SINE tRNA; do
    awk -v type=${i} -v OFS="\t" '{ print $0, type }' ${prefix}.${i}.count.by_transposon
done | sort -k3,3 -k1,1 > ${prefix}.count.by_transposon_member

mkdir ${prefix}.debugging
for i in DNA LINE LTR rRNA Satellite Simple_repeat SINE tRNA; do
    mv ${prefix}.${i}.count.by_transposon ${prefix}.${i}.count.by_transposon_class ${prefix}.${i}.count ${prefix}.debugging/
done

