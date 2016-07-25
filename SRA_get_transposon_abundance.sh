if [[ "$#" < 1 ]]; then
    echo -e "Usage: \nfor dUTP RNA-seq where read1 is antisense"
    echo -e "\tRSQ_get_transposon_rpm.sh your.bed2"
    exit 0
fi

input=$1
BEDTOOLS=~/data/piPipes/bin/bedtools_piPipes
prefix=${input%.bed2}
genome=mm10


awk 'BEGIN{ print "#transposon\ttype\tsense_all\tantisense_all\tsense_uniq\tantisense_uniq"; }' > ${prefix}.count.by_transposon
if [[ ${genome} == "mm10" ]]; then
    for i in DNA LINE LTR rRNA Satellite Simple_repeat SINE tRNA; do
	${BEDTOOLS} intersect -wo -f 0.5 -a ${input} -b ~/data/shared/mm10/UCSC.rmsk.${i}.bed | awk -v OFS="\t" -v type=${i} '{ c=$4/$5/$NF; if($6==$13){ sense_all[$11]+=c; if($5==1) { sense_uniq[$11]+=c } } else { antisense_all[$11]+=c; if($5==1) { antisense_uniq[$11]+=c } }; all[$11]=0; } END{ for(i in all) { s_all=sense_all[i]!=0?sense_all[i]:0; s_uniq=sense_uniq[i]!=0?sense_uniq[i]:0; as_all=antisense_all[i]!=0?antisense_all[i]:0; as_uniq=antisense_uniq[i]!=0?antisense_uniq[i]:0; print i, type, s_all, as_all, s_uniq, as_uniq  } }' 
    done >> ${prefix}.count.by_transposon
fi


# cat ${cmd} | parallel -j ${cpu}
# cat ${cmd2} | parallel -j ${cpu}
# cat ${cmd3} | parallel -j ${cpu}
# cat ${cmd4} | parallel -j ${cpu}

