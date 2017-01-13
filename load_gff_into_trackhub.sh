


# Annotations-2016-11-08.gff3 is the gff file exported from WebApollo on 2016-11-08
input=Annotations-2016-11-08.gff3
chromInfo=../../../genome/hi5_v0.8.2.ChromInfo.txt
prefix=${input%.gff3}
prefix=${prefix%.gff}

gffread ${input} -T -o ${prefix}.gtf

gtfToGenePred -geneNameAsName2 -infoOut=${prefix}.infoOut.txt -genePredExt ${prefix}.gtf ${prefix}.gp
awk -v OFS="\t" '{ $1=$12; print }' ${prefix}.gp > tmp.gp
genePredToBed tmp.gp stdout | sort -k1,1 -k2,2n > ${prefix}.bed12
bedToBigBed -type=bed12 -extraIndex=name ${prefix}.bed12 ${chromInfo} ${prefix}.bigBed

grep -v "^#" ${prefix}.infoOut.txt | awk '{printf "%s\t%s,%s,%s,%s,%s\n", $1,$2,$3,$8,$9,$10}' > ${prefix}.nameIndex.txt
ixIxx ${prefix}.nameIndex.txt ${prefix}.ix ${prefix}.ixx
