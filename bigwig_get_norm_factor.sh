
cat ~/data/piPipes/common/mm9/mm9.ChromInfo.txt | awk 'BEGIN{OFS="\t"} { print $1, "1", $2, $1}' > mm9.bed
cat ~/data/shared/mm10/mm10.chrom.sizes | grep -v "random" | grep -v "chrUn" | awk 'BEGIN{OFS="\t"} { print $1, "1", $2, $1}' > mm10.bed

for i in ../../../data/*signal.bigWig; do
    bn=$(basename ${i})
    output=${bn%.bigWig}.norm_factor
    if [[ ${i} =~ mm9 ]]; then
	echo "bigWigAverageOverBed ${i} mm9.bed stdout | awk '{ s+=\$4 } END{ print s / 1e6 }' > ${output}"
    elif [[ ${i} =~ mm10 ]]; then
	echo "bigWigAverageOverBed ${i} mm10.bed stdout | awk '{ s+=\$4 } END{ print s / 1e6 }' > ${output}"
    fi
done
