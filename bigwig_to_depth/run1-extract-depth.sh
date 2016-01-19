
# for i in $(cat regions_of_interest.bed); do
SRA1=$(find /home/fuy2/data/pachytene/data/Xuan/SRA -maxdepth 1 -type d -name "Zamore*")
SRA2=$(find /home/fuy2/data/pachytene/data/Xuan/SRA2 -maxdepth 1 -type d -name "Zamore*")

for l in $SRA1 $SRA2; do
    # echo "#" ${l}
    watson=$(readlink -e ${l}/bigWig_normalized_by_unique/*.all.piRNA.sorted.Watson.bigWig)
    crick=$(readlink  -e ${l}/bigWig_normalized_by_unique/*.all.piRNA.sorted.Crick.bigWig)
    watson_bg=$(basename ${watson%.bigWig}).bedGraph
    crick_bg=$(basename ${crick%.bigWig}).bedGraph    
    # echo "bigWigToBedGraph ${watson} ${watson_bg}"
    # echo "bigWigToBedGraph ${crick} ${crick_bg}"
    while IFS='' read -r i || [[ -n "$line" ]]; do
    	chrom=$(echo "${i}" | cut -f1)
    	start=$(echo "${i}" | cut -f2)
    	end=$(echo "${i}" | cut -f3)
    	echo "bedgraph_to_depth.chrom.py ${watson_bg} ${chrom} ${start} ${end} > ${watson_bg%.bedGraph}.${chrom}.depth && \
    	bedgraph_to_depth.chrom.py ${crick_bg} ${chrom} ${start} ${end} > ${crick_bg%.bedGraph}.${chrom}.depth && \
	paste ${watson_bg%.bedGraph}.${chrom}.depth ${crick_bg%.bedGraph}.${chrom}.depth | cut -f1,2,3,6 > ${watson_bg%.Watson.bedGraph}.${chrom}.depth"
    done < regions_of_interest.bed
done
