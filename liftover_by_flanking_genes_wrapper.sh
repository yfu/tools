
for i in prepachytene hybrid pachytene; do
    cluster=~/data/shared/mm9ToMm10/piRNA.cluster.${i}.bed12
    mkdir -p ${i}
    while read line; do
	prefix=$(echo $line | awk '{ print $4 }')
	bed12=${i}/${prefix}.bed12
	echo "${line}" > ${bed12}
    echo "bash ./liftover_by_flanking_genes.sh ${bed12}"
    done < ${cluster}
done
