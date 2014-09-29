for i in h3k27ac  h3k27me3  h3k4me1  h3k4me3  h3k9ac  h3k9me3; do
    cd ${i}
    for j in 0 4 8 12 16 20; do
	jj=$((j + 4))
	wig_file="*E${j}-${jj}_*_density.wig.gz"
	echo $wig_file
	bw_file=${j}to${jj}.${i}.bw
	output=${j}to${jj}.${i}.tab
	echo $bw_file
	echo $output
	
	zcat ${wig_file} | grep -v "mitochondrion" | wigToBigWig stdin ../dm3.chrom.sizes.modified ${bw_file}
	bigWigAverageOverBed ${bw_file} ~/data/piPipes/common/dm3/Brennecke.piRNAcluster.bed6.gz ${output}
    done
    cd ..
done
