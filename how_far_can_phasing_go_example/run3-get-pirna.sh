#!/usr/bin/env bash



root=/data/fuy2/gfp_pirna/data

# grep 'chr2R' ~/data/piPipes/common/dm3_chengjian/dm3_chengjian.ChromInfo.txt > dm3_chengjian.chr2R.ChromInfo.txt
## get bed2 files for unique mappers
# for i in 1 2 3; do
#     echo "awk -v OFS=\"\t\" '{ if(\$5==1) { print } }' /data/fuy2/gfp_pirna/data/SRS.GS1BT1.br${i}.tr1.ox.ovary/genome_mapping/SRS.GS1BT1.br${i}.tr1.ox.ovary.x_rRNA.x_hairpin.dm3_chengjianv0.all.piRNA.bed2 > SRS.GS1BT1.br${i}.tr1.ox.ovary.x_rRNA.x_hairpin.dm3_chengjianv0.uniq.piRNA.bed2"
# done




print_intersect_cmd() {
    my_source_bed2=$1
    my_target_bed=$2
    my_prefix=${source_bed2%.bed2}
    # echo "bedtools intersect -a ${my_source_bed2} -b ${my_target_bed} | awk -v OFS=\"\t\" '{ if(\$5==1) {print} }' > ${my_prefix}.uniq.${my_target_bed}"
    echo "gppc_3_to_5.sh dm3_chengjian.chr2R.ChromInfo.txt ${my_prefix}.uniq.${my_target_bed}"
}



for i in SRS.GS1BT1.br1.tr1.ox.ovary  SRS.GS1BT1.br2.tr1.ox.ovary  SRS.GS1BT1.br3.tr1.ox.ovary; do
    for j in $(seq 0 500 10000); do
	offset=-${j}
	target_bed_1=GSV6_upstream_in_zip.dm3.2k.offset${offset}.bed
	offset=${j}
	target_bed_2=GSV6_downstream_in_zip.dm3.2k.offset${offset}.bed
	source_bed2=${i}.x_rRNA.x_hairpin.dm3_chengjianv0.uniq.piRNA.bed2
	print_intersect_cmd ${source_bed2} ${target_bed_1}	
	print_intersect_cmd ${source_bed2} ${target_bed_2}	
    done
done

