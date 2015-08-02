#!/usr/bin/env bash

# Given a prefix, this script outputs the stats
genome=hg19
prefix=$1
sep="------------------------------------------------------------------------"
echo "Stats for ${prefix}"
echo ${sep}
if [[ ${FQ##*.} == gz ]]; then
    total_input=$(zcat ${prefix} | wc -l | awk '{ print $1/4 }')
else
    total_input=$(cat ${prefix} | wc -l | awk '{ print $1/4 }')
fi

echo "Raw reads"
printf "Total # reads:\t%.1f\n" ${total_input}
total_unmapped=$(cat ${prefix}.hg19.Unmapped.out.mate1 | wc -l | awk '{ print $1/4 }' )
printf "Total # unmapped reads to start with:\t%.1f\n" ${total_unmapped}

echo ${sep}
echo "5' intron mapping:"
total_5intron_reads=$(samtools view ${prefix}.hg19.map_to_5intron.bam | wc -l)
total_5intron_s_reads=$(cat ${prefix}.hg19.map_to_5intron.s.sam | wc -l)
total_5intron_as_reads=$(cat ${prefix}.hg19.map_to_5intron.as.sam | wc -l)
printf "# reads containing 5' intron sequence:\t%.1f\n " ${total_5intron_reads}
printf "    sense:\t%.1f\n " ${total_5intron_s_reads}
printf "    antisense:\t%.1f\n" ${total_5intron_as_reads}

echo ${sep}
echo "Usable reads:"
mappable_bef=$(cat ${prefix}.before.map_to_genome.bed | wc -l)
mappable_ali=$(cat ${prefix}.align.map_to_genome.bed  | wc -l)
mappable_join_bef_ali=$(cat ${prefix}.before+after.map_to_genome.bed | wc -l)
printf "# reads with bef mappable:\t%.1f\n" ${mappable_bef}
printf "# reads with ali mappable:\t%.1f\n" ${mappable_ali}
printf "# reads with both bef and ali mappable:%.1f\n" ${mappable_join_bef_ali}
first_filter=$(cat ${prefix}.before+align.map_to_genome.f1.bed | wc -l)
second_filter=$(cat ${prefix}.before+align.map_to_genome.f2.bed | wc -l)
printf "# reads with both bef and ali on the same strand on the same chromosome:\t%.1f\n" ${first_filter}
printf "# reads with both bef and ali in the same intron:\t%.1f\n" ${second_filter}

echo ${sep}
echo "Extract lariat supporting reads"
put_lar_supporting=$(cat ${prefix}.before+align.map_to_genome.f2.before.tab | wc -l)
put_lar_count=$(grep '>' ${prefix}.put_lar.fa | wc -l)
lariat_mapping=$(samtools view ${prefix}.map_to_put_lar.bam | wc -l)
crsbp_count=$(samtools view ${prefix}.map_to_put_lar.crsbp.bam | wc -l)
printf "# lariat supporting reads:\t%.1f\n" ${put_lar_supporting}
printf "# lariats:\t%.1f\n" ${put_lar_count}

printf "Map all unmapped reads to lariats, # mappable reads:\t%.1f\n" ${lariat_mapping}
printf "Out of the above reads, # reads crossing the BP with at least 5nt on left and right of BP:\t%.1f\n" ${crsbp_count}
echo ${sep}
echo "Final:"
final_lar_count=$(cat ${prefix}.final.lar.bp.bed | wc -l)
printf "# lariats found:\t%.1f\n" ${final_lar_count}



