#!/usr/bin/env bash
# Input must contain only 1 record
# Input file must have the suffix of bed12
input=$1
prefix=${input%.bed12}

gtf=$HOME/data/shared/mm10/gencode.vM4.corrected_chrom_names.protein_coding.gtf
bed12=mm10.protein_coding.bed12
# gtfToGenePred -ignoreGroupsWithoutExons ${gtf} stdout -genePredExt -geneNameAsName2 | awk 'BEGIN{OFS="\t"} { $1=$1"_"$12; print }' | genePredToBed stdin stdout > ${bed12}

## liftOver -bedPlus=12 ~/data/shared/mm9ToMm10/piRNA.cluster.prepachytene.bed12 ~/data/shared/mm10/mm10ToHg19.over.chain piRNA.cluster.prepachytene.mm10ToHg19.bed12 piRNA.cluster.prepachytene.mm10ToHg19.bed12.unMapped
# Ignore overlapping features
## prefix=test
## input=${prefix}.bed12
output_u_tmp=${prefix}.closest.u.bed12.tmp
output_u=${prefix}.closest.u.bed12
output_d_tmp=${prefix}.closest.d.bed12.tmp
output_d=${prefix}.closest.d.bed12

# Use genome as the definition for up/downstream
bedtools closest -a ${input} -b ${bed12} -id -io -D ref -t first > ${output_u_tmp}
bedtools closest -a ${input} -b ${bed12} -iu -io -D ref -t first > ${output_d_tmp}
# Remember to remove the last column: distance
cat ${output_u_tmp} | awk '{ for(i=13; i<NF; i++) { printf $i "\t" }; print "" }' > ${output_u}
cat ${output_d_tmp} | awk '{ for(i=13; i<NF; i++) { printf $i "\t" }; print "" }' > ${output_d}

# lo: liftOver
lo=mm10ToHg19
output_u_lo=${prefix}.closest.u.${lo}.bed12
output_u_lo_unmapped=${prefix}.closest.u.${lo}.bed12.unMapped
output_d_lo=${prefix}.closest.d.${lo}.bed12
output_d_lo_unmapped=${prefix}.closest.d.${lo}.bed12.unMapped
chain=$HOME/data/shared/mm10/mm10ToHg19.over.chain
mm=0.1				# Use a lenient minMatch for interspecies liftOver
mb=0.1
liftOver -fudgeThick -minMatch=${mm} -minBlocks=${mb} ${output_u} ${chain} ${output_u_lo} ${output_u_lo_unmapped}
liftOver -fudgeThick -minMatch=${mm} -minBlocks=${mb} ${output_d} ${chain} ${output_d_lo} ${output_d_lo_unmapped}

# New start
start=$(cut -f3 ${output_u_lo})
start_chr=$(cut -f1 ${output_u_lo})
end=$(cut -f2 ${output_d_lo})
end_chr=$(cut -f1 ${output_d_lo})
name=$(cut -f4 ${input})_syntenic_flanking_region
if [[ -z "$start_chr" ]] || [[ -z "$end_chr" ]] || [[ "$start_chr" !=  "$end_chr" ]]; then
## if [[ "$start_chr" !=  "$end_chr" ]]; then    
    echo "Flanking genes are not in the same chromosome in the target genome. Skipping" > /dev/stderr
else
    echo -e "$start_chr\t$start\t$end\t$name" > ${prefix}.syntenic_flanking_region.bed6
fi
