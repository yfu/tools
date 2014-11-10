#!/usr/bin/env bash

# Given a insert file, bowite index, this script will output a sorted bam file, a sorted bed file and a sorted bed2 file.
# Usage: eval_mapping_rate.sh genome.idx insert_file
# Author: Yu Fu

if [ "$#" -ne 3 ]; then
    echo "Usage: eval_smallrna_mapping_rate.sh genome.idx insert_file suffix"
    exit 1
fi

genome_MM=0
CPU=8
BOWTIE_INDEX=$1
INPUT=$2
# Just a string for the suffix of the output filename
SUFFIX=$3

PREFIX=$(basename $INPUT)
PREFIX=${PREFIX%.insert}
PREFIX=${PREFIX%.inserts}

BAM_FILE=${PREFIX}_${genome_MM}mm.${SUFFIX}.sorted.bam
BED_FILE=${PREFIX}_${genome_MM}mm.${SUFFIX}.sorted.bed
BED2_FILE=${PREFIX}_${genome_MM}mm.${SUFFIX}.sorted.bed2

bowtie -r -v $genome_MM -a --best --strata -p $CPU \
    --al  ${PREFIX}.v${genome_MM}a.al.insert \
    --un  ${PREFIX}.v${genome_MM}a.un.insert \
    -S \
    $BOWTIE_INDEX \
    ${INPUT} \
    2> ${PREFIX}_${genome_MM}mm.${SUFFIX}.log | samtools view -b - -u | samtools sort - -m 5G -O bam -T ${RANDOM} > ${BAM_FILE}

bedtools bamtobed -i ${BAM_FILE} > ${BED_FILE}
~/data/piPipes/bin/piPipes_insertBed_to_bed2 ${INPUT} ${BED_FILE} > ${BED2_FILE}

echo -ne "Total number of reads: "
n_reads=$(awk 'BEGIN{sum=0} {sum+=$2} END{print sum}' < $INPUT)
echo $n_reads

echo -ne "Total number of reads mapped: "
n_reads_mapped=$(awk '{sum+=$4/$5} END{print sum}' $BED2_FILE)
echo $n_reads_mapped

echo -ne "Mapping rate of reads: "
awk -v n_reads_mapped=$n_reads_mapped -v n_reads=$n_reads 'BEGIN{ print n_reads_mapped/n_reads }'

echo ""

echo -ne "Total number of species: "
n_species=$(wc -l < $INPUT)
echo $n_species

echo -ne "Total number of species mapped: "
n_species_mapped=$(awk '{if($7 in d){} else {d[$7]=1; sum+=1} } END{print sum}' < $BED2_FILE)
echo $n_species_mapped

echo -ne "Mapping rate of species: "
awk -v n_species_mapped=$n_species_mapped -v n_species=$n_species 'BEGIN{ print n_species_mapped/n_species }'


