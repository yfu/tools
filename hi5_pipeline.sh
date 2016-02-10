#!/usr/bin/env bash

source /data/fuy2/piPipes/bin/piPipes_bash_functions
BEDTOOLS=/home/fuy2/data/piPipes/bin/bedtools_piPipes
FASTQ2INSERT=/data/fuy2/piPipes/bin/piPipes_fastq_to_insert
INSERT2BED2=/data/fuy2/piPipes/bin/piPipes_insertBed_to_bed2

while getopts "hi:Rc:g:o:sD" OPTION; do
    case $OPTION in
	h)usage && exit 0 ;;
	i)FQ=$(readlink -e $OPTARG) ;;
	c)CPU=$OPTARG ;;
	g)GENOME=$OPTARG ;;
	o)OUTDIR=$OPTARG ;;
	D)DEBUG=1 ;; ## DEBUG mode: skipping some steps to make testing easier
    esac
done

if [ $GENOME == "hi5" ]; then
    genome_MM=0
    GENOME_bowtie_idx=/data/fuy2/cl/shared/draft_genome/Hi5_69mer.2pe_2mp.fa
    rRNA_bowtie_idx=/data/fuy2/cl/shared/rrna.candidate.fasta
    hairpin_bowtie_idx=/data/fuy2/cl/shared/mirna/hairpin.fa
    hairpin_MM=0
fi

if [ -z "${GENOME}" ]; then
    echo "Please provide a genome assembly" >&2
    exit 1
fi

[ -z "${CPU}" ] && CPU=8

PREFIX=${FQ%.gz}
PREFIX=${PREFIX%.fastq}
PREFIX=${PREFIX%.fq}
PREFIX=$(basename ${PREFIX})

MM=0

[ -z "${OUTDIR}" ] && OUTDIR=.
mkdir -p ${OUTDIR} && cd ${OUTDIR}
table=${PREFIX}.basic_stats
mkdir -p input_files
mkdir -p genome_mapping

STEP=1

insert=${PREFIX}.insert
genome_bam=${PREFIX}.${genome_MM}mm.bam
log=${PREFIX}.${genome_MM}mm.log

[ ! -f .status.${STEP}.fq2insert ] && \
    ${FASTQ2INSERT} ${FQ} input_files/${PREFIX}.insert && \
    touch .status.${STEP}.fq2insert
[ ! -f .status.${STEP}.fq2insert ] && echo "fq2insert failed" "error"

STEP=$((STEP+1))

# Exclude rRNA reads:
rRNA_MM=2
echo "Mapping to rRNA, with $rRNA_MM mismatch(es) allowed"
INSERT=input_files/${PREFIX}.insert
x_rRNA_INSERT=input_files/${PREFIX}.x_rRNA.insert
rRNA_BED_LOG=input_files/${PREFIX}.rRNA.log

[ ! -f .status.${STEP}.rRNA_mapping ] && \
    totalReads=`awk '{a+=$2}END{printf "%d", a}' ${INSERT}` && echo $totalReads > .totalReads && \
    bowtie -r -S -v $rRNA_MM -k 1 -p $CPU \
	   --un $x_rRNA_INSERT \
	   ${rRNA_bowtie_idx} \
	   ${INSERT} \
	   1> /dev/null \
	   2> $rRNA_BED_LOG && \
    nonrRNAReads=`awk '{a+=$2}END{printf "%d", a}' ${x_rRNA_INSERT}` && echo $nonrRNAReads > .nonrRNAReads && \
    rRNAReads=$((totalReads-nonrRNAReads)) && \
    echo $rRNAReads > .rRNAReads && \
    touch .status.${STEP}.rRNA_mapping
[ ! -f .status.${STEP}.rRNA_mapping ] && echo "mapping to rRNA failed"
STEP=$((STEP+1))
# reading values from file, this is for resuming the job, which won't run the previous step
totalReads=`cat .totalReads`
rRNAReads=`cat .rRNAReads`
nonrRNAReads=`cat .nonrRNAReads`

# Hairpin mapping
# MIRNA_DIR=hairpins_mapping
PREFIX=$(basename ${x_rRNA_INSERT%.insert})
mkdir -p hairpin_mapping
echo "Mapping to microRNA Hairpin, with $hairpin_MM mismatch(es) allowed; only keep unique mappers"
x_rRNA_HAIRPIN_INSERT=hairpin_mapping/${PREFIX}.hairpin.insert # insert file storing reads that nonmappable to rRNA and mappable to hairpin
x_rRNA_x_hairpin_INSERT=hairpin_mapping/${PREFIX}.x_hairpin.insert # reads that nonmappable to rRNA or hairpin
x_rRNA_HAIRPIN_BED2=hairpin_mapping/${PREFIX}.hairpin.v${hairpin_MM}m1.bed2 # bed2 format with hairpin mapper, with the hairpin as reference
x_rRNA_HAIRPIN_LOG=hairpin_mapping/${PREFIX}.hairpin.v${hairpin_MM}m1.log # bed2 format with hairpin mapper, with the hairpin as reference
x_rRNA_HAIRPIN_BED2_LENDIS=hairpin_mapping/${PREFIX}.hairpin.v${hairpin_MM}m1.lendis # length distribution for hairpin mapper
x_rRNA_HAIRPIN_GENOME_BED2=genome_mapping/${PREFIX}.hairpin.${GENOME}v${genome_MM}a.bed2 # bed2 format with hairpin mapper, with genome as reference
x_rRNA_HAIRPIN_GENOME_LOG=genome_mapping/${PREFIX}.hairpin.${GENOME}v${genome_MM}a.log # log file for hairpin mapping
[ ! -f .status.${STEP}.hairpin_mapping ] && \
    bowtie -r --norc -v $hairpin_MM -m 1 --best --strata -p $CPU -S \
	   --al $x_rRNA_HAIRPIN_INSERT \
	   --un $x_rRNA_x_hairpin_INSERT \
	   $hairpin_bowtie_idx \
	   $x_rRNA_INSERT 2> ${x_rRNA_HAIRPIN_LOG} \
	| \
	samtools view -bSF 0x4 - 2>/dev/null | \
	$BEDTOOLS bamtobed -i - > ${PREFIX}.x_rRNA.hairpin.v${hairpin_MM}m1.bed && \
    ${INSERT2BED2} $x_rRNA_INSERT ${PREFIX}.x_rRNA.hairpin.v${hairpin_MM}m1.bed > $x_rRNA_HAIRPIN_BED2 && \
    rm -rf ${PREFIX}.x_rRNA.hairpin.v${hairpin_MM}m1.bed && \
    bed2lendis $x_rRNA_HAIRPIN_BED2 > $x_rRNA_HAIRPIN_BED2_LENDIS && \
    bowtie -r -v $genome_MM -a --best --strata -p $CPU \
	   -S \
	   ${GENOME_bowtie_idx} \
	   $x_rRNA_HAIRPIN_INSERT \
	   2> $x_rRNA_HAIRPIN_GENOME_LOG | \
	samtools view -uS -F0x4 - 2>/dev/null | \
	${BEDTOOLS} bamtobed -i - > ${PREFIX}.x_rRNA.hairpin.${GENOME}v${genome_MM}a.bed && \
    ${INSERT2BED2} $x_rRNA_HAIRPIN_INSERT ${PREFIX}.x_rRNA.hairpin.${GENOME}v${genome_MM}a.bed > $x_rRNA_HAIRPIN_GENOME_BED2 && \
    rm -rf ${PREFIX}.x_rRNA.hairpin.${GENOME}v${genome_MM}a.bed && \
    hairpinReads=`bedwc $x_rRNA_HAIRPIN_BED2` && echo $hairpinReads > .hairpinReads && \
    touch .status.${STEP}.hairpin_mapping
STEP=$((STEP+1))
hairpinReads=`cat .hairpinReads`

# Genome mapping
# INSERT=$x_rRNA_INSERT
INSERT=$(basename $x_rRNA_x_hairpin_INSERT)
# bed2 format storing all mappers for genomic mapping
GENOME_ALLMAP_BED2=genome_mapping/${INSERT%.insert}.${GENOME}v${genome_MM}.all.bed2 # all mapper in bed2 format
GENOME_ALLMAP_LENDIS=genome_mapping/${INSERT%.insert}.${GENOME}v${genome_MM}.all.lendis # all mapper in bed2 format
GENOME_ALLMAP_LOG=genome_mapping/${INSERT%.insert}.${GENOME}v${genome_MM}.all.log # log file
# bed2 format storing unique mappers for genomic mapping
GENOME_UNIQUEMAP_BED2=genome_mapping/${INSERT%.insert}.${GENOME}v${genome_MM}.unique.bed2
GENOME_UNIQUEMAP_LENDIS=genome_mapping/${INSERT%.insert}.${GENOME}v${genome_MM}.unique.lendis
# bed2 format storing unique mappers for genomic mapping and miRNA hairpin mapper
GENOME_UNIQUEMAP_HAIRPIN_BED2=genome_mapping/${INSERT%.insert}.${GENOME}v${genome_MM}.unique.+hairpin.bed2
# mapping insert file to genome
BEDTOOLS=/data/fuy2/piPipes/bin/bedtools_piPipes

echo "Mapping to genome, with ${genome_MM} mismatch(es) allowed"
[ ! -f .status.${STEP}.genome_mapping ] && \
    bowtie -r -v $genome_MM -a --best --strata -p $CPU \
	   --al  genome_mapping/${INSERT%.insert}.${GENOME}v${genome_MM}a.al.insert \
	   --un  genome_mapping/${INSERT%.insert}.${GENOME}v${genome_MM}a.un.insert \
	   -S \
	   $GENOME_bowtie_idx \
	   ${x_rRNA_x_hairpin_INSERT} \
	   2> $GENOME_ALLMAP_LOG | \
	samtools view -uS -F0x4 - 2>/dev/null | \
	${BEDTOOLS} bamtobed -i - > ${INSERT%.insert}.${GENOME}v${genome_MM}a.insert.bed && \
    ${INSERT2BED2} $x_rRNA_x_hairpin_INSERT ${INSERT%.insert}.${GENOME}v${genome_MM}a.insert.bed > ${GENOME_ALLMAP_BED2} && \
    rm -rf ${INSERT%.insert}.${GENOME}v${genome_MM}a.insert.bed && \
    touch .status.${STEP}.genome_mapping
[ ! -f .status.${STEP}.genome_mapping ] && echo "Genome mapping failed" "error"
STEP=$((STEP+1))

# Separate unique and multiple mappers
echo "Separating unique and multiple mappers"
[ ! -f .status.${STEP}.separate_unique_and_multiple ] && \
    awk 'BEGIN{OFS="\t"}{if ($5==1) print $0}' ${GENOME_ALLMAP_BED2} \
	1> ${GENOME_UNIQUEMAP_BED2}&& \
    totalMapCount=`bedwc ${GENOME_ALLMAP_BED2}` && echo $totalMapCount > .totalMapCount && \
    uniqueMapCount=`bedwc ${GENOME_UNIQUEMAP_BED2}` && echo $uniqueMapCount > .uniqueMapCount && \
    multipMapCount=$((totalMapCount-uniqueMapCount)) && echo $multipMapCount > .multipMapCount && \
    cat $x_rRNA_HAIRPIN_GENOME_BED2 ${GENOME_UNIQUEMAP_BED2} > $GENOME_UNIQUEMAP_HAIRPIN_BED2 && \
    touch .status.${STEP}.separate_unique_and_multiple
STEP=$((STEP+1))
totalMapCount=`cat .totalMapCount`
uniqueMapCount=`cat .uniqueMapCount`
multipMapCount=`cat .multipMapCount`

# TODO: uncommet the following line for the next version
awk '{ l[length($7)] += $4/$5 } END{ for(i in l) { print i "\t" l[i] } }' ${GENOME_ALLMAP_BED2} > ${GENOME_ALLMAP_LENDIS}
awk '{ if($5==1) { l[length($7)] += $4/$5 } } END{ for(i in l) { print i "\t" l[i] } }' ${GENOME_UNIQUEMAP_BED2} > ${GENOME_UNIQUEMAP_LENDIS}
plot_lendis.R ${GENOME_ALLMAP_LENDIS} ${GENOME_ALLMAP_LENDIS}.pdf
plot_lendis.R ${GENOME_UNIQUEMAP_LENDIS} ${GENOME_UNIQUEMAP_LENDIS}.pdf

# Generating stats
echo -e "Total reads:\t$totalReads" > ${table}
echo -e "rRNA reads:\t$rRNAReads" >> ${table}
echo -e "Hairpin-mapping reads:\t$hairpinReads" >> ${table}
echo -e "Total genome mappers:\t$totalMapCount" >> ${table}
echo -e "Total genome unique mappers:\t$uniqueMapCount" >> ${table}
echo -e "Total genome multiple mappers:\t$multipMapCount" >> ${table}
