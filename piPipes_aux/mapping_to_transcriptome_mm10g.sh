#!/usr/bin/env bash

source ~/repo/tools/piPipes_aux/functions.sh
source ~/repo/tools/piPipes_aux/global.sh
GENOME=mm10g
COMMON_FOLDER=~/data/piPipes/common/${GENOME}

while getopts "hd:i:c:o:g:B:xvLD" OPTION; do
    case $OPTION in
	h)usage && exit 0 ;;
	d)PR=`readlink -f ${OPTARG}`;; # the directory of pipeline results
	## g)export GENOME=${OPTARG};;
	o)OUTDIR=`readlink -f ${OPTARG}` ;;
	c)CPU=$OPTARG ;;
    esac
done


xrRNA_LEFT_FQ=$(readlink -e ${PR}/input_read_files/*.x_rRNA.1.fq)
xrRNA_RIGHT_FQ=$(readlink -e ${PR}/input_read_files/*.x_rRNA.2.fq)
PREFIX=$(basename $PR)
PREFIX=${PREFIX%.x_rRNA.1.fq}
# xrRNA_LEFT_FQ=test.r1.fq
# xrRNA_RIGHT_FQ=test.r2.fq
echo ${PREFIX}
echo ${xrRNA_LEFT_FQ}
echo ${xrRNA_RIGHT_FQ}

[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
[ ! -z ${OUTDIR} ] || OUTDIR=$PWD # if -o is not specified, use current directory


PREFIX=`echo -e "$(basename ${xrRNA_LEFT_FQ})\n$(basename ${xrRNA_RIGHT_FQ})" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'` && export PREFIX=${PREFIX%.*}
INDIR=
LIGATIONLIB=0
if [ "${LIGATIONLIB}" == 1 ]; then
    LIBRARY_TYPE="fr-secondstrand";
    END_TO_REVERSE_STRAND=2
    SENSE_HTSEQ_OPT="yes";
    ANTISENSE_HTSEQ_OPT="reverse";
    EXPRESS_OPTION_PE="--fr-stranded"
    EXPRESS_OPTION_SE="--f-stranded"
else
    LIBRARY_TYPE="fr-firststrand"
    END_TO_REVERSE_STRAND=1
    SENSE_HTSEQ_OPT="reverse";
    ANTISENSE_HTSEQ_OPT="yes";
    EXPRESS_OPTION_PE="--rf-stranded"
    EXPRESS_OPTION_SE="--r-stranded"
fi

echo2 "Determining the version of fastQ using SolexaQA"

# determine version of fastq used, using a modified SolexaQA.pl
PHRED_SCORE=`Solexa_QA.pl "${xrRNA_LEFT_FQ}"`
case ${PHRED_SCORE} in
    solexa)bowtie2PhredOption="--solexa-quals" && STARoutQSconversionAdd="-31" ;; # Solexa+64, raw reads typically (-5, 40)
    illumina)bowtie2PhredOption="--phred64" && STARoutQSconversionAdd="-31" ;; # Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
    sanger)bowtie2PhredOption="--phred33" && STARoutQSconversionAdd="0" ;; # Phred+33,  raw reads typically (0, 40) (http://en.wikipedia.org/wiki/FASTQ_format)
    *)echo2 "unable to determine the fastq version. Using sanger..." "warning";;
esac

echo2 "Mapping to genes and transposon directly with Bowtie2"
TRANSCRIPTOME_INDEX="/data/fuy2/piPipes/common/mm10g/Bowtie2Index/gene+cluster+repBase"
TRANSCRIPTOME_NAME="gene+cluster+repBase"
# TRANSCRIPTOME_SIZES=$BOWTIE2_INDEXES/${TRANSCRIPTOME_INDEX}.sizes
DIRECTMAPPING_DIR=${OUTDIR}/transcriptome_mapping
mkdir -p ${DIRECTMAPPING_DIR}
[ ! -f ${OUTDIR}/.status.direct_mapping ] && \
    bowtie2 -x ${TRANSCRIPTOME_INDEX} \
	    -1 ${xrRNA_LEFT_FQ} \
	    -2 ${xrRNA_RIGHT_FQ} \
	    -q \
	    $bowtie2PhredOption \
	    -a \
	    -X 800 \
	    --no-mixed \
	    --quiet \
	    -p $CPU \
	    2> ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_NAME}.log | \
	samtools view -bS - > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_NAME}.bam && \
    samtools sort -o -@ $CPU ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_NAME}.bam ${DIRECTMAPPING_DIR}/foo | \
	${BEDTOOLS} bamtobed -i - | \
	awk -v etr=$END_TO_REVERSE_STRAND -v MAPQ=10 'BEGIN{FS=OFS="\t"}{l=split($4,arr,""); if (arr[l]==etr) $6=($6=="+"?"-":"+"); if ($5 > MAPQ) print $0}' > ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_NAME}.sorted.unique.bed && \
    touch .status.direct_mapping

# PR=/data/fuy2/pachytene/data/Xuan/Zamore.RSQ.B6.0dpp.rep1
# PREFIX=Zamore.RSQ.B6.0dpp.rep1

CUFFLINKS_DIR=${PR}/cufflinks_output
MapMass=$(grep 'Normalized Map Mass' $CUFFLINKS_DIR/${PREFIX}.cufflinks.log | cut -d' ' -f4)
export NormScale=`echo $MapMass | awk '{printf "%f",1000000.0/$1}'`

eXpressBATCH=10
[ ! -f .${OUTDIR}/.status.eXpress_quantification ] && \
express \
    $EXPRESS_OPTION_PE \
    -B $eXpressBATCH \
    -o $DIRECTMAPPING_DIR \
    --no-update-check \
    --library-size ${MapMass%.*} \
    $COMMON_FOLDER/${GENOME}.${TRANSCRIPTOME_NAME}.fa \
    ${DIRECTMAPPING_DIR}/${PREFIX}.${TRANSCRIPTOME_NAME}.bam \
    1>&2 2> $DIRECTMAPPING_DIR/${PREFIX}.${TRANSCRIPTOME_NAME}.eXpress.log && \
    awk -v depth=$NormScale 'BEGIN{OFS="\t"; getline; print}{$8*=depth; print}' $DIRECTMAPPING_DIR/results.xprs > \
	$DIRECTMAPPING_DIR/results.xprs.normalized && \
    touch .${JOBUID}.status.${STEP}.eXpress_quantification
