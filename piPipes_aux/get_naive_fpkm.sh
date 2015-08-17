#!/usr/bin/bash

# Given a output dir from piPipes, this script creates a folder called naive_counts and output the raw counts and them FPKM into them

REPO=/home/fuy2/repo/tools
source ${REPO}/piPipes_aux/functions.sh
COUNTER=${REPO}/htseq_bam_to_count.sh

while getopts "hd:i:c:o:g:B:xvLD" OPTION; do
    case $OPTION in
	h)usage && exit 0 ;;
	d)PR=`readlink -f ${OPTARG}`;; # the directory of pipeline results
	c)CPU=$OPTARG ;;
	g)GENOME=$OPTARG ;;
    esac
done

if [[ -z "$PR" ]]; then
    echo "Usage: get_naive_fpkm.sh my_pipipes_dir"
    exit 0
fi

if [[ -z "$CPU" ]]; then
    CPU=8
fi

total_uniq=$(grep 'genom.*unique_mapper_reads' ${PR}/*.basic_stats | cut -f2)
nf=$(echo ${total_uniq} | awk '{ print 1e6/$1 }')
if [[ -d 'naive_counts' ]]; then
    echo2 "naive_counts dir already exists. Exit."
    exit 1
fi

bam=$(readlink -e ${PR}/genome_mapping/*sorted.bam)
echo2 "Naive count for ${bam}..."

if [[ ${GENOME} == "hg19" ]]; then
    GTF=~/data/shared/hg19/gencode_r19/gencode.v19.annotation.gtf
    med_len=/data/fuy2/shared/hg19/hg19.gencode_r19.med_gene_len
elif [[ ${GENOME} == "mm10" ]]; then
    GTF=~/data/shared/mm10/gencode.vM4.corrected_chrom_names.gtf
    med_len=/data/fuy2/shared/mm10/mm10.gencode_vM4.med_gene_len
fi
########
# test #
########
# GTF=/home/fuy2/data/shared/mm10/mm10.test10000.gtf
${COUNTER} ${GTF} ${bam}

mkdir -p ${PR}/naive_counting
dest=${PR}/naive_counting
mv ${PR}/genome_mapping/*htseq* ${PR}/genome_mapping/*sorted_by_name* ${PR}/genome_mapping/*.ok ${dest}

raw_count=$(readlink -e ${dest}/*.htseq)
count=${raw_count}.naive_counting
grep -v '^__' ${raw_count} | awk -v nf=${nf} 'BEGIN{print "# gene\traw\tppm\tfpkm\tmed_len"} { if(FNR==NR) { l[$1]=$2 } else { print $1 "\t" $2 "\t" $2*nf "\t" $2*nf/l[$1]*1000 "\t" l[$1] } }' ${med_len} - > ${count}
