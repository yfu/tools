#!/usr/bin/env bash
while getopts "hf:i:u:a:g:" OPTION; do
    case $OPTION in
	h)usage && exit 0 ;;
	i)raw_count=`readlink -f ${OPTARG}`;; # the directory of pipeline results
	u)nf_uniq="${OPTARG}";; # normalization factor based on # unique mappers
	a)nf_all="${OPTARG}";; # nf based on # all mappers	
	g)GENOME=$OPTARG ;;
    esac
done

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

count=${raw_count}.naive_counting
grep -v '^__' ${raw_count} | awk -v nf_uniq=${nf_uniq} -v nf_all=${nf_all} 'BEGIN{print "# gene\traw\tppm_uniq\tppm_all\tfpkm_uniq\tfpkm_all\tmed_len"} { if(FNR==NR) { l[$1]=$2 } else { print $1 "\t" $2 "\t" $2*nf_uniq "\t" $2*nf_all "\t" $2*nf_uniq/l[$1]*1000 "\t" $2*nf_all/l[$1]*1000 "\t" l[$1] } }' ${med_len} - > ${count}
