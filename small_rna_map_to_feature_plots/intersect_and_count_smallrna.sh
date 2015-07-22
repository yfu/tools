#!/usr/bin/env bash

source ~/repo/tools/piPipes_aux/global.sh
source ~/repo/tools/piPipes_aux/functions.sh

# Usage:
# bash run1.sh -N allxmirna -a aaa -o aaaout -g dm3_chengjian -A aaa_name
# intersect_and_count_smallrna.sh -N mirna -g mm10g -a Zamore.SRA.pi7M_D3.59dpp -A Zamore.SRA.pi7M_D3.59dpp.name -o test

lib_path=$HOME/repo/tools/piPipes_aux

source ${lib_path}/global.sh
source ${lib_path}/functions.sh

#############################
# ARGS reading and checking #
#############################
while getopts "hva:b:g:c:o:A:B:N:c" OPTION; do
    case $OPTION in
	h)usage && exit 0 ;;
	v)echo2 "SMALLRNA2_VERSION: v$SMALLRNA2_VERSION" && exit 0 ;;
	a)SAMPLE_A_DIR=`readlink -f $OPTARG` ;;
	# b)SAMPLE_B_DIR=`readlink -f $OPTARG` ;;
	o)PARA_OUT=`readlink -f $OPTARG` ;;
	g)export GENOME=${OPTARG};;
	A)export SAMPLE_A_NAME=$OPTARG ;;
	# B)  export SAMPLE_B_NAME=$OPTARG ;;
	N)export NORMMETHOD=`echo ${OPTARG} | tr '[A-Z]' '[a-z]'` ;;
	c)CPU=$OPTARG ;;
	*)usage && exit 1 ;;
    esac
done

# The script is meant to be run in parallel so default CPU is 1
if [ -z "$CPU" ]; then
    CPU=1
fi
source ${lib_path}/${GENOME}.parameters.sh

# ALL_BED2_A=`ls $SAMPLE_A_DIR/intersect_genomic_features/*${GENOME}*all.x_rpmk_MASK.bed2`
ALL_BED2_A=${SAMPLE_A_DIR}/genome_mapping/*.all.piRNA.bed2
# Use the bed2 without Hsp70Bb_in_UASt_248
# ALL_BED2_A=${SAMPLE_A_DIR}/genome_mapping/*.no_hsp5primeutr.bed2

echo "***"
echo "Using this file as input: ${ALL_BED2_A}"
echo "***"
# NORMMETHOD=allxmirna
echo "Normalization method" ${NORMMETHOD}

case "$NORMMETHOD" in
    unique)
	SAMPLE_A_NORMFACTOR=`head -6 $SAMPLE_A_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	# SAMPLE_B_NORMFACTOR=`head -6 $SAMPLE_B_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	;;
    uniquexmirna)
	SAMPLE_A_NORMFACTOR=`head -7 $SAMPLE_A_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	# SAMPLE_B_NORMFACTOR=`head -7 $SAMPLE_B_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	;;
    all)
	SAMPLE_A_NORMFACTOR=`head -4 $SAMPLE_A_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	# SAMPLE_B_NORMFACTOR=`head -4 $SAMPLE_B_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	;;
    allxmirna)
	SAMPLE_A_NORMFACTOR=`head -5 $SAMPLE_A_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	# SAMPLE_B_NORMFACTOR=`head -5 $SAMPLE_B_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	;;
    mirna)
	SAMPLE_A_NORMFACTOR=`head -3 $SAMPLE_A_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	# SAMPLE_B_NORMFACTOR=`head -3 $SAMPLE_B_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	;;
    sirna)
	case $GENOME in
	    dm3)
		SAMPLE_A_NORMFACTOR=`grep -P "structural_loci|cisNATs" $SAMPLE_A_DIR/summaries/*siRNA.sum | awk '{a+=$9}END{print 1000000/a}'`
		# SAMPLE_B_NORMFACTOR=`grep -P "structural_loci|cisNATs" $SAMPLE_B_DIR/summaries/*siRNA.sum | awk '{a+=$9}END{print 1000000/a}'`
		;;
	    *)
		echo2 "The annotation for siRNA in ${GENOME} is poor. Please choose a different normalization method. \nIf unox library, choose \"miRNA\". If ox, choose \"uniquexmirna\"" "error"
		;;
	esac
	;;
    42ab)
	case $GENOME in
	    dm3)
		SAMPLE_A_NORMFACTOR=`grep -P "piRNA_Cluster_42AB" $SAMPLE_A_DIR/summaries/*piRNA.sum | awk '{a+=$9}END{print 1000000/a}'`
		# SAMPLE_B_NORMFACTOR=`grep -P "piRNA_Cluster_42AB" $SAMPLE_B_DIR/summaries/*piRNA.sum | awk '{a+=$9}END{print 1000000/a}'`
		;;
	    *)
		echo2 "this normalization method is not supported for this organism" "error"
		;;
	esac
	;;
    flam)
	case $GENOME in
	    dm3)
		SAMPLE_A_NORMFACTOR=`grep -P "piRNA_Cluster_flam" $SAMPLE_A_DIR/summaries/*piRNA.sum | awk '{a+=$9}END{print 1000000/a}'`
		# SAMPLE_B_NORMFACTOR=`grep -P "piRNA_Cluster_flam" $SAMPLE_B_DIR/summaries/*piRNA.sum | awk '{a+=$9}END{print 1000000/a}'`
		;;
	    *)
		echo2 "this normalization method is not supported for this organism" "error"
		;;
	esac
	;;
    *)
	echo2 "unrecognized normalization option: $NORMMETHOD; using the default method" "warning"
	SAMPLE_A_NORMFACTOR=`head -6 $SAMPLE_A_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	# SAMPLE_B_NORMFACTOR=`head -6 $SAMPLE_B_DIR/*basic_stats | tail -1 | cut -f2 | awk '{print 1000000/$0}'`
	;;
esac

COMMON_FOLDER=/home/fuy2/data/piPipes/common/${GENOME}
# I have to specify it this way since the genome names are different
if [[ ${GENOME} == 'dm3' ]]; then
    Trn=~/data/piPipes/common/dm3/Zamore.transposon.bed.gz
    piRNA_Cluster=~/data/piPipes/common/dm3/Brennecke.piRNAcluster.bed6.gz
elif [[ ${GENOME} == 'mm10g' ]]; then
    echo "Nothing special about ${GENOME}"
    # source ~/repo/tools/piPipes_aux/${GENOME}.parameters.sh
else
    echo "Unrecognized genome. Exiting..."
    exit 1
fi

if [[ -f ${PARA_OUT} ]]; then
    echo "${PARA_OUT} exists. I will delete it."
    rm ${PARA_OUT}
fi
if [[ ${GENOME} == 'dm3' ]]; then
    #########################################
    # draw piRNA abundance for genic mapper #
    #########################################
    echo2 "Calculating piRNA abundances for each genic transcripts"
    threePrimeUTR=$COMMON_FOLDER/UCSC.refSeq.3UTR.bed.gz
    echo "awk '\$3-\$2>=$piRNA_bot' $ALL_BED2_A | ${BEDTOOLS} intersect -split -wo -f 0.99 -a stdin -b $threePrimeUTR | awk -v depth=$SAMPLE_A_NORMFACTOR '{if (\$6==\$13) a[\$11]+=\$4/\$NF/\$5; else b[\$11]+=\$4/\$NF/\$5; total[\$11]=1}END{for (c in total) { printf \"%s\t%.2f\t%.2f\n\", c, (a[c]?a[c]:0)*depth, (b[c]?b[c]:0)*depth}}' > ${SAMPLE_A_NAME}.genic.abundance.normalized_by_$NORMMETHOD " >> ${PARA_OUT}
    STEP=$((STEP+1))
    
    ###################################################
    # draw piRNA abundance for each transposon family #
    ###################################################
    echo2 "Calculating piRNA abundances for each transposon family"
    echo "awk '\$3-\$2>=$piRNA_bot' $ALL_BED2_A | ${BEDTOOLS} intersect -split -wo -f 0.99 -a stdin -b $Trn | awk -v depth=$SAMPLE_A_NORMFACTOR '{split(\$11,name,\".\"); if (\$6==\$13) {a[name[1]]+=\$4/\$NF/\$5; sa[name[1]]+=(\$4/\$NF/\$5)*(\$3-\$2);} else {b[name[1]]+=\$4/\$NF/\$5; sb[name[1]]+=(\$4/\$NF/\$5)*(\$3-\$2);}; total[name[1]]=1; class[name[1]]=\$12;}END{for (c in total) { printf \"%s\t%s\t%.2f\t%.2f\n\", c, class[c], (a[c]?a[c]:0)*depth, (b[c]?b[c]:0)*depth; printf \"%s\t%s\t%.2f\t%.2f\n\", c, class[c], a[c]==0?0:( (sa[c]?sa[c]:0)/a[c] ), b[c]==0?0: ( (sb[c]?sb[c]:0)/b[c] ) >> \"/dev/stderr\";}}' > ${SAMPLE_A_NAME}.transposon.abundance.normalized_by_$NORMMETHOD 2> ${SAMPLE_A_NAME}.transposon.mean_len.normalized_by_$NORMMETHOD" >> $PARA_OUT
    # awk '$3>0&&$4>0' ${SAMPLE_A_NAME}.transposon.mean_len.normalized_by_$NORMMETHOD > $TRN_DIR/${SAMPLE_A_NAME}.transposon.mean_len.normalized_by_${NORMMETHOD}.no_zero
    STEP=$((STEP+1))

    ###############################################
    # draw piRNA abundance for each piRNA cluster #
    ###############################################
    echo2 "Calculating piRNA abundances for each piRNA cluster family"
    echo "awk '\$3-\$2>=$piRNA_bot' $ALL_BED2_A | ${BEDTOOLS} intersect -split -wo -f 0.99 -a stdin -b $piRNA_Cluster | awk -v depth=$SAMPLE_A_NORMFACTOR '{if (\$6==\$13) a[\$11]+=\$4/\$NF/\$5; else b[\$11]+=\$4/\$NF/\$5; total[\$11]=1}END{for (c in total) { printf \"%s\t%.2f\t%.2f\n\", c, (a[c]?a[c]:0)*depth, (b[c]?b[c]:0)*depth}}' > ${SAMPLE_A_NAME}.cluster.abundance.normalized_by_$NORMMETHOD " >> $PARA_OUT
    echo2 "Done with small RNA pipeline dual-sample mode"

    ###########################################
    # draw piRNA abundance for each construct #
    ###########################################
    # echo2 "Calculating piRNA abundances for each construct"
    # echo "grep -E 'GSV6|nos-gal4-vp16|UASp-EGFP' ${ALL_BED2_A} | awk -v depth=depth '{if(\$3-\$2>=23) { total[\$1] = 1; if(\$6==\"+\") {a[\$1] += \$4/\$5} else {b[\$1] += \$4/\$5} } } END{ for(i in total) { print i \"\t\" a[i]*depth \"\t\" b[i]*depth } }' > ${SAMPLE_A_NAME}.construct.abundance.normalized_by_$NORMMETHOD" >> ${PARA_OUT}
elif [[ ${GENOME} == 'mm10g' ]]; then
    ##############################################
    # calculate piRNA abundance for genic mapper #
    ##############################################
    gtf=~/data/shared/mm10/gencode.vM4.corrected_chrom_names.gtf
    bed2=${SAMPLE_A_DIR}/genome_mapping/*.all.piRNA.bed2
    echo2 "Using this file: $(ls ${bed2})"
    echo "count_reads_for_gene_ntm.py ${gtf} ${ALL_BED2_A} | awk -v depth=$SAMPLE_A_NORMFACTOR 'BEGIN{OFS=\"\t\"} { \$2=\$2*depth; \$3=\$3*depth; print }' > ${SAMPLE_A_NAME}.gene.abundance.normalized_by_$NORMMETHOD" >> ${PARA_OUT}

    ###############################################
    # draw piRNA abundance for each piRNA gene #
    ###############################################
    for i in Prepachytene Hybrid Pachytene; do
	target=piRNA_Genes_${i}
	echo "Current target:" $target
	echo "awk '\$3-\$2>=$piRNA_bot' $ALL_BED2_A | ${BEDTOOLS} intersect -split -wo -f 0.99 -a stdin -b ${!target} | awk -v depth=$SAMPLE_A_NORMFACTOR '{if (\$6==\$13) a[\$11]+=\$4/\$NF/\$5; else b[\$11]+=\$4/\$NF/\$5; total[\$11]=1}END{for (c in total) { printf \"%s\t%.2f\t%.2f\n\", c, (a[c]?a[c]:0)*depth, (b[c]?b[c]:0)*depth}}' | awk '{ match(\$1, /^(.+)\.[0-9]+\.[0-9]+\$/, ary); i=ary[1]; a[i] += \$2; b[i] += \$3; total[i]=1 } END{ for(i in total) { s=a[i]?a[i]:0; as=b[i]?b[i]:0; print i \"\t\" s \"\t\" as } }' | sort -k1,1 > ${SAMPLE_A_NAME}.piGenes.${i}.abundance.normalized_by_$NORMMETHOD " >> $PARA_OUT
    done

    ##################################################
    # transposon piRNA abundance for each piRNA gene #
    ##################################################
    
    # transposon_MM=2 # defined in /home/fuy2/repo/tools/piPipes_aux
    idx=/data/fuy2/shared/mm10/bowtie_index/mm10.repBase
    t="repBase"
    INSERT=${SAMPLE_A_DIR}/input_read_files/*.x_rRNA.x_hairpin.insert
    PREFIX=${SAMPLE_A_NAME}
    echo "
bowtie -r -v ${transposon_MM} -a --best --strata -p $CPU \
       -S \
       ${idx} \
       ${INSERT} 2> ${t}.log | \
    samtools view -uS -F0x4 - 2>/dev/null | \
    samtools sort -T ${RANDOM}${RANDOM} -O bam -@ $CPU - | \
    ~/data/piPipes/bin/bedtools_piPipes bamtobed -i - > ${PREFIX}.${t}.a${transposon_MM}.insert.bed && \
    ~/data/piPipes/bin/piPipes_insertBed_to_bed2 $x_rRNA_INSERT ${PREFIX}.${t}.a${transposon_MM}.insert.bed > ${PREFIX}.${t}.a${transposon_MM}.insert.bed2 && awk -v depth=${SAMPLE_A_NORMFACTOR} '{ total[\$1]=1; if(\$6==\"+\") { a[\$1] += \$4/\$5 } else { b[\$1] += \$4/\$5 } } END{ for(i in total) { print i \"\t\" a[i]*depth \"\t\" b[i]*depth } }' ${PREFIX}.${t}.a${transposon_MM}.insert.bed2 > ${SAMPLE_A_NAME}.transposon.abundance.normalized_by_$NORMMETHOD" >> ${PARA_OUT}
fi
    
