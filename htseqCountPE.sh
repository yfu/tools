#! /bin/bash -x
case ${3} in
## mouse specific intersectBed target files
"mouse")
## TODO: change those to BED
	FOLDER=/home/hanb/nearline/small_RNA_Pipeline/common_files/mm9/UCSC_BEDS
	# piRNA clusters annoated in Li, et al, Mol Cell, 2013
	MOUSE_PREPACHYTENE_PIRNA_CLUSTER=$FOLDER/piRNA.cluster.prepachytene.bed
	MOUSE_HYBRID_PIRNA_CLUSTER=$FOLDER/piRNA.cluster.hybrid.bed
	MOUSE_PACHYTENE_PIRNA_CLUSTER=$FOLDER/piRNA.cluster.pachytene.bed
	# NM transcripts annotated in Zamore lab, 2013
	MOUSE_ZAMORE_NM=$FOLDER/x100910.rnaseq.transcripts.NM.all.final.bed.c
	# NR transcripts annotated in Zamore lab, 2013
	MOUSE_ZAMORE_NR=$FOLDER/x100910.rnaseq.transcripts.NR.all.final.bed.c
	# transcriptome combined from piRNA cluster, NM and NR
	MOUSE_ZAMORE_TRANSCRIPTOME=$FOLDER/zamore.transcriptome.bed
	# take the col1-6 of $MOUSE_ZAMORE_TRANSCRIPTOME (including intron)
	MOUSE_ZAMORE_ANNOTATION=$FOLDER/zamore.annotation.bed 
	
	MOUSE_refSeq_GENE=$FOLDER/UCSC.refSeq.Genes.bed
	MOUSE_refSeq_EXON=$FOLDER/UCSC.refSeq.Exon.bed
	MOUSE_refSeq_INTRON=$FOLDER/UCSC.refSeq.Intron.bed
	MOUSE_refSeq_5UTR=$FOLDER/UCSC.refSeq.5UTR.bed
	MOUSE_refSeq_3UTR=$FOLDER/UCSC.refSeq.3UTR.bed
	
	MOUSE_REPEATMASKER_DNA=$FOLDER/repeat_mask.bed.DNA
	MOUSE_REPEATMASKER_SINE=$FOLDER/repeat_mask.bed.SINE
	MOUSE_REPEATMASKER_LINE=$FOLDER/repeat_mask.bed.LINE
	MOUSE_REPEATMASKER_LTR=$FOLDER/repeat_mask.bed.LTR
	MOUSE_REPEATMASKER_SATELLITE=$FOLDER/repeat_mask.bed.Satellite
	MOUSE_REPEATMASKER_SIMPLEREPEATS=$FOLDER/repeat_mask.bed.Simple_repeat
	MOUSE_REPEATMASKER_tRNA=$FOLDER/repeat_mask.bed.tRNA
	
	declare -a TARGETS=( \
	"MOUSE_PIRNACLUSTER" \
    "MOUSE_PREPACHYTENE_PIRNA_CLUSTER" \
    "MOUSE_HYBRID_PIRNA_CLUSTER" \
    "MOUSE_PACHYTENE_PIRNA_CLUSTER" \
    "MOUSE_ZAMORE_NM" \
    "MOUSE_ZAMORE_NR" \
    "MOUSE_ZAMORE_TRANSCRIPTOME" \
    "MOUSE_refSeq_GENE" \
    "MOUSE_refSeq_EXON" \
    "MOUSE_refSeq_INTRON" \
    "MOUSE_refSeq_5UTR" \
    "MOUSE_refSeq_3UTR" \
    "MOUSE_REPEATMASKER_DNA" \
    "MOUSE_REPEATMASKER_SINE" \
    "MOUSE_REPEATMASKER_LINE" \
    "MOUSE_REPEATMASKER_LTR" \
    "MOUSE_REPEATMASKER_SATELLITE" \
    "MOUSE_REPEATMASKER_SIMPLEREPEATS" \
    "MOUSE_REPEATMASKER_tRNA" )
;;
## fly specific intersectBed target files
"fly")
	FOLDER=/data/fuy2/shared/annotations
	
	FLY_PIRNA_CLUSTER=$FOLDER/Brennecke.pirnaCluster.gtf 
	FLY_flyBase_GENE=$FOLDER/UCSC.flyBase.xnostrand.Genes.gtf
	FLY_flyBase_EXON=$FOLDER/UCSC.flyBase.xnostrand.Exon.gtf
	FLY_flyBase_INTRON=$FOLDER/UCSC.flyBase.xnostrand.Intron.gtf
	FLY_flyBase_5UTR=$FOLDER/UCSC.flyBase.xnostrand.5UTR.gtf
	FLY_flyBase_3UTR=$FOLDER/UCSC.flyBase.xnostrand.3UTR.gtf
	FLY_REPEATMASKER_DNA=$FOLDER/repeat_mask.gtf.DNA
	FLY_REPEATMASKER_LINE=$FOLDER/repeat_mask.gtf.LINE
	FLY_REPEATMASKER_LTR=$FOLDER/repeat_mask.gtf.LTR
	FLY_REPEATMASKER_RNA=$FOLDER/repeat_mask.gtf.RNA
	FLY_REPEATMASKER_rRNA=$FOLDER/repeat_mask.gtf.rRNA
	FLY_REPEATMASKER_Satellite=$FOLDER/repeat_mask.gtf.Satellite
	FLY_REPEATMASKER_Simple_repeat=$FOLDER/repeat_mask.gtf.Simple_repeat
	FLY_REPEATMASKER_Unknown=$FOLDER/repeat_mask.gtf.Unknown
	
	FLY_TRANSPOSON_ALL=$FOLDER/transposon.gtf
	FLY_TRANSPOSON_GROUP1=$FOLDER/Zamore.group1.gtf
	FLY_TRANSPOSON_GROUP2=$FOLDER/Zamore.group2.gtf
	FLY_TRANSPOSON_GROUP3=$FOLDER/Zamore.group3.gtf
	FLY_TRANSPOSON_GROUP0=$FOLDER/Zamore.group0.gtf
	
	FLY_TRANSPOSON_OUTCLUSTER=$FOLDER/transposon.outCluster.gtf
	FLY_TRANSPOSON_INCLUSTER=$FOLDER/transposon.inCluster.gtf
	
	FLY_TRANSPOSON_OUTCLUSTER_INGENE=$FOLDER/transposon.outCluster.inGene.gtf
	FLY_TRANSPOSON_OUTCLUSTER_OUTGENE=$FOLDER/transposon.outCluster.outGene.gtf
	
	FLY_TRANSPOSON_INCLUSTER_INGENE=$FOLDER/transposon.inCluster.inGene.gtf
	FLY_TRANSPOSON_INCLUSTER_OUTGENE=$FOLDER/transposon.inCluster.outGene.gtf
	
	FLY_GENE_TRANSPOSON_ALL=$FOLDER/UCSC.flyBase.xnostrand.Genes.transposon.gtf #for statistic purpose
	
	declare -a TARGETS=("FLY_PIRNA_CLUSTER" \
	"FLY_TRANSPOSON_ALL" \
	"FLY_TRANSPOSON_GROUP1" \
	"FLY_TRANSPOSON_GROUP2" \
	"FLY_TRANSPOSON_GROUP3" \
	"FLY_TRANSPOSON_GROUP0" \
	"FLY_TRANSPOSON_OUTCLUSTER" \
	"FLY_TRANSPOSON_INCLUSTER" \ 
	"FLY_TRANSPOSON_OUTCLUSTER_INGENE" \
	"FLY_TRANSPOSON_OUTCLUSTER_OUTGENE" \
	"FLY_TRANSPOSON_INCLUSTER_INGENE" \
	"FLY_TRANSPOSON_INCLUSTER_OUTGENE" \
     "FLY_flyBase_GENE" \
     "FLY_flyBase_EXON" \
     "FLY_flyBase_INTRON" \
     "FLY_flyBase_5UTR" \
     "FLY_flyBase_3UTR" \
     "FLY_REPEATMASKER_DNA" \
     "FLY_REPEATMASKER_LINE" \
     "FLY_REPEATMASKER_LTR" \
     "FLY_REPEATMASKER_Satellite" \
     "FLY_REPEATMASKER_Simple_repeat" \
     "FLY_REPEATMASKER_RNA" \
     "FLY_REPEATMASKER_rRNA" \
    "FLY_REPEATMASKER_Unknown" \
    "FLY_GENE_TRANSPOSON_ALL" )
;;
*)
	echo "unknown orgnanism... currently only mouse/fly is supported"
	exit 2
;;
esac

SAM=$1
SAM_NAME=${1##*/}
LIBTYPE=$2

# -s tag of htseq-count:
## yes: \1 is on the same strand
## reverse: \2 is on the same strand 
## no: unstranded

		#F) CAGE
			#	STRAND=$((STRAND+2))
		#;;
		#R) RNASeq
			#	STRAND=$((STRAND+1))

case $LIBTYPE in	
1) # RNASeq/ChIPSeq
	SENSE_STRAND="reverse"
	ANTI_SENSE_STRAND="yes"
;;
2) # CAGE/Degradome 
	SENSE_STRAND="yes" # for CAGE/DEG, \1 is on the same strand, thus "yes" tag
	ANTI_SENSE_STRAND="reverse" # for CAGE/DEG, \2 is on the reverse strand, thus "reverse" tag
;;
esac

COUNTOUTDIR=$4
CPU=$5

#${INTERSECTOUTDIR}/${UNIQUE_BED_NAME%bam}
#${INTERSECTOUTDIR}/${MULTIP_BED_NAME%bam}
# @ running pairtobed 
paraFile=${COUNTOUTDIR}/htseqCount.${RANDOM}.para && \
rm -rf $paraFile
for t in ${TARGETS[@]}

do \
	echo "[ ! -s ${COUNTOUTDIR}/${SAM_NAME%sam}${t}.htseqcount.S.out ] && htseq-count -m intersection-strict -s $SENSE_STRAND -t exon -i transcript_id -q $SAM ${!t} > ${COUNTOUTDIR}/${SAM_NAME%sam}${t}.htseqcount.S.out "	>>  $paraFile
	echo "[ ! -s ${COUNTOUTDIR}/${SAM_NAME%sam}${t}.htseqcount.AS.out ] && htseq-count -m intersection-strict -s $ANTI_SENSE_STRAND -t exon -i transcript_id -q $SAM ${!t} > ${COUNTOUTDIR}/${SAM_NAME%sam}${t}.htseqcount.AS.out "	>>  $paraFile
done

[ ! -f ${paraFile}.completed ] && \
	ParaFly -c $paraFile -CPU $CPU -failed_cmds $paraFile.failed_commands

