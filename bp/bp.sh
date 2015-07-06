
GENOME=hg19
CPU=16
GTF=/data/fuy2/shared/${GENOME}/${GENOME}.gencode.gtf
GENOME_FA=/data/fuy2/shared/${GENOME}/${GENOME}.fa
PIPELINE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

OUTPUT_DIR=test_output
COMMON_DIR=${PIPELINE_DIR}/common/${GENOME}
mkdir -p ${COMMON_DIR}

# Build STAR index
if [ ! -f ${COMMON_DIR}/.starindex.ok ]; then
    mkdir -p ${COMMON_DIR}/STARIndex && cd ${COMMON_DIR}/STARIndex
    echo "Building STAR index..."
    STAR --runMode genomeGenerate --runThreadN $CPU --genomeDir . --genomeFastaFiles ${GENOME_FA} --sjdbGTFfile ${GTF} --sjdbOverhang 100 &> STAR.build_log && touch ../.starindex.ok
else
    echo "STAR index already exists. If you want to rebuild it, please delete .starindex.ok and STARIndex."
fi

cd ${COMMON_DIR}
INTRON_FIVE_PRIME=23
if [ ! -f .intronprefix.ok ]; then
    # gtf_extract_intron_to_gtf.sh can be found at https://github.com/yfu/tools/blob/master/gtf_extract_intron_to_gtf.sh
    cat ${GTF} | gtf_extract_intron_to_gtf.sh > ${GENOME}.gencode.intron.gtf

    # Build the index for unique introns (first 23nt)
    # If multiple introns have the same 5' ends, keep just one of them
    awk '{ k=$4 $7; if(k in d) {} else { d[k]=$0 } } END{ for(k in d) { print d[k] } }' ${GENOME}.gencode.intron.gtf | sort -k1,1 -k4,4n > ${GENOME}.gencode.intron.uniq.gtf
    cat ${GENOME}.gencode.intron.uniq.gtf | tr -d '";' | awk '{ print $1 "\t" $4-1 "\t" $5 "\t" $12 "." $14 "\t" 0 "\t" $7 }' | awk -v fl=${INTRON_FIVE_PRIME} 'BEGIN{OFS="\t"} { if($3-$2>=fl) { $3=$2+fl; print } }' > ${GENOME}.gencode.intron.uniq.${INTRON_FIVE_PRIME}nt.bed
    bedtools getfasta -s -fi ${GENOME_FA} -bed ${GENOME}.gencode.intron.uniq.${INTRON_FIVE_PRIME}nt.bed -fo - > ${GENOME}.gencode.intron.uniq.${INTRON_FIVE_PRIME}nt.fa
    bowtie2-build ${GENOME}.gencode.intron.uniq.${INTRON_FIVE_PRIME}nt.fa ${GENOME}.gencode.intron.uniq.${INTRON_FIVE_PRIME}nt

    # Build the index for unique introns (not just the frist 23nt)
    cat ${GENOME}.gencode.intron.uniq.gtf | tr -d '";' | awk '{ print $1 "\t" $4-1 "\t" $5 "\t" $12 "." $14 "\t" 0 "\t" $7 }' | awk -v fl=${INTRON_FIVE_PRIME} 'BEGIN{OFS="\t"} { if($3-$2>=fl) { print } }' > ${GENOME}.gencode.intron.uniq.bed
    bedtools getfasta -s -fi ${GENOME_FA} -bed ${GENOME}.gencode.intron.uniq.bed -fo - > ${GENOME}.gencode.intron.uniq.fa
    bowtie2-build ${GENOME}.gencode.intron.uniq.fa ${GENOME}.gencode.intron.uniq
    touch .intronprefix.ok
else
    echo "Intron prefix has already been generated."
fi

FQ=$1
# for FQ in Mattick.CSQ.K562.BP.rep{1,2,3}.fq; do

OUTPUT_PREFIX=${FQ%.fq}.${GENOME}.
STAR \
    --runMode alignReads \
    --genomeDir ${COMMON_DIR}/STARIndex \
    --readFilesIn ${FQ}\
    --runThreadN $CPU \
    --outFilterScoreMin 0 \
    --outFilterScoreMinOverLread 0.72 \
    --outFilterMatchNmin 0 \
    --outFilterMatchNminOverLread 0.72 \
    --outFilterMultimapScoreRange 1 \
    --outFilterMultimapNmax -1 \
    --outFilterMismatchNmax 10 \
    --outFilterMismatchNoverLmax 0.05 \
    --alignIntronMax 0 \
    --alignIntronMin 21 \
    --outFilterIntronMotifs None \
    --genomeLoad NoSharedMemory \
    --outFileNamePrefix $OUTPUT_PREFIX \
    --outSAMunmapped None \
    --outReadsUnmapped Fastx \
    --outSJfilterReads Unique \
    --seedSearchStartLmax 20 \
    --seedSearchStartLmaxOverLread 1.0 \
    --chimSegmentMin 0 2>&1 1> ${OUTPUT_PREFIX}STAR.log
    
    
