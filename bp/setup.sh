
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
    # If multiple introns have the same 5' end, 3' end and strand, keep just one of them
    awk '{ k=$1 $4 $5 $7; if(k in d) {} else { d[k]=$0 } } END{ for(k in d) { print d[k] } }' ${GENOME}.gencode.intron.gtf | sort -k1,1 -k4,4n > ${GENOME}.gencode.intron.uniq.gtf
    cat ${GENOME}.gencode.intron.uniq.gtf | tr -d '";' | awk '{ print $1 "\t" $4-1 "\t" $5 "\t" $12 "." $14 "\t" 0 "\t" $7 }' | awk -v fl=${INTRON_FIVE_PRIME} 'BEGIN{OFS="\t"} { if($3-$2>=fl) { if($6=="+") {$3=$2+fl; print } else {$2=$3-fl; print} } }' | awk '{ if($6=="+") {watson[$1 $2]=$0} else { crick[$1 $2]=$0 } } END{ for (w in watson){print watson[w]} for (c in crick) {print crick[c]}  }' > ${GENOME}.gencode.intron.uniq.${INTRON_FIVE_PRIME}nt.bed
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

