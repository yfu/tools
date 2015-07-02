
GENOME=hg19
CPU=16
ROOT=/data/fuy2/shared/${GENOME}
GTF=${ROOT}/${GENOME}.gencode.gtf
GENOME_FA=${ROOT}/hg19.fa
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
if [ ! -f .intronprefix.ok ]; then
    # gtf_extract_intron_to_gtf.sh can be found at https://github.com/yfu/tools/blob/master/gtf_extract_intron_to_gtf.sh
    cat ${GTF} | gtf_extract_intron_to_gtf.sh > ${GENOME}.gencode.intron.gtf
    # If multiple introns have the same 5' ends, keep just one of them
    awk '{ k=$4 $7; if(k in d) {} else { d[k]=$0 } } END{ for(k in d) { print d[k] } }' ${GENOME}.gencode.intron.gtf | sort -k1,1 -k4,4n > ${GENOME}.gencode.intron.uniq.gtf 
    bedtools getfasta -tab -s -fi ${GENOME_FA} -bed ${GENOME}.gencode.intron.gtf > ${GENOME}.gencode.intron.fa
    touch .intronprefix.ok
else
    echo "Intron prefix has already been generated."
fi
    
INTRON_PREFIX_L=23

