while getopts "f:g:c:h" opt; do
    case "$opt" in
	h|\?)
	    echo "Usage: TODO"
	    exit 0
	    ;;
	g)  genome=$OPTARG
	    ;;
	f)  f=$OPTARG
	    ;;
	c)  cpu=$OPTARG
	    ;;
    esac
done

if [[ -z "$cpu" ]]; then
    cpu=16
fi

if [[ -z "${f}" || -z "${genome}" ]]; then
    echo "No arguments provided! Exit!"
    exit 1
else
    echo "Input file: $f"
    echo "Genome: $genome"
    echo "CPU: ${cpu}"
fi

# f=$1
# input=test1000000.fq
input=$(basename ${f})
ln -s ${f} ${input}

CIRCPATH=/home/fuy2/repo/tools/circ/find_circ
PATH=${CIRCPATH}:${PATH}

if [[ $genome == "hg19" ]]; then
    bowtie2_index=~/data/shared/hg19/bowtie2index/hg19
    genome_info=~/data/shared/hg19/by_chrom
elif [[ $genome == "mm10" ]]; then
    bowtie2_index=~/data/shared/mm10/bowtie2index/mm10
    genome_info=~/data/shared/mm10/by_chrom
else
    echo "Unrecognized genome! Exit!"
    exit 1
fi

prefix=${input%.fq}
prefix=${prefix%.fastq}
gm_log=${prefix}.genome_mapping.log
# gm_sens_opt="-D 20 -R 3 -i C,1,0"
gm_sens_opt="-D 20 -R 4 -i C,2,0 -L 20"
# gm_sens_opt="--very-sensitive"
gm_bam=${prefix}.genome_mapping.bam
unmapped_bam=${prefix}.unmapped.bam
# Note that "-M option is now deprecated. It will be removed in subsequent versions. What used to be called -M mode is still the default mode, and -k and -a are still there alternatives to the default mode, but adjusting the -M setting is deprecated. Use the -D and -R options to adjust the effort expended to find valid alignments."
bowtie2 -p ${cpu} ${gm_sens_opt} --phred33 --mm --score-min=C,-15,0 -x ${bowtie2_index} -q -U ${input} 2> ${gm_log} | samtools view -hbuS - | samtools sort -@ ${cpu} -O bam -T ${RANDOM}${RANDOM} - > ${gm_bam}
samtools view -hf 4 ${gm_bam} | samtools view -Sb - > ${unmapped_bam}
anchors=${prefix}.anchors.fq.gz
unmapped2anchors.py ${unmapped_bam} | gzip > ${anchors}

d=${prefix}_circ
mkdir ${d}
anchor_sens_opt=${gm_sens_opt}

circ_prefix=${prefix}_circ_
circ_log=${d}/${prefix}.sites.log
circ_bed=${d}/${prefix}.sites.bed
circ_reads=${d}/${prefix}.sites.reads
bowtie2_circ_log=${prefix}_bowtie2_circ.log
bowtie2 -p ${cpu} --reorder --mm ${anchor_sens_opt} --score-min=C,-15,0 -q -x ${bowtie2_index} -U ${anchors} 2> ${bowtie2_circ_log} | find_circ.py -G ${genome_info} -p ${circ_prefix} -s ${circ_log} > ${circ_bed} 2> ${circ_reads}
