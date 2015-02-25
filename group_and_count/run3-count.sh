# This script will count given one output from run1-group-bed2.sh

rmsk_dir=/home/fuy2/data/shared/mm10
rmsk=${rmsk_dir}/UCSC.rmsk.sorted.bed

rmsk_dna=${rmsk_dir}/UCSC.rmsk.DNA.sorted.bed
rmsk_line=${rmsk_dir}/UCSC.rmsk.LINE.sorted.bed
rmsk_ltr=${rmsk_dir}/UCSC.rmsk.LTR.sorted.bed
rmsk_rrna=${rmsk_dir}/UCSC.rmsk.rRNA.sorted.bed
rmsk_satellite=${rmsk_dir}/UCSC.rmsk.Satellite.sorted.bed
rmsk_simplerepeats=${rmsk_dir}/UCSC.rmsk.Simple_repeat.sorted.bed
rmsk_sine=${rmsk_dir}/UCSC.rmsk.SINE.sorted.bed
rmsk_trna=${rmsk_dir}/UCSC.rmsk.tRNA.sorted.bed

cluster=piRNA.cluster.bed12
cluster_exon=piRNA.cluster.bed6
prepachytene=piRNA.cluster.prepachytene.bed12
prepachytene_exon=piRNA.cluster.prepachytene.bed6
hybrid=piRNA.cluster.hybrid.bed12
hybrid_exon=piRNA.cluster.hybrid.bed6
pachytene=piRNA.cluster.pachytene.bed12
pachytene_exon=piRNA.cluster.pachytene.bed6

intersect() {
    bed2=$1
    anno=$2
    output_prefix=$3
    echo "bedtools intersect -sorted -a ${bed2} -b ${anno} -u -s | awk '{ s+=\$4/\$5 } END{ print s }' > ${output_prefix}.s.count"
    echo "bedtools intersect -sorted -a ${bed2} -b ${anno} -u -S | awk '{ s+=\$4/\$5 } END{ print s }' > ${output_prefix}.as.count"
}

prefix=$1
prefix=02dpp
bed2=${prefix}.rmsk.bed2
rep=$rmsk

intersect $bed2 $rep ${bed2%.bed2}

for i in dna line ltr rrna satellite simplerepeats sine trna; do
    # different repeats usually do not overlap, simply do the following
    eval "rep=\${rmsk_${i}}"
    eval "bed2=${prefix}.rmsk_${i}.bed2"
    intersect ${bed2} $rep ${bed2%.bed2}
done

for i in protein_coding lincRNA pseudogene; do
    eval ${i}=gene.${i}.bed
    eval ${i}_exon=exon.${i}.bed
    eval ${i}_upstream5k=gene.${i}.upstream5k.bed
    eval ${i}_downstream5k=gene.${i}.downstream5k.bed    
done

targets=(cluster lincRNA pseudogene protein_coding protein_coding_upstream5k protein_coding_downstream5k lincRNA_upstream5k lincRNA_downstream5k pseudogene_upstream5k pseudogene_downstream5k)
# -rw-r--r--  1 fuy2 UNIXLDAPUsers  524 Feb 25 09:51 02dpp.cluster.bed2
# -rw-r--r--  1 fuy2 UNIXLDAPUsers  169 Feb 25 09:51 02dpp.lincRNA.bed2
# -rw-r--r--  1 fuy2 UNIXLDAPUsers    0 Feb 25 09:51 02dpp.lincRNA_downstream5k.bed2
# -rw-r--r--  1 fuy2 UNIXLDAPUsers    0 Feb 25 09:51 02dpp.lincRNA_upstream5k.bed2
# -rw-r--r--  1 fuy2 UNIXLDAPUsers 2.9K Feb 25 09:51 02dpp.other.bed2
# -rw-r--r--  1 fuy2 UNIXLDAPUsers 1.9K Feb 25 09:51 02dpp.protein_coding.bed2
# -rw-r--r--  1 fuy2 UNIXLDAPUsers  636 Feb 25 09:51 02dpp.protein_coding_downstream5k.bed2
# -rw-r--r--  1 fuy2 UNIXLDAPUsers  527 Feb 25 09:51 02dpp.protein_coding_upstream5k.bed2
# -rw-r--r--  1 fuy2 UNIXLDAPUsers    0 Feb 25 09:51 02dpp.pseudogene.bed2
# -rw-r--r--  1 fuy2 UNIXLDAPUsers    0 Feb 25 09:51 02dpp.pseudogene_downstream5k.bed2
# -rw-r--r--  1 fuy2 UNIXLDAPUsers    0 Feb 25 09:51 02dpp.pseudogene_upstream5k.bed2
for i in ${targets[@]}; do
    eval t=\$$i
    echo "Processing target: $t" >&2
    bed2=${prefix}.$i.bed2
    intersect ${bed2} ${t} ${bed2%.bed2}
done

# For those features that have exons:
targets=(cluster lincRNA pseudogene protein_coding)
for i in ${targets[@]}; do
    echo "Processing exon target: $t" >&2
    eval t=\$${i}_exon
    intersect ${bed2} ${t} ${bed2%.bed2}.exon
done

# For the unnannoted, I do not know if it is sense or antisnese...

awk '{ s+=$4/$5 } END{ print s }' ${prefix}.other.bed2 > ${prefix}.other.count
