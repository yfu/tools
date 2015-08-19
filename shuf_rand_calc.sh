#!/usr/bin/env bash

# Similarity score: chop the given seq into 10-mers and get the frequency of each 10-mer according to the JellyFish database; then the sum of them is the S-Score
# Given a JellyFish database and a single-entry fasta file, calculate the S-scores of the input seq, shuffled input seq, random seq of the same size and shuffled random seq

iter=50
jfdb=$1
seq=$2
prefix=${seq%.fa}
prefix=${prefix%.fasta}
scores=${prefix}.scores
chrominfo=~/data/shared/mm10/mm10.chrom.sizes
genome_fa=~/data/shared/mm10/mm10.fa
#### Upregulated genes

echo "# Part 1: Score for the input sequence"
jellyfish query ${jfdb} -s ${seq} | awk '{ s+=$2 } END{ print "uniq_sequence\t" s }' > ${scores}

echo "# Part 2: Shuffled input sequence"
# Get the # bases
l=$(fasta_lengths.pl ${seq} | awk '{ s+=$2 } END{ print s }')
cmd_1=${prefix}.step1.${RANDOM}${RANDOM}.cmd
cmd_2=${prefix}.step2.${RANDOM}${RANDOM}.cmd
for i in $(seq 1 ${iter}); do
    echo "fasta-dinucleotide-shuffle.py -f ${seq} -s ${RANDOM} 2>/dev/null > ${prefix}.shuffled.${i}.tmp.fa" >> ${cmd_1}
    echo "jellyfish query ${jfdb} -s ${prefix}.shuffled.${i}.tmp.fa | awk '{ s+=\$2 } END{ print s }' > ${prefix}.shuffled.${i}.tmp.score" >> ${cmd_2}
done
cat ${cmd_1} | parallel --progress -j 64
cat ${cmd_2} | parallel --progress -j 64
for i in $(seq 1 ${iter}); do
    cat ${prefix}.shuffled.${i}.tmp.score | awk -v i=${i} '{ print "shuffled_" i "\t" $1 }' >> ${scores}
done

echo "# Part 3: Random regions of the same size"
for i in $(seq 1 ${iter}); do
    bedtools random -g ${chrominfo} -l ${l} -n 1 > ${prefix}.random.${i}.tmp.bed
    bedtools getfasta -fi ${genome_fa} -bed ${prefix}.random.${i}.tmp.bed -fo - > ${prefix}.random.${i}.tmp.fa
    jellyfish query ${jfdb} -s ${prefix}.random.${i}.tmp.fa | awk -v i=${i} '{ s+=$2 } END{ print "random_" i "\t" s }' >> ${scores}
done

echo "# Part 4: Shuffled random"
cmd_3=${prefix}.step3.${RANDOM}${RANDOM}.cmd
cmd_4=${prefix}.step4.${RANDOM}${RANDOM}.cmd
for idx in $(seq 1 5); do
    for i in $(seq 1 ${iter}); do
	echo "fasta-dinucleotide-shuffle.py -f ${prefix}.random.${idx}.tmp.fa -s ${RANDOM} 2>/dev/null > ${prefix}.random.${idx}.shuffled.${i}.tmp.fa" >> ${cmd_3}
	echo "jellyfish query ${jfdb} -s ${prefix}.random.${idx}.shuffled.${i}.tmp.fa | awk '{ s+=\$2 } END{ print s }' > ${prefix}.random.${idx}.shuffled.${i}.tmp.score" >> ${cmd_4}
    done
done
cat ${cmd_3} | parallel --progress -j 64
cat ${cmd_4} | parallel --progress -j 64
for idx in $(seq 1 5); do
    for i in $(seq 1 ${iter}); do
	cat ${prefix}.random.${idx}.shuffled.${i}.tmp.score | awk -v idx=${idx} -v i=${i} '{ print "random_" idx "_shuffled_" i "\t" $1 }' >> ${scores}
    done
done

cat ${scores} | awk '{ if($2 ~ /[0-9]+/) {print} }' > ${scores}.valid

# rm *tmp*
