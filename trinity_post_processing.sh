TRINITY_HOME=~/src/trinityrnaseq_r20140717
## ${TRINITY_HOME}/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left ../Ponten.RSQ.testis_7e_1.fastq --right ../Ponten.RSQ.testis_7e_2.fastq --est_method eXpress --aln_method bowtie2 --SS_lib_type RF --thread_count 16 --max_ins_size 800 --output_dir . --trinity_mode &> align_and_estimate_abundance.log

${TRINITY_HOME}/util/misc/count_features_given_MIN_FPKM_threshold.pl results.xprs > cumul_counts.txt

awk '{ if($11>1 && NR>1) {print $2} }' results.xprs > fpkm1.namelist
cut -d' ' -f1 Trinity.fasta > Trinity.fasta.trimmed_header
seqtk subseq Trinity.fasta.trimmed_header fpkm1.namelist -l 60 > Trinity.fpkm1.trimmed_header.fasta
rm Trinity.fasta.trimmed_header

bowtie-build Trinity.fpkm1.trimmed_header.fasta Trinity.fpkm1.trimmed_header.fasta
