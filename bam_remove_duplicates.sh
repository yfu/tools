picard="java -jar /home/fuy2/src/picard-tools-1.138/picard.jar"
bam=$1
CPU=4
prefix=${bam%.bam}
prefix=${prefix%.sorted}
picard_log=${prefix}.dup_marked.log
dup_marked=${prefix}.dup_marked.bam
metric=${prefix}.dup_marked.metric
dup=${prefix}.dup.bam
x_dup=${prefix}.x_dup.bam
${picard} MarkDuplicates I=${bam} O=${dup_marked} M=${metric} &> ${picard_log}
samtools view -b -f0x400 ${dup_marked} > ${dup}
samtools view -b -F0x400 ${dup_marked} | samtools sort -@ ${CPU} -T ${RANDOM}${RANDOM} -O bam > ${x_dup}

dup_count_log=${prefix}.dup_count.log
samtools view ${bam} | wc -l | awk '{ print "total_mapped_reads\t" $1 }' > ${dup_count_log}
samtools view ${x_dup} | wc -l | awk '{ print "total_mapped_reads_wo_dup\t" $1}' >> ${dup_count_log}
samtools view ${dup} | wc -l | awk '{ print "dup_reads\t" $1}' >> ${dup_count_log}
