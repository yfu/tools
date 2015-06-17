# fuy2@ummsres18:~/data/pachytene/results/2015-06-17$ samtools view -h ../../data/Xuan/Zamore.RSQ.pi9H_pi17H_RD13.73dpp.031415.rep1/genome_mapping/Zamore.RSQ.pi9H_pi17H_RD13.73dpp.031415.rep1.x_rRNA.mm10g.sorted.bam | head -n100000 | samtools view -Sb - > test1.bam
# fuy2@ummsres18:~/data/pachytene/results/2015-06-17$ samtools view -h ../../data/Xuan/Zamore.RSQ.pi9M_pi17M_RD13.73dpp_031415.rep1/genome_mapping/Zamore.RSQ.pi9M_pi17M_RD13.73dpp_031415.rep1.x_rRNA.mm10g.sorted.bam | head -n100000 | samtools view -Sb - > test2.bam

# bam1=test1.bam
# bam2=test2.bam
bam1=../../data/Xuan/Zamore.RSQ.pi9H_pi17H_RD13.73dpp.031415.rep1/genome_mapping/Zamore.RSQ.pi9H_pi17H_RD13.73dpp.031415.rep1.x_rRNA.mm10g.sorted.bam
bam2=../../data/Xuan/Zamore.RSQ.pi9M_pi17M_RD13.73dpp_031415.rep1/genome_mapping/Zamore.RSQ.pi9M_pi17M_RD13.73dpp_031415.rep1.x_rRNA.mm10g.sorted.bam
label1=pi9H_pi17H
label2=pi9M_pi17M
output=cuffdiff_double_het_mut
mkdir -p ${output}
log=${output}/${output}.log
cuffdiff -o ${output} -L pi9H_pi17H,pi9M_pi17M -p 32 --compatible-hits-norm --frag-bias-correct ~/data/piPipes/common/mm10g/mm10g.fa --multi-read-correct --FDR 0.05 --library-type fr-firststrand --library-norm-method geometric --no-update-check ~/data/piPipes/common/mm10g/mm10g.genes.gtf  ${bam1} ${bam2} &> ${log}

