
# idx=../../../data/Ponten/testis/Ponten.RSQ.testis_7e.trinity/Trinity.fasta
idx=/data/fuy2/pachytene/data/Ponten/testis/Ponten.RSQ.testis_7e.trinity/Trinity.fpkm1.trimmed_header.fasta
insert=../../../data/Xing/Xing.SRA.human.testis.rep1/input_read_files/Xing.SRA.human.testis.rep1.trimmed.x_rRNA.x_hairpin.insert

prefix=Xing.SRA.human.testis.rep1.trimmed.x_rRNA.x_hairpin.to.Ponten.RSQ.testis_7e.trinity
al=${prefix}.al
un=${prefix}.un
log=${prefix}.bowtie.log
output=${prefix}.bed
bowtie -p 16 -r --best -strata --all -v 1 -S --al ${al} --un ${un} ${idx} ${insert} 2>${log}| samtools view -F 0x4 -b - | bedtools bamtobed -i - > ${output}
