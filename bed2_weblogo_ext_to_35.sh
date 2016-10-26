input=$1
# input=from_animal/Zamore.SRA.ovary.ox.rep1.ID4/genome_mapping/Zamore.SRA.ovary.ox.rep1.ID4.x_rRNA.x_hairpin.hi5_v0.8.2v0.all.piRNA.bed2

cat ${input} | awk -v OFS="\t" '{ if($6=="+") { print $1, $2, $2+1, $4/$5, 0, $6 } else { print $1, $3-1, $3, $4/$5, "0", $6 } }' | bedtools slop -i - -g ~/data/cl/shared/draft_genome/hi5_v0.8.2.ChromInfo.txt -s -r 34 -l 0 | awk '$3-$2==35' | bedtools getfasta -s -bed - -fi ~/data/cl/shared/draft_genome/hi5_v0.8.2.fasta -fo - -name -tab | awk '{ d[$2]+=$1 } END{ for(i in d) { printf "%s\t%d\n", i, d[i]+0.5 } }' > ${input}.weblogo.tmp
awk '{ for(i=0; i<$2; i++) { print ">" ++j; print $1 } }' ${input}.weblogo.tmp | weblogo -P "" -F pdf -U bits > ${input}.weblogo.entropy.pdf
awk '{ for(i=0; i<$2; i++) { print ">" ++j; print $1 } }' ${input}.weblogo.tmp | weblogo -P "" -F pdf -U probability > ${input}.weblogo.probability.pdf
