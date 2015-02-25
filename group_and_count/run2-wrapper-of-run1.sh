for i in 00 02 04 07 10 12 14 17 20 42; do
    dpp=${i}dpp
    bed2=~/data/prepachytene/results/2014-11-21/SRA_mm10g/${dpp}/genome_mapping/Zamore.SRA.total.${dpp}.testis.trimmed.x_rRNA.x_hairpin.mm10gv1.all.piRNA.bed2
    # head -n1000 ${bed2} > test.bed
    # bed2=test.bed
    echo "bash run1-group-bed2.sh ${bed2} ${dpp}"
    bed2=/data/fuy2/prepachytene/results/2014-11-21/SRA_mm10g/${dpp}_unox/genome_mapping/Zamore.SRA.total.unox.${dpp}.testis.trimmed.x_rRNA.x_hairpin.mm10gv1.all.piRNA.bed2
    echo "bash run1-group-bed2.sh ${bed2} ${dpp}_unox"
done
