#!/usr/bin/env bash

annotation=$1
bed2=$2
nbin=100

## /data/fuy2/piPipes/common/mm9/UCSC.refSeq.Genes.bed12.gz
bn2=$(basename $annotation)
prefix2=${bn2%.*}
my_anno=${prefix2}.${RANDOM}${RANDOM}${RANDOM}.bed
if [ ${annotation##*.} == gz ] || [ ${annotation##*.} == gzip ]; then
    zcat ${annotation} | cut -f1-6 > ${my_anno}
else
    cat ${annotation} | cut -f1-6 > ${my_anno}
fi

bn=$(basename $bed2)
prefix=${bn%.*}

intersect_tmp_s=${prefix}.vs.${prefix2}.intersect.s.bed
intersect_tmp_as=${prefix}.vs.${prefix2}.intersect.as.bed
awk 'BEGIN{FS=OFS="\t"} { if($6=="+") {$3=$2+1} else{$2=$3-1}; $8=$4/$5; print }' ${bed2} | sort -k1,1 -k2,2n | bedtools groupby -g 1,2,3,6 -c 8 -o sum | awk '{ print $1 "\t" $2 "\t" $3 "\t" "." "\t" $5 "\t" $4 }'|  bedtools intersect -wa -wb -a - -b $my_anno -S | bedtools intersect -wa -wb -a - -b $my_anno -s > ${intersect_tmp_s}
awk 'BEGIN{FS=OFS="\t"} { if($6=="+") {$3=$2+1} else{$2=$3-1};  $8=$4/$5; print }' ${bed2} | sort -k1,1 -k2,2n | bedtools groupby -g 1,2,3,6 -c 8 -o sum | awk '{ print $1 "\t" $2 "\t" $3 "\t" "." "\t" $5 "\t" $4 }' | bedtools intersect -wa -wb -a - -b $my_anno -S > ${intersect_tmp_as}

signal_tmp_s=${prefix}.vs.${prefix2}.fair.signal.s.bed
signal_tmp_as=${prefix}.vs.${prefix2}.fair.signal.as.bed
# $5 is copy/ntm
awk '{ if(NR==FNR) {ntm[$1"\t"$2"\t"$3"\t"$4"\t"$5] += 1} else {print $0 "\t" ntm[$1"\t"$2"\t"$3"\t"$4"\t"$5]} } ' ${intersect_tmp_s} ${intersect_tmp_s} | awk 'BEGIN{FS=OFS="\t"} {print $1 "\t" $2 "\t" $3 "\t" $(NF-3) "\t" $5/$(NF) "\t" $4}' > ${signal_tmp_s}
awk '{ if(NR==FNR) {ntm[$1"\t"$2"\t"$3"\t"$4"\t"$5] += 1} else {print $0 "\t" ntm[$1"\t"$2"\t"$3"\t"$4"\t"$5]} } ' ${intersect_tmp_as} ${intersect_tmp_as} | awk 'BEGIN{FS=OFS="\t"} {print $1 "\t" $2 "\t" $3 "\t" $(NF-3) "\t" $5/$(NF) "\t" $4}' > ${signal_tmp_as}

# Remember to filter out with bin numbers > 100. These are caused by wrongly annotated transcripts.
cat ${my_anno} | awk -v nbin=${nbin} '{ if(NR==FNR) {l[$4]=$3-$2; start[$4]=$2;} else { t=($2-start[$4])/l[$4]*nbin; if($6=="+") {print $4 "\t" t "\t" $5} else {print $4 "\t" nbin-1-t "\t" $5} } }' - ${signal_tmp_s} | awk -v nbin=$nbin '{ bin=int($2); if(bin>nbin || bin<0) {next}; sig[bin]+=$3 } END{ for (i in sig) { print "sense" "\t" i "\t" sig[i] } }'

cat ${my_anno} | awk -v nbin=${nbin} '{ if(NR==FNR) {l[$4]=$3-$2; start[$4]=$2;} else { t=($2-start[$4])/l[$4]*nbin; if($6=="+") {print $4 "\t" t "\t" $5} else {print $4 "\t" nbin-1-t "\t" $5} } }' - ${signal_tmp_as} | awk -v nbin=$nbin '{ bin=int($2); if(bin>nbin || bin<0) {next}; sig[bin]+=$3 } END{ for (i in sig) { print "antisense" "\t" i "\t" sig[i] } }'

rm $intersect_tmp_s $intersect_tmp_as $signal_tmp_s $signal_tmp_as
