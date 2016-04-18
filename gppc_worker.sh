

# cp ../2016-04-15-paralogs/piRNA.cluster.pachytene.bed12 ./

i=$1
tail -n +${i} piRNA.cluster.pachytene.bed12 | head -n1 > ${i}.bed12
bedtools bed12tobed6 -i ${i}.bed12 > ${i}.bed6

DEG=DEG.amyb_het.uniq.bed12
SRA=SRA.amyb_het.uniq.piRNA.bed2 
bedtools intersect -sorted -a ${DEG} -b ${i}.bed6 -wa > ${DEG}.${i}
bedtools intersect -sorted -a ${SRA} -b ${i}.bed6 -wa > ${SRA}.${i}
bedtools intersect -v -a ${DEG}.${i} -b ~/data/shared/mm10/UCSC.rmsk.sorted.bed -sorted > ${DEG}.x_rmsk.${i}
bedtools intersect -v -a ${SRA}.${i} -b ~/data/shared/mm10/UCSC.rmsk.sorted.bed -sorted > ${SRA}.x_rmsk.${i} 

get_5_watson() {
    awk -v OFS="\t" '{ if($6=="+") { print $1, $2, $2+1, $4, $5, $6 } }' $1
}

get_5_crick() {
    awk -v OFS="\t" '{ if($6=="-") { print $1, $3-1, $3, $4, $5, $6 } }' $1
}

get_5_watson ${DEG}.${i} > ${DEG}.${i}.Watson
get_5_watson ${DEG}.x_rmsk.${i} > ${DEG}.x_rmsk.${i}.Watson

get_5_crick  ${DEG}.${i} > ${DEG}.${i}.Crick
get_5_crick  ${DEG}.x_rmsk.${i} > ${DEG}.x_rmsk.${i}.Crick

get_5_watson ${SRA}.${i} > ${SRA}.${i}.Watson
get_5_watson ${SRA}.x_rmsk.${i} > ${SRA}.x_rmsk.${i}.Watson

get_5_crick  ${SRA}.${i} > ${SRA}.${i}.Crick
get_5_crick  ${SRA}.x_rmsk.${i} > ${SRA}.x_rmsk.${i}.Crick

# ## PP for SRA
gppc_internal.py -r 100 -a ${SRA}.x_rmsk.${i}.Watson -b ${SRA}.x_rmsk.${i}.Crick -p 32 -c mm10.ChromInfo.txt > ${SRA}.x_rmsk.${i}.Watson.vs.${SRA}.x_rmsk.${i}.Crick
pps_plot.R ${SRA}.x_rmsk.${i}.Watson.vs.${SRA}.x_rmsk.${i}.Crick
cat ${SRA}.x_rmsk.${i}.Watson.vs.${SRA}.x_rmsk.${i}.Crick | awk '{ if($1 >= 0 && $1 <=100 && $1 != 9) { sum += $2; sumsq += ($2)^2; n+=1 } if($1==9) { tenth=$2 } } END{ mean=sum/n; std=sqrt((sumsq-sum^2/n)/n); print (tenth-mean)/std  }' > ${SRA}.x_rmsk.${i}.Watson.vs.${SRA}.x_rmsk.${i}.Crick.z10

# ## PP for SRA (excluding RMSK mappers)
gppc_internal.py -r 100 -a ${SRA}.${i}.Watson -b ${SRA}.${i}.Crick -p 32 -c mm10.ChromInfo.txt > ${SRA}.${i}.Watson.vs.${SRA}.${i}.Crick
pps_plot.R ${SRA}.${i}.Watson.vs.${SRA}.${i}.Crick
cat ${SRA}.${i}.Watson.vs.${SRA}.${i}.Crick | awk '{ if($1 >= 0 && $1 <=100 && $1 != 9) { sum += $2; sumsq += ($2)^2; n+=1 } if($1==9) { tenth=$2 } } END{ mean=sum/n; std=sqrt((sumsq-sum^2/n)/n); print (tenth-mean)/std  }' > ${SRA}.${i}.Watson.vs.${SRA}.${i}.Crick.z10

## DEG vs SRA
## 5' DEG-5' SRA distance on Watson
gppc_internal.py -r 100 -a ${DEG}.x_rmsk.${i}.Watson -b ${SRA}.x_rmsk.${i}.Watson -p 32 -c mm10.ChromInfo.txt > ${DEG}.x_rmsk.${i}.Watson.vs.${SRA}.x_rmsk.${i}.Watson
pps_plot.R ${DEG}.x_rmsk.${i}.Watson.vs.${SRA}.x_rmsk.${i}.Watson

## 5' DEG-5' SRA distance on Crick
gppc_internal.py -r 100 -a ${DEG}.x_rmsk.${i}.Crick -b ${SRA}.x_rmsk.${i}.Crick -p 32 -c mm10.ChromInfo.txt > ${DEG}.x_rmsk.${i}.Crick.vs.${SRA}.x_rmsk.${i}.Crick
pps_plot.R ${DEG}.x_rmsk.${i}.Crick.vs.${SRA}.x_rmsk.${i}.Crick

## 5' DEG-5' SRA on different strands
gppc_internal.py -r 100 -a ${DEG}.x_rmsk.${i}.Watson -b ${SRA}.x_rmsk.${i}.Crick -p 32 -c mm10.ChromInfo.txt > ${DEG}.x_rmsk.${i}.Watson.vs.${SRA}.x_rmsk.${i}.Crick
pps_plot.R ${DEG}.x_rmsk.${i}.Watson.vs.${SRA}.x_rmsk.${i}.Crick

gppc_internal.py -r 100 -a ${DEG}.x_rmsk.${i}.Crick -b ${SRA}.x_rmsk.${i}.Watson -p 32 -c mm10.ChromInfo.txt > ${DEG}.x_rmsk.${i}.Crick.vs.${SRA}.x_rmsk.${i}.Watson
pps_plot.R ${DEG}.x_rmsk.${i}.Crick.vs.${SRA}.x_rmsk.${i}.Watson

