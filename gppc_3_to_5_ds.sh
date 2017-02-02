#!/usr/bin/env bash

# Calculate 3'-to-5' distance on different strands
if [[ "$#" -ne 2 ]]; then
    echo "Usage: gppc_pp.sh XX.ChromInfo.txt your.bed2"
    exit 0
fi

chrominfo=$1
bed2=$2

watson3=${bed2}.Watson3
watson5=${bed2}.Watson5
crick3=${bed2}.Crick3
crick5=${bed2}.Crick5

cat ${bed2} | awk '{ if($6=="+") { print } }' | awk -v OFS="\t" '{ print $1, $2, $2+1, $4, $5, $6, $7 }' > ${watson5}
cat ${bed2} | awk '{ if($6=="+") { print } }' | awk -v OFS="\t" '{ print $1, $3-1, $3, $4, $5, $6, $7 }' > ${watson3}
cat ${bed2} | awk '{ if($6=="-") { print } }' | awk -v OFS="\t" '{ print $1, $3-1, $3, $4, $5, $6, $7 }' > ${crick5}
cat ${bed2} | awk '{ if($6=="-") { print } }' | awk -v OFS="\t" '{ print $1, $2, $2+1, $4, $5, $6, $7 }' > ${crick3}

# There are 16 'chromosomes' in fly genome
gppc_internal.py -a ${watson5} -b ${crick3} -c ${chrominfo} -p 16 2>${bed2}.Watson3toCrick5.ds.log >${bed2}.Watson3toCrick5.ds.phasing_hist
gppc_internal.py -a ${watson3} -b ${crick5} -c ${chrominfo} -p 16 2>${bed2}.Crick5toWatson3.ds.log >${bed2}.Crick5toWatson3.ds.phasing_hist

cat ${bed2}.Watson3toCrick5.ds.phasing_hist | awk '{ if($1==-2) {signal=$2;} if($1>=1 && $1<=29) { sum += $2; sumsq += ($2)^2} } END{ mean=sum/29; std=sqrt((sumsq-sum^2/29)/29);  if(std==0) {print 0} else {print (signal-mean) / std} }' > ${bed2}.Watson3toCrick5.ds.phasing_z-2
cat ${bed2}.Crick5toWatson3.ds.phasing_hist | awk '{ if($1==-2) {signal=$2;} if($1>=1 && $1<=29) { sum += $2; sumsq += ($2)^2} } END{ mean=sum/29; std=sqrt((sumsq-sum^2/29)/29);  if(std==0) {print 0} else {print (signal-mean) / std} }' > ${bed2}.Crick5toWatson3.ds.phasing_z-2

