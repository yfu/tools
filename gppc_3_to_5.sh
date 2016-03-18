#!/usr/bin/env bash

# Calculate 3'-to-5' distance
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
gppc_internal.py -a ${watson3} -b ${watson5} -c ${chrominfo} -p 16 2>${bed2}.phasing.log >${bed2}.Waston3to5.phasing_hist
gppc_internal.py -a ${crick5} -b ${crick3} -c ${chrominfo} -p 16 2>${bed2}.phasing.log >${bed2}.Crick5to3.phasing_hist
