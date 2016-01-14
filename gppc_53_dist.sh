#!/usr/bin/env bash

# Separate a given bed2 file, such that
# a. the first bed2 file contains signals on the Watson strand and only records the 5' position of reads
# b. the second bed2 file contains signals on the Crick strand and only records the 5' position of reads
# Then this script runs gppc_internal.py to get the Ping-Pong score
if [[ "$#" -ne 2 ]]; then
    echo "Usage: gppc_pp.sh XX.ChromInfo.txt your.bed2"
    exit 0
fi

chrominfo=$1
bed2=$2
watson=${bed2}.Watson
crick=${bed2}.Crick

watson_5=${bed2}.Watson5
watson_3=${bed2}.Watson3
crick_5=${bed}.Crick5
crick_3=${bed2}.Crick3

cat ${bed2} | grep "+" | awk -v OFS="\t" '{ print $1, $2, $2+1, $4, $5, $6, $7 }' > ${watson_5}
cat ${bed2} | grep "-" | awk -v OFS="\t" '{ print $1, $3-1, $3, $4, $5, $6, $7 }' > ${crick_5}

cat ${bed2} | grep "+" | awk -v OFS="\t" '{ print $1, $3-1, $3, $4, $5, $6, $7 }' > ${watson_3}
cat ${bed2} | grep "-" | awk -v OFS="\t" '{ print $1, $2, $2+1, $4, $5, $6, $7 }' > ${crick_3}

gppc_internal.py -a ${watson_3} -b ${watson_5} -c ${chrominfo} -p 16 2>${bed2}.Watson.53_dist.log >${bed2}.Watson.53_dist_hist
gppc_internal.py -a ${crick_5} -b ${crick_3} -c ${chrominfo} -p 16 2>${bed2}.Crick.53_dist.log >${bed2}.Crick.53_dist_hist

