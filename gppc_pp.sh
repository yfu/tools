#!/usr/bin/env bash

# Just a wrapper around gppc_internal.py
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
cat ${bed2} | grep "+" | awk -v OFS="\t" '{ print $1, $2, $2+1, $4, $5, $6, $7 }' > ${watson}
cat ${bed2} | grep "-" | awk -v OFS="\t" '{ print $1, $3-1, $3, $4, $5, $6, $7 }' > ${crick}
# There are 16 'chromosomes' in fly genome
gppc_internal.py -a ${watson} -b ${crick} -c ${chrominfo} -p 16 2>${bed2}.pp.log >${bed2}.pp_hist
