#!/usr/bin/env bash

# For cases where we need to calculate ping-pongs between two different samples
# Just a wrapper around gppc_internal.py
# Separate a given bed2 file, such that
# a. the first bed2 file contains signals on the Watson strand and only records the 5' position of reads
# b. the second bed2 file contains signals on the Crick strand and only records the 5' position of reads
# Then this script runs gppc_internal.py to get the Ping-Pong score
if [[ "$#" -ne 4 ]]; then
    echo "Usage: gppc_pp.sh XX.ChromInfo.txt one.bed2 two.bed2 out_prefix" 
    exit 0
fi

chrominfo=$1
x=$2
y=$3
out_prefix=$4

RANDOM=$$
ID=${RANDOM}${out_prefix}${RANDOM}

x_watson=${ID}.${x}.Watson
x_crick=${ID}.${x}.Crick

y_watson=${ID}.${y}.Watson
y_crick=${ID}.${y}.Crick

cat ${x} | awk -v OFS="\t" '{ if($6=="+") { print } }' | awk -v OFS="\t" '{ print $1, $2, $2+1, $4, $5, $6, $7 }' > ${x_watson}
cat ${x} | awk -v OFS="\t" '{ if($6=="-") { print } }' | awk -v OFS="\t" '{ print $1, $3-1, $3, $4, $5, $6, $7 }' > ${x_crick}

cat ${y} | awk -v OFS="\t" '{ if($6=="+") { print } }' | awk -v OFS="\t" '{ print $1, $2, $2+1, $4, $5, $6, $7 }' > ${y_watson}
cat ${y} | awk -v OFS="\t" '{ if($6=="-") { print } }' | awk -v OFS="\t" '{ print $1, $3-1, $3, $4, $5, $6, $7 }' > ${y_crick}

# There are 16 'chromosomes' in fly genome
gppc_internal.py -a ${x_watson} -b ${y_crick} -c ${chrominfo} -p 4 2>${out_prefix}.pp.p1.log >${out_prefix}.p1.pp_hist.tmp
gppc_internal.py -a ${y_watson} -b ${x_crick} -c ${chrominfo} -p 4 2>${out_prefix}.pp.p2.log >${out_prefix}.p2.pp_hist.tmp

paste ${out_prefix}.p1.pp_hist.tmp ${out_prefix}.p2.pp_hist.tmp | awk -v OFS="\t" '{ print $1, $2+$4 }' > ${out_prefix}.pp_hist && rm ${x_watson} ${x_crick} ${y_watson} ${y_crick}


