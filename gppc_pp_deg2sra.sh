#!/usr/bin/env bash

# For cases where we need to calculate ping-pongs between two different samples
# Just a wrapper around gppc_internal.py
# Separate a given bed2 file, such that
# a. the first bed2 file contains signals on the Watson strand and only records the 5' position of reads
# b. the second bed2 file contains signals on the Crick strand and only records the 5' position of reads
# Then this script runs gppc_internal.py to get the Ping-Pong score
if [[ "$#" -ne 4 ]]; then
    echo "Usage: gppc_pp.sh XX.ChromInfo.txt deg.bed2 sra.bed2 out_prefix" 
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

# Always use deg as the referece (-a)
gppc_internal.py -a ${x_watson} -b ${y_crick} -c ${chrominfo} -p 4 2>${out_prefix}.pp.p1.log >${out_prefix}.p1.tmp
gppc_internal.py -a ${x_crick} -b ${y_watson} -c ${chrominfo} -p 4 2>${out_prefix}.pp.p2.log >${out_prefix}.p2.tmp

# Need to flip the sign of column one in p2
awk '{ if(NR==FNR) { d[$1]=$2 } else { d[-$1]+=$2 } } END{ for(i in d) { print i "\t" d[i] } }' ${out_prefix}.p1.tmp ${out_prefix}.p2.tmp | sort -k1,1n > ${out_prefix}.sra2deg.hist
pps_plot.R ${out_prefix}.sra2deg.hist &> ${out_prefix}.sra2deg.hist.log



