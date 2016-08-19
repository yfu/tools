#!/usr/bin/env bash

<<<<<<< HEAD
# Based on gppc_pp.sh
# This script calculates A.watson vs B.Crick and also B.Watson and A.Crick
# This is for cases when A and B are different
# If A and B are the same, i.e. they come from the same library, two results are the same
# Just a wrapper around gppc_internal.py
# Process bed2 files 
=======
# For cases where we need to calculate ping-pongs between two different samples
# Just a wrapper around gppc_internal.py
>>>>>>> af71a6bdd0209f2bd029f8f4f269efe7b606f919
# Separate a given bed2 file, such that
# a. the first bed2 file contains signals on the Watson strand and only records the 5' position of reads
# b. the second bed2 file contains signals on the Crick strand and only records the 5' position of reads
# Then this script runs gppc_internal.py to get the Ping-Pong score
<<<<<<< HEAD
if [[ "$#" -ne 2 ]]; then
    echo "Usage: gppc_pp.sh XX.ChromInfo.txt your.bed2"
=======
if [[ "$#" -ne 4 ]]; then
    echo "Usage: gppc_pp.sh XX.ChromInfo.txt one.bed2 two.bed2 out_prefix" 
>>>>>>> af71a6bdd0209f2bd029f8f4f269efe7b606f919
    exit 0
fi

chrominfo=$1
<<<<<<< HEAD
A=$2
B=$3
Awatson=${bed2}.Watson
Acrick=
Bcrick=${bed2}.Crick

cat ${bed2} | awk -v OFS="\t" '{ if($6=="+") { print } }' | awk -v OFS="\t" '{ print $1, $2, $2+1, $4, $5, $6, $7 }' > ${watson}
# cat ${bed2} | grep "+" | awk -v OFS="\t" '{ print $1, $2, $2+1, $4, $5, $6, $7 }' > ${watson}
cat ${bed2} | awk -v OFS="\t" '{ if($6=="-") { print } }' | awk -v OFS="\t" '{ print $1, $3-1, $3, $4, $5, $6, $7 }' > ${crick}
# There are 16 'chromosomes' in fly genome
gppc_internal.py -a ${watson} -b ${crick} -c ${chrominfo} -p 16 2>${bed2}.pp.log >${bed2}.pp_hist
gppc_internal.py -a ${watson} -b ${crick} -c ${chrominfo} -p 16 2>${bed2}.pp.log >${bed2}.pp_hist
=======
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


>>>>>>> af71a6bdd0209f2bd029f8f4f269efe7b606f919
