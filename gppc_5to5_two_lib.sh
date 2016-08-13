#!/usr/bin/env bash

# For cases where we need to calculate 5' to 5' distance on the same strand between two different samples
# Just a wrapper around gppc_internal.py
# Separate a given bed2 file, such that
# a. the first bed2 file contains signals on the Watson strand and only records the 5' position of reads
# b. the second bed2 file contains signals on the Crick strand and only records the 5' position of reads
# Then this script runs gppc_internal.py to get the Ping-Pong score
if [[ "$#" -ne 4 ]]; then
    echo "Usage: gppc_5to5_two_lib.sh XX.ChromInfo.txt one.bed2 two.bed2 out_prefix" 
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
gppc_internal.py -a ${x_watson} -b ${y_watson} -c ${chrominfo} -p 4 2>${out_prefix}.5to5_watson.log >${out_prefix}.5to5_watson.hist.tmp
gppc_internal.py -a ${x_crick} -b ${y_crick} -c ${chrominfo} -p 4 2>${out_prefix}.5to5_crick.log >${out_prefix}.5to5_crick.hist.tmp

# Notice how I transform the distance: you need to plot them on a piece of paper to understand why
awk '{ if(NR==FNR) { d[$1]=$2 } else { d[-$1]+=$2 } } END{ for(i in d) { print i "\t" d[i] } }' ${out_prefix}.5to5_watson.hist.tmp ${out_prefix}.5to5_crick.hist.tmp | sort -k1,1n > ${out_prefix}.5to5.hist


# A possible wrapper
# pp2() {
#     echo "gppc_5to5_two_lib.sh /nfs/fuy2@zlab4/piPipes/common/dm3_chengjian/dm3_chengjian.ChromInfo.txt $1 $2 ${1%.bed2}.${2%.bed2} && pps_plot.R ${1%.bed2}.${2%.bed2}.5to5.hist 0 28 >${1%.bed2}.${2%.bed2}.pp_hist.pdf.log"
# }
