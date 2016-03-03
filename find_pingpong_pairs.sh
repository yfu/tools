#!/usr/bin/env bash


watson=$1
crick=$2
awk '{ if(NR==FNR) {all_w5[$2]=$0} else { c5=$3-1; tmp=c5-9; if(tmp in all_w5) { print all_w5[tmp] "\t" $0 } } }' ${watson} ${crick} | sort -k2,2n


