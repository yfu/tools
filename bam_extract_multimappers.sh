#!/usr/bin/env bash

bam=$1
samtools view -h ${bam} | awk '{ if($0 ~ /^@/) {print} else { match($0, /NH:i:([0-9]+)/, a); nh=a[1]; if(length(nh)==0) { print "Error: the following line does not contain NH tag!"; print $0 } else { if( nh > 1) { print $0 } } } }' | samtools view -b - 
