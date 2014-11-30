#!/usr/bin/env bash

# Convert fastq to fasta. This script does not consider quality score

awk '{ if(NR%4==1) {name=substr($0, 2); print ">" name } if(NR%4==2) {print} }'
