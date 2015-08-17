#!/usr/bin/env bash

# Given a one-line bed file, output the conservation score as floats joined by commas.
bed=$1
# You need to have the wib and the wig file in the current dir
ln -s /home/fuy2/data/shared/hg19/phyloP/phyloP100way.wib
ln -s /home/fuy2/data/shared/hg19/phyloP/phyloP100wayAll.wig

hgWiggle -bedFile=${bed} phyloP100wayAll.wig | grep -v "^#" | grep -v "variableStep" > ${bed}.wig.tmp
cut -f2 ${bed}.wig.tmp | awk '{ printf $1 "," }' > ${bed}.csv.tmp
paste ${bed} ${bed}.csv.tmp > ${bed}.with_phylop


