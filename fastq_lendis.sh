#!/usr/bin/env bash

# cat $1 | awk '{ l[$3-$2] += $4/$5 } END{ for(ll in l) { print ll "\t" l[ll] } }'
cat ${i} | awk '{ if(NR%4==2) { a[length($0)]+=1 } } END{ for(i in a) { print i "\t" a[i] } }' 
