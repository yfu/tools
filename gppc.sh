#/usr/bin/env bash

# Seperate a bed2 file by chromosomes and run gppc_core.py for them in parallel
# to reduce memory requirement and running time

cl_fn=$1
bed2=$2
cat $bed2 | awk -v bed2=$bed2 '{ d[$1]=0; chrom=$1; print $0 > "tmp."bed2"."chrom } END{ for(i in d) { print i > "tmp.filelist" } }' 

for i in 
for i in $(cat tmp.filelist); do
    echo "gppc_core.py ${cl_fn} tmp.${bed2}.${i}"
done
