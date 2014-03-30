#!/usr/bin/env BASH

# Download all papers of the Nature Review Genetics
root="http://www.nature.com"
all_cate=$(curl -s "http://www.nature.com/nrg/series/index.html" | grep -E -o "nrg/series/[a-zA-Z/]+/index.html" | tr '\n' ' ')
for a in $all_cate; do
    echo "Searching PDFs on this page:" $a
    cate_name=$(echo $a | awk -F'/' '{print $3}')
    mkdir -p $cate_name
    paths=$(curl -s ${root}/${a} | grep -E -o  "nrg.*pdf" | tr "\n" " ")
    cd $cate_name
    for p in $paths; do
	# echo ${root}/${a}/${p}
	# echo ${root}/${p} >> tmp/papers.list
	curl -OL ${root}/${p}
    done;
    cd -
done;

