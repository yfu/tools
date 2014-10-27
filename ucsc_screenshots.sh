#!/usr/bin/env bash

# Author: Yu Fu
# Usage ucsc_screenshots.sh hgsid slop
# Usage ucsc_screenshots.sh 393202989_aSfr1NFH0Tgwt9jhDFJzj3Ng8sBR 5000
# Obtain the hgsid from the address bar (the hgsid is associated with all the tracks you have loaded)
# hgsid=393202989_aSfr1NFH0Tgwt9jhDFJzj3Ng8sBR
# https://genome.ucsc.edu/cgi-bin/hgTracks?db=dm3&position=chr2R%3A2083755-2447312&hgsid=393202989_aSfr1NFH0Tgwt9jhDFJzj3Ng8sBR
hgsid="$1"
slop="$2"

# db is the genome
db=dm3

ranges=( $(bedtools slop -i ~/data/piPipes/common/dm3/Brennecke.piRNAcluster.bed6.gz -g ~/data/piPipes/common/dm3/dm3.ChromInfo.txt -b ${slop} | cut -f1-4 | grep -v 'X-TAS' | sed -E 's/([^\s]+)\s([^\s]+)\s([^\s]+)\s([^\s]+)/\1%3A\2-\3.\4/') )

for i in "${ranges[@]}"; do
    # position="chr2R%3A2166393-2222282"
    position=$(echo $i | cut -d. -f1)
    dl_fn=$(echo $i | cut -d. -f2).pdf
    curl -s "https://genome.ucsc.edu/cgi-bin/hgTracks?db=$db&position=$position&hgsid=$hgsid" >/dev/null
    dl_page="https://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=${hgsid}&hgt.psOutput=on"
    pdf_url=$(curl -s "$dl_page" | grep "the current browser graphic in PDF" | grep -E -o "\".+\"" | tr -d "\"" | sed 's/../https:\/\/genome.ucsc.edu/')
    # dl_fn=$(echo ${i} | sed -E 's/%3A//')
    # dl_fn=$dl_fn.pdf
    echo "Saving $dl_fn from $pdf_url"
    curl -s -o ${dl_fn} "$pdf_url"
    # https://genome.ucsc.edu/cgi-bin/hgTracks?db=dm3&hgsid=393202989_aSfr1NFH0Tgwt9jhDFJzj3Ng8sBR&position=chr2R%3A2173344-2229233
done
