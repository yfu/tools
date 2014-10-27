#!/usr/bin/env bash

# Obtain the hgsid from the address bar (the hgsid is associated with all the tracks you have loaded)
hgsid=393202989_aSfr1NFH0Tgwt9jhDFJzj3Ng8sBR
# db is the genome
db=dm3


for position in "chr2R%3A2166393-2222282" "chr2R%3A2566393-2622282"; do
	# position="chr2R%3A2166393-2222282"
	curl -s "https://genome.ucsc.edu/cgi-bin/hgTracks?db=$db&position=$position&hgsid=$hgsid" >/dev/null
	dl_page="https://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=${hgsid}&hgt.psOutput=on"
	pdf_url=$(curl -s "$dl_page" | grep "the current browser graphic in PDF" | grep -E -o "\".+\"" | tr -d "\"" | sed 's/../https:\/\/genome.ucsc.edu/')
	
	curl -O "$pdf_url"
	# https://genome.ucsc.edu/cgi-bin/hgTracks?db=dm3&hgsid=393202989_aSfr1NFH0Tgwt9jhDFJzj3Ng8sBR&position=chr2R%3A2173344-2229233
done
