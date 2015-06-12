#!/usr/bin/env awk
# This is used to parse a bed12 file and generate the length of CDS for each transcript

BEGIN{ OFS="\t"; } { chr=$1; start=$2; end=$3; name=$4; thickStart=$7; thickEnd=$8; split($11, exon_len, ","); split($12, exon_start, ","); n=length(exon_len); for(i=1; i<n; i++) { s=exon_start[i]+start; e=exon_start[i]+start+exon_len[i]; if(s<thickStart) { s=thickStart }; if(e>thickEnd) { e=thickEnd }; if(e>s) { l[name] += e-s } } } END{ for (i in l) { print i, l[i] } }
