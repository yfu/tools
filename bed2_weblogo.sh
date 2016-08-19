#!/usr/bin/env bash

awk '{ d[$7]+=$4/$5 } END{ for( i in d ) { for(ii=1; ii<=d[i]; ii++) { print ">" i "_" ii; print i } } }' ${1} | weblogo -F pdf -U bits > ${1}.entropy.pdf

awk '{ d[$7]+=$4/$5 } END{ for( i in d ) { for(ii=1; ii<=d[i]; ii++) { print ">" i "_" ii; print i } } }' ${1} | weblogo -F pdf -U probability > ${1}.prob.pdf
