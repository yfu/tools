#!/usr/bin/env bash

insert=$1

# Plot the first 18 nt in case the seqeunces are not of the same length
l=18

cat ${insert} | awk -v l=18 '{ for(i=1; i<=$2; i++) { print ">" ++c; print substr($1, 1, l) } }' | weblogo -F pdf -U probability > ${insert}.prob.pdf
cat ${insert} | awk -v l=18 '{ for(i=1; i<=$2; i++) { print ">" ++c; print substr($1, 1, l) } }' | weblogo -F pdf -U bits > ${insert}.entropy.pdf
