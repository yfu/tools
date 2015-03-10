#!/usr/bin/env bash

cat $1 | awk '{ l[$3-$2] += $4/$5 } END{ for(ll in l) { print ll "\t" l[ll] } }'
