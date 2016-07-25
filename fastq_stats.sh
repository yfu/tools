#!/usr/bin/env bash

n_read=$(zless ${1} | wc -l | awk '{ print $1/4 }')
echo -e "n_read\t${n_read}"

n_base=$(zless ${1} | awk '{ if(NR%4==2) { s+=length($0) } } END{ print s }')
echo -e "n_base\t${n_base}"
