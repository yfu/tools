#!/usr/bin/env bash

# "Columnize" a fasta file for easier processing
# Author: Yu Fu

command -v gawk >/dev/null 2>&1 || { echo >&2 "I require gawk but it's not installed. Aborting."; exit 1; }
cat $1 | gawk '{ if($1 ~ "^>") { match($1, ">(.+)", array); printf "\n" array[1] "\t"} else {printf $1}}' | grep -v -E "^$"
