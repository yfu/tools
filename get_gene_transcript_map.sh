#!/usr/bin/env bash

if [[ -z "$1" ]]; then
    echo "Usage: get_gene_transcript_map.sh trinity.fasta"
    exit 1
else
    grep '>' $1 | cut -f1 -d' ' | sed 's/^>//' | sed -r 's/^(TR[0-9]+\|c[0-9]+_g[0-9]+)_i[0-9]+$/\1\t\0/'
fi
