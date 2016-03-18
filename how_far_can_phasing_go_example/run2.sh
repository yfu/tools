#!/usr/bin/env bash

for i in $(seq 0 500 10000); do
    offset=-${i}
    awk -v OFS="\t" -v offset=${offset} '{ $2+=offset; $3+=offset; print $0 }' GSV6_upstream_in_zip.dm3.2k.bed > GSV6_upstream_in_zip.dm3.2k.offset${offset}.bed

    offset=${i}
    awk -v OFS="\t" -v offset=${offset} '{ $2+=offset; $3+=offset; print $0 }' GSV6_downstream_in_zip.dm3.2k.bed > GSV6_downstream_in_zip.dm3.2k.offset${offset}.bed
done


