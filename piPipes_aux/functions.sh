#!/usr/bin/env bash

function echo2 {
    ISO_8601='%Y-%m-%d %H:%M:%S %Z'
    case $2 in
	error)echo -e $COLOR_RED_BOLD"[`date "+$ISO_8601"`] Error: $1${COLOR_END}" && exit 1 ;;
	warning)echo -e $COLOR_MAGENTA"[`date "+$ISO_8601"`] Warning: $1${COLOR_END}" ;;
	*)echo -e $COLOR_GREEN"[`date "+$ISO_8601"`] $1${COLOR_END}";;
    esac
}
export -f echo2

# produce length distribution for bed2 format
function bed2lendis {
    awk '{l=$3-$2; if (l>m) m=l; c[l]+=$4/$5;}END{for (d=1;d<=m;++d) {printf "%d\t%.0f\n", d, (c[d]?c[d]:0)} }' $1
}
export -f bed2lendis
export bedtools_piPipes=/data/fuy2/piPipes/bin/bedtools_piPipes
