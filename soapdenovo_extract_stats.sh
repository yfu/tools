#!/usr/bin/env bash

# Usage: soapdenovo_extract_stats.sh soapdenovo.stderr
i=$1

fn=${i}
grep 'contig(s) longer than' $fn | grep 'There are'
grep 'The longest length is' $fn
grep 'contig(s) longer than' $fn | grep 'output'

grep 'Scaffold number' $fn

grep ' In-scaffold contig number' $fn
cat $fn | grep 'Total scaffold length'
cat $fn | grep 'Average scaffold length'
cat $fn | grep ' Filled gap number'
cat $fn | grep 'Longest scaffold'
cat $fn | grep 'Scaffold and singleton number'
cat $fn | grep 'Scaffold and singleton length'
cat $fn | grep ' Average length'
cat $fn | grep -E '^ N50'
cat $fn | grep -E '^ N90'

cat $fn | grep -E '^ Average insert size'
echo

