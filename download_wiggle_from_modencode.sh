#!/bin/bash

# Given a url like this: http://data.modencode.org/cgi-bin/cloud_list.pl?accessions=834,835,836,837,838,839&urls=1
# this script will download wiggle files.

url="$1"
curl -s "$url" | grep 'wig' | grep -E -o 'ftp.+' | wget -i -
