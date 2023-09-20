#!/usr/bin/env bash

FILE=$1
COLNAME=$2
OUT=$3

less ${FILE} | head -n100 | grep -A1 '%' > temp_header
less ${FILE} | grep -v '%' | tail -n +2 | sort -k1,1n -k2,2n -t' ' | cat temp_header - | gzip > ${OUT}
