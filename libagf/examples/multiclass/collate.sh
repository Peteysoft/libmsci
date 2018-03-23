#!/bin/bash

DATASETS="letter pendigits sat segment shuttle usps"

#rm -f results.txt
#for $SET in $DATASETS; do
#  ../statlog/calc_stats < $SET.txt >> results.txt
#done;

./write_results 0 12 $DATASETS < results.txt
echo
./write_results 12 4 $DATASETS < results.txt
echo
./write_results 16 6 $DATASETS < results.txt

