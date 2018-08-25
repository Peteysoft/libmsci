#!/bin/bash

DATASETS="letter pendigits usps segment sat urban shuttle humidity"

rm -f results.txt
for SET in $DATASETS; do
  ../statlog/sum_col 0 1 / 2 3 / 4 5 / 6 7 / 8 9 / 10 11 / 12 13 / 14 15 / 16 17 / 18 19 / 20 21 / 22 23 / 24 25 / 26 27 / 29 / 31 / 33 / 35 / 37 / 39 / 41 / 43 / 44 / 45 / 46 / 47 / 48 / 49 / 50 / 51 / 52 / 53 / 54 / 55 / 56 / 57 < $SET.txt | ../statlog/calc_stats >> results.txt
done;

./write_results -F 0 6 $DATASETS < results.txt
echo
./write_results -F 6 8 $DATASETS < results.txt
echo
./write_results 14 8 $DATASETS < results.txt
echo
./write_results 22 6 $DATASETS < results.txt
echo
./write_results 28 8 $DATASETS < results.txt

