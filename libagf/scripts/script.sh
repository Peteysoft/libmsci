#!/bin/bash

NTRIAL=1
N="1 2 3 4 5 6 7 8 9 10 12 14 16 18 20"

rm -f train.txt
rm -f test.txt
rm -f stats.txt

for n in $N; do
  for ((I=0; I<NTRIAL; I++)); do
 
    (time svm_accelerate -s $n dum.svmmod sctrain dum.agf) 2>> train.txt
    (time classify_s dum.agf sctest.vec dum > junk.txt ) 2>> test.txt
    echo -n "$n " >> stats.txt
    cls_comp_stats -Hb sctest.cls dum >> stats.txt
  done
done

