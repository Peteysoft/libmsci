#!/bin/bash

set -xe
#set -e

#numbers:
NTRIAL=20
NTEST=10
NTRAIN="10 12 14 16 18 20 25 30 35 40 45 50 60 70 80 100"

OUTFILE=accvntrain.txt

TEST=test
BINBASE=/home/lenovo/my_software/libmsci/libagf/examples/sample_classes/

TRAIN_COMMAND=svm-train
CLASSIFY_COMMAND=svm-predict

TRAIN=train
MODEL=model.svmmod
OUTPUT=dum

for ((J=0; J<NTRIAL; J++)); do
  ${BINBASE}sample_class -R 1000 2000 ${TEST}$I > stuff.txt
  agf2ascii -M ${TEST}$I > ${TEST}$J.txt
done

for nt in $NTRAIN; do
  ${BINBASE}sample_class $nt $((nt*2)) $TRAIN 
  agf2ascii -M $TRAIN ${TRAIN}.txt

  for ((J=0; J<NTRIAL; J++)); do
    echo "(time ${TRAIN_COMMAND} ${TRAIN}.txt ${MODEL}) 2>> train.log"
    (time ${TRAIN_COMMAND} ${TRAIN}.txt ${MODEL}) 2>> train.log
    echo "(time ${CLASSIFY_COMMAND} ${TEST}$J.txt ${MODEL} ${OUTPUT}.txt > junk.txt ) 2>> test.log"
    (time ${CLASSIFY_COMMAND} ${TEST}$J.txt ${MODEL} ${OUTPUT}.txt > junk.txt ) 2>> test.log
    svmout2agf ${OUTPUT}.txt $OUTPUT
    echo -n "$n " >> ${OUTFILE}
    echo "cls_comp_stats -Hb ${TEST}.cls ${OUTPUT} >> ${OUTFILE}"
    cls_comp_stats -Hb ${TEST}.cls ${OUTPUT} >> ${OUTFILE}
  done
  echo $NTRIAL
done

