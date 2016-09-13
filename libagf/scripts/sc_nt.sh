# #!/bin/bash

#set -xe
set -e

#numbers:
NTRIAL=20
NTEST=10
NTRAIN="10 12 14 16 18 20 25 30 35 40 45 50 60 70 80 90 100 120 140 160 180 200"
#NTRAIN="50 100"

OUTFILE=accvntrain.txt

TEST=test
#BINBASE=/home/lenovo/my_software/libmsci/libagf/examples/sample_classes/
BINBASE=/mnt/sdc1/home2/pete/libmsci/libagf/examples/sample_classes/

TRAIN_COMMAND="svm-train -h 0 -c 200"
CLASSIFY_COMMAND=svm-predict

TRAIN=train
MODEL=model.svmmod
OUTPUT=dum

for ((J=0; J<NTRIAL; J++)); do
  ${BINBASE}sample_class -R 1000 2000 ${TEST}$J
  agf2ascii -M ${TEST}$J > ${TEST}$J.txt
done

rm -f $OUTFILE

for nt in $NTRAIN; do

  for ((J=0; J<NTRIAL; J++)); do
    #generate training data:
    ${BINBASE}sample_class $nt $((nt*2)) $TRAIN 
    agf2ascii -M $TRAIN ${TRAIN}.txt
    #train the model:
    echo "${TRAIN_COMMAND} ${TRAIN}.txt ${MODEL}"
    ${TRAIN_COMMAND} ${TRAIN}.txt ${MODEL}
    #make the predictions:
    echo "${CLASSIFY_COMMAND} ${TEST}$J.txt ${MODEL} ${OUTPUT}.txt"
    ${CLASSIFY_COMMAND} ${TEST}$J.txt ${MODEL} ${OUTPUT}.txt
    #convert results to AGF and measure the accuracy:
    echo "svmout2agf ${OUTPUT}.txt $OUTPUT"
    svmout2agf ${OUTPUT}.txt $OUTPUT
    echo -n "$nt " >> ${OUTFILE}
    echo "cls_comp_stats -Hb ${TEST}$J.cls ${OUTPUT} >> ${OUTFILE}"
    cls_comp_stats -Hb ${TEST}$J.cls ${OUTPUT} >> ${OUTFILE}
  done
  echo $NTRIAL
done

for ((J=0; J<NTRIAL; J++)); do
  rm -f ${TEST}$J.*
done

rm -f $TRAIN.*
rm -f $MODEL
rm -f $OUTPUT.*
