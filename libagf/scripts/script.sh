#!/bin/bash

#set -xe
set -xe

#numbers:
NTRIAL=20
NTEST=10
MAXSAMPLE=50

while getopts 'gMN:q:s:' DUM
do
  case $DUM in
    g) LOGFLAG=1
      ;;
    M) SVMFLAG=1
      TRAIN_COMMAND=svm_accelerate
      CLASSIFY_COMMAND=classify_c
      ;;
    N) NTRIAL=$OPTARG
      ;;
    q) NTEST=$OPTARG
      ;;
    s) MAXSAMPLE=$OPTARG
      ;;
  esac
done

shift $((OPTIND - 1))

#output file:
OUTFILE=stats.txt

echo $#

if [ $# -eq 3 ]; then
  CONTROL=$1
  TRAINING_DATA=$2
  OUTFILE=$3
  if [ -z $SVMFLAG ]; then
    TRAIN_COMMAND=multi_borders
    CLASSIFY_COMMAND=classify_m
  fi
elif [ $# -eq 2 ]; then
  TRAINING_DATA=$1
  OUTFILE=$2
  if [ -z $SVMFLAG ]; then
    TRAIN_COMMAND=class_borders
    CLASSIFY_COMMAND=classify_b
    echo $TRAIN_COMMAND
  fi
else
  echo "script.sh [-g] [-M] [-N ntrial] [-q ntest] [-s maxsample]\\"
  echo "        [model] train output"
  echo
  echo "  model     = LIBSVM model or multi-borders control file"
  echo "  train     = training data"
  echo "  output    = output skill scores"
  echo
  echo "  g         = use logarithmic progression"
  echo "  M         = use LIBSVM model"
  echo "  ntrial    = number of repeated trials"
  echo "  ntest     = number of sample sizes"
  echo "  maxsample = maximum number of border samples"
  exit 0
fi

MODEL=dum.agf
TEST=testdata
#TEST_POINTS=${TEST}.vec
OUTPUT=dum
#TEST_CLASSES=${TEST}.cls
TRAIN="$CONTROL $TRAINING_DATA"


rm -f ${OUTFILE}
rm -f train.log
rm -f test.log

for ((I=0; I<NTRIAL; I++)); do
  /home/lenovo/my_software/libmsci/libagf/examples/sample_classes/sample_class -R 1000 2000 ${TEST}$I > stuff.txt
done

for ((I=0; I<NTEST; I++)); do
  if [ $LOGFLAG ]; then
    n=$(echo "scale=2; e(($I+1)*l($MAXSAMPLE)/$NTEST)" | bc -l)
  else
    n=$(((I+1)*MAXSAMPLE/NTEST))
  fi
  echo $NTRIAL
  for ((J=0; J<NTRIAL; J++)); do
    echo "(time ${TRAIN_COMMAND} -s $n ${TRAIN} ${MODEL}) 2>> train.log"
    (time ${TRAIN_COMMAND} -s $n ${TRAIN} ${MODEL}) 2>> train.log
    (time ${CLASSIFY_COMMAND} ${MODEL} ${TEST}$I.vec ${OUTPUT} > junk.txt ) 2>> test.log
    echo -n "$n " >> ${OUTFILE}
    cls_comp_stats -Hb ${TEST}$I.cls ${OUTPUT} >> ${OUTFILE}
  done
  echo $NTRIAL
done


