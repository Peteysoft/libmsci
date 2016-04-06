#!/bin/bash

#set -xe
set -e

#numbers:
NTRIAL=20
NTEST=10
MAXSAMPLE=50
FRAC=0.2

while getopts 'gKMZN:f:k:q:s:t:W:' DUM
do
  case $DUM in
    g) LOGFLAG=1
      ;;
    K) KFLAG=1
      ;;
    M) SVMFLAG=1
      TRAIN_COMMAND=svm_accelerate
      CLASSIFY_COMMAND=classify_c
      ;;
    Z) ZFLAG=1
      TRAIN_COMMAND="multi_borders -Z"
      CLASSIFY_COMMAND=classify_m
      ;;
    N) NTRIAL=$OPTARG
      ;;
    f) FRAC=$OPTARG
      ;;
    q) NTEST=$OPTARG
      ;;
    s) MAXSAMPLE=$OPTARG
      ;;
    t) ALLOPT="$ALLOPT -t $OPTARG"
      ;;
    k) AGFOPT="$AGFOPT -k $OPTARG"
      ;;
    W) AGFOPT="$AGFOPT -W $OPTARG"
      ;;
  esac
done

shift $((OPTIND - 1))

#all temporary files:
BASE=mb$RANDOM.tmp

if [ $# -eq 3 ]; then
  CONTROL=$1
  TRAINING_DATA=$2
  OUTFILE=$3
  if [ -z $SVMFLAG ]; then
    if [ -z $ZFLAG ]; then
      TRAIN_COMMAND="multi_borders $AGFOPT"
      CLASSIFY_COMMAND=classify_m
      MODELFILEBASE=$BASE
    fi
  fi
elif [ $# -eq 2 ]; then
  TRAINING_DATA=$1
  OUTFILE=$2
  if [ -z $SVMFLAG ]; then
    if [ -z $ZFLAG ]; then
      TRAIN_COMMAND="class_borders $AGFOPT"
      CLASSIFY_COMMAND=classify_b
    fi
  fi
else
  echo "   MULTI-BORDERS NUMBER OF BORDERS"
  echo
  echo "mb_nb.sh [-g] [-M] [-K] [-N ntrial] [-q ntest] [-s maxsample]\\"
  echo "        [-f frac] [model] train output"
  echo
  echo "  model     = LIBSVM model or multi-borders control file"
  echo "  train     = training data"
  echo "  output    = output skill scores"
  echo
  echo "  g         = use logarithmic progression"
  echo "  K         = keep temporary files"
  echo "  M         = use LIBSVM model--native multi-class"
  echo "  Z         = use LIBSVM model--libAGF multi-class"
  echo
  echo "  ntrial    = number of repeated trials [$NTRIAL]"
  echo "  ntest     = number of sample sizes [$NTEST]"
  echo "  maxsample = maximum number of border samples [$MAXSAMPLE]"
  echo "  frac      = fraction to reserve for testing [$FRAC]"
  exit 0
fi

MODEL=$BASE.agf
OUTPUT=$BASE

if [ $ZFLAG ]; then
  MODELFILEBASE=$BASE
fi

rm -f ${OUTFILE}
rm -f train.log
rm -f test.log

for ((I=0; I<NTRIAL; I++)); do
  #/home/lenovo/my_software/libmsci/libagf/examples/sample_classes/sample_class -R 1000 2000 ${TEST}$I > stuff.txt
  agf_preprocess -zf $FRAC ${TRAINING_DATA} $BASE.$I.trn $BASE.$I.tst
done

for ((I=0; I<NTEST; I++)); do
  if [ $LOGFLAG ]; then
    n=$(echo "e(($I+1)*l($MAXSAMPLE)/$NTEST)" | bc -l)
  else
    n=$(((I+1)*MAXSAMPLE/NTEST))
  fi

  for ((J=0; J<NTRIAL; J++)); do
    TRAIN="$CONTROL $BASE.$J.trn"
    TEST="$BASE.$J.tst"
    echo "(time ${TRAIN_COMMAND} -s $n ${ALLOPT} ${TRAIN} ${MODELFILEBASE} ${MODEL}) 2>> train.log"
    (time ${TRAIN_COMMAND} -s $n ${ALLOPT} ${TRAIN} ${MODELFILEBASE} ${MODEL}) 2>> train.log
    echo "(time ${CLASSIFY_COMMAND} ${MODEL} ${TEST}.vec ${OUTPUT} > junk.txt ) 2>> test.log"
    (time ${CLASSIFY_COMMAND} ${MODEL} ${TEST}.vec ${OUTPUT} > junk.txt ) 2>> test.log
    echo -n "$n " >> ${OUTFILE}
    echo "cls_comp_stats -Hb ${TEST}.cls ${OUTPUT} >> ${OUTFILE}"
    cls_comp_stats -Hb ${TEST}.cls ${OUTPUT} >> ${OUTFILE}
  done
  echo $NTRIAL
done

#clean up, remove temporary files:
if [ -z $KFLAG ]; then
  rm -f $BASE.*
fi

