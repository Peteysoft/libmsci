#!/bin/bash

#set -xe
set -e

#keep this a bash script (cygwin bash version of time doesn't support -o option):
TIME=/usr/bin/time

#numbers:
NTRIAL=20
NTEST=10
MINSAMPLE=5
MAXSAMPLE=50
FRAC=0.2

#if we want to do it with the synthetic test classes:
SAMPLE_DIR=/home/lenovo/my_software/libmsci/libagf/examples/sample_classes
#fixed numbers (have to fix that...)
NSCTEST=5000
#NSCTRAIN=10000

while getopts 'gKMZN:f:k:q:s:S:t:W:' DUM
do
  case $DUM in
    g) LOGFLAG=1
      ;;
    K) KFLAG=1
      ;;
    M) SVMFLAG=1
      TRAIN_COMMAND=svm_accelerate
      CLASSIFY_COMMAND=classify_s
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
    s) MINSAMPLE=$OPTARG
      ;;
    S) MAXSAMPLE=$OPTARG
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
elif [ $# -eq 1 ]; then
  OUTFILE=$1
  if [ -z $SVMFLAG ]; then
    if [ -z $ZFLAG ]; then
      CLASSIFY_COMMAND=classify_m
    fi
  fi
else
  echo "   MULTI-BORDERS NUMBER OF BORDERS"
  echo
  echo "mb_nb.sh [-g] [-M] [-K] [-N ntrial] [-q ntest] [-s minsample]\\"
  echo "        [-S maxsample] [-f frac] [[model] train] output"
  echo
  echo "  model     = LIBSVM model or multi-borders control file"
  echo "  train     = training data--if omitted uses synthetic test classes"
  echo "  output    = output skill scores"
  echo
  echo "  -g        = use logarithmic progression"
  echo "  -K        = keep temporary files"
  echo "  -M        = use LIBSVM model--native multi-class"
  echo "  -Z        = use LIBSVM model--libAGF multi-class"
  echo
  echo "  ntrial    = number of repeated trials [$NTRIAL]"
  echo "  ntest     = number of sample sizes [$NTEST]"
  echo "  minsample = maximum number of border samples [$MINSAMPLE]"
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

#generate test and training data:
for ((I=0; I<NTRIAL; I++)); do
  if [ -z ${TRAINING_DATA} ]; then
    ${SAMPLE_DIR}/sample_class -R $NSCTEST $((2*NSCTEST)) $BASE.$I.tst > $BASE.$I.tst.txt
  else
    agf_preprocess -zf $FRAC ${TRAINING_DATA} $BASE.$I.trn $BASE.$I.tst
  fi
done

#we need a control file for the sample clases because we want to measure overhead
#in "multi-borders" classification:
if [ -z ${TRAINING_DATA} ]; then
  echo "$BASE 0 / 1; {0 1}" > $MODEL
fi

for ((I=0; I<NTEST; I++)); do
  if [ $LOGFLAG ]; then
    n=$(echo "a=e(l($MINSAMPLE)+$I*(l($MAXSAMPLE)-l($MINSAMPLE))/($NTEST-1)); scale=0; (a+1)/1" | bc -l)
  else
    n=$((MINSAMPLE+I*(MAXSAMPLE-MINSAMPLE)/(NTEST-1)))
  fi

  for ((J=0; J<NTRIAL; J++)); do
    TEST="$BASE.$J.tst"
    TIMING="$BASE.$J.tm"
    TRAIN="$CONTROL $BASE.$J.trn"
    if [ -z ${TRAINING_DATA} ]; then
      echo "${SAMPLE_DIR}/sc_borders -s $n $BASE"
      ${SAMPLE_DIR}/sc_borders -s $n $BASE
    else
      echo $TRAIN
      echo ${TRAIN_COMMAND}
      echo "$TIME -o $TIMING ${TRAIN_COMMAND} -s $n ${ALLOPT} ${TRAIN} ${MODELFILEBASE} ${MODEL}"
      $TIME -o $TIMING ${TRAIN_COMMAND} -s $n ${ALLOPT} ${TRAIN} ${MODELFILEBASE} ${MODEL}
    fi
    #echo "(time ${CLASSIFY_COMMAND} ${MODEL} ${TEST}.vec ${OUTPUT} > junk.txt ) 2>> test.log"
    echo "$TIME -o ${TEST}.tm ${CLASSIFY_COMMAND} ${MODEL} ${TEST}.vec ${OUTPUT}"
    $TIME -o ${TEST}.tm \
	    ${CLASSIFY_COMMAND} ${MODEL} ${TEST}.vec ${OUTPUT}

    echo -n "$n " >> ${OUTFILE}
    echo "cls_comp_stats -Hb ${TEST}.cls ${OUTPUT} >> ${OUTFILE}"
    cls_comp_stats -Hb ${TEST}.cls ${OUTPUT} | tr -d "\n" >> ${OUTFILE}
    echo -n " " >> $OUTFILE
    if [ ${TRAINING_DATA} ]; then
      grep -o "[0-9]*\.[0-9]*" $TIMING | sed '1q;d' | tr -d '\n' >> $OUTFILE
      echo -n " " >> $OUTFILE
      grep -o "[0-9]*\.[0-9]*" $TIMING | sed '2q;d' | tr -d '\n' >> $OUTFILE
      echo -n " " >> $OUTFILE
    fi
    grep -o "[0-9]*\.[0-9]*" $TEST.tm | sed '1q;d' | tr -d '\n' >> $OUTFILE
    echo -n " " >> $OUTFILE
    grep -o "[0-9]*\.[0-9]*" $TEST.tm | sed '2q;d' >> $OUTFILE
  done
  echo $NTRIAL
done

#clean up, remove temporary files:
if [ -z $KFLAG ]; then
  rm -f $BASE.*
fi

