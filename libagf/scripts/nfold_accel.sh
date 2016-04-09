#!/bin/bash

set -e

#BASE=/home/lenovo/my_software/libmsci/libagf/examples/sample_classes

N="10 20 30 40 50 75 100 200 300 400 500 750 1000 2000 3000 4000 5000 7500 10000"
#N="10 50 100 500 1000 5000 10000"

#$BASE/sample_class -R 1000 2000 sctest > sctest.lvq
#agf2ascii -M sctest sctest.svm

#numbers:
METHOD=0
NTRIAL=10
NSAMPLE=100

while getopts 'GKMq:s:c:d:g:h:Q:r:t:' DUM
do
  case $DUM in
    K) KFLAG=1
      ;;
    M) SVMFLAG=1
      ;;
    d) NTRIAL=$OPTARG
      ;;
    s) NSAMPLE=$OPTARG
      ;;
    c) SVMOPT="$SVMOPT -c $OPTARG"
      ;;
    d) SVMOPT="$SVMOPT -d $OPTARG"
      ;;
    g) SVMOPT="$SVMOPT -g $OPTARG"
      ;;
    h) SVMOPT="$SVMOPT -h $OPTARG"
      ;;
    s) METHOD=$OPTARG
      ;;
    r) SVMOPT="$SVMOPT -r $OPTARG"
      ;;
    t) SVMOPT="$SVMOPT -t $OPTARG"
      ;;
  esac
done

shift $((OPTIND - 1))

#base name for intermediate files:
BASE=accel$RANDOM.tmp

echo $#

if [ $# -eq 2 ]; then
  CONTROL=$1
  TRAINING_DATA=$2
  RESULTBASE=$BASE.results
elif [ $# -eq 3 ]; then
  CONTROL=$1
  TRAINING_DATA=$2
  RESULTBASE=$3
else
  echo "   N-FOLD FOR ACCELERATED SVM"
  echo
  echo "nfold_accel.sh [-M] [-K] [-s nsamples] [-d ntrial] [-Q method]\\"
  echo "           [svm-options] control train out"
  echo
  echo "  svm-options  = options to pass to svm-train"
  echo "  control    = input control file"
  echo "  train     = training data"
  echo "  out       = base name for output files"
  echo
  echo "  K         = keep temporary files"
  echo "  M         = data is in LIBSVM format"
  echo
  echo "  nsamples  = number border samples [$NSAMPLE]"
  echo "  ntrial    = number of folds [$NTRIAL]"
  echo "  method    = method of solving for conditional probabilities [$METHOD]"
  exit 0
fi


if test $SVMFLAG; then
  #this is kind of dum:
  echo "svm2agf ${TRAINING_DATA} $BASE.0"
  svm2agf ${TRAINING_DATA} $BASE.0
  TRAINING_DATA=$BASE.0
fi

echo "agf_preprocess -zd $NTRIAL ${TRAINING_DATA} $BASE.tst"
agf_preprocess -zd $NTRIAL ${TRAINING_DATA} $BASE.tst

for ((I=0; I<NTRIAL; I++)); do
  ISTRING=$(printf "%2.2d" $I)
  echo "agf2ascii -M $BASE.tst-$ISTRING $BASE.tst-$ISTRING.svm"
  agf2ascii -M $BASE.tst-$ISTRING $BASE.tst-$ISTRING.svm
done

for ((I=0; I<NTRIAL; I++)); do
  ISTRING=$(printf "%2.2d" $I)
  rm -f $BASE.trn-$ISTRING.svm 
  for ((J=0; J<NTRIAL; J++)); do
    if [ $J -ne $I ]; then
      echo "cat $BASE.tst-$(printf "%2.2d" $J).svm >> $BASE.trn-$ISTRING.svm"
      cat $BASE.tst-$(printf "%2.2d" $J).svm >> $BASE.trn-$ISTRING.svm
    fi
  done
  echo "svm2agf $BASE.trn-$ISTRING.svm $BASE.trn-$ISTRING"
  svm2agf $BASE.trn-$ISTRING.svm $BASE.trn-$ISTRING
done

rm -f svm_train.log
rm -f svm_test.log

rm -f border_train.log
rm -f border_test.log

for ((I=0; I<NTRIAL; I++)); do
  ISTRING=$(printf "%2.2d" $I)
  TRAINBASE=$BASE.trn-$ISTRING
  TESTBASE=$BASE.tst-$ISTRING
  OUTBASE=$BASE.$ISTRING

  #Step 1: make a fully native LIBSVM model to:
  #	a) check the time
  #	b) compare against hybrid model
  echo "(time svm-train -b 1 $SVMOPT $TRAINBASE.svm $OUTBASE.svmmod) 2>> svm_train.log"
  (time svm-train -b 1 $SVMOPT $TRAINBASE.svm $OUTBASE.svmmod) 2>> svm_train.log
  echo "(time svm-predict -b 1 $TESTBASE.svm $OUTBASE.svmmod $OUTBASE.0.svmout) 2>> svm_test.log"
  (time svm-predict -b 1 $TESTBASE.svm $OUTBASE.svmmod $OUTBASE.0.svmout) 2>> svm_test.log
  echo "svmout2agf $OUTBASE.0.svmout $BASE.libsvm.$ISTRING >> $BASE.junk.txt"
  svmout2agf $OUTBASE.0.svmout $BASE.libsvm.$ISTRING >> $BASE.junk.txt

  #Step 2: build the hybrid borders/libsvm model:
  echo "multi_borders -M -- \"svm-train -b 1 $SVMOPT\" $CONTROL $TRAINBASE.svm $BASE $OUTBASE.svc"
  multi_borders -M -- "svm-train -b 1 $SVMOPT" $CONTROL $TRAINBASE.svm $OUTBASE $OUTBASE.svc
  echo "classify_m -M -O \"svm-predict -b 1\" -Q $METHOD $OUTBASE.svc $TESTBASE.svm $OUTBASE.svm >> $BASE.txt"
  classify_m -M -O "svm-predict -b 1" -Q $METHOD $OUTBASE.svc $TESTBASE.svm $OUTBASE.svmout >> $BASE.txt
  echo "svmout2agf $OUTBASE.svmout $BASE.svm.$ISTRING >> $BASE.txt"
  svmout2agf $OUTBASE.svmout $BASE.svm.$ISTRING >> $BASE.txt

  #Step 3: accelerate the hybrid model:
  echo "(time (multi_borders -Z -s $NSAMPLE $OUTBASE.svc $TRAINBASE $OUTBASE $OUTBASE.mbc 2>> $BASE.jnk)) 2>> border_train.log"
  (time (multi_borders -Z -s $NSAMPLE $OUTBASE.svc $TRAINBASE $OUTBASE $OUTBASE.mbc 2>> $BASE.jnk)) 2>> border_train.log
  echo "(time (classify_m -Q $METHOD $OUTBASE.mbc $TESTBASE.vec $OUTBASE &>> $BASE.txt)) 2>> border_test.log"
  (time (classify_m -Q $METHOD $OUTBASE.mbc $TESTBASE.vec $OUTBASE &>> $BASE.txt)) 2>> border_test.log
done

echo "cat $BASE.libsvm.??.cls > $RESULTBASE.libsvm.cls"
cat $BASE.libsvm.??.cls > $RESULTBASE.libsvm.cls
echo "cat $BASE.libsvm.??.con > $RESULTBASE.libsvm.con"
cat $BASE.libsvm.??.con > $RESULTBASE.libsvm.con

echo "cat $BASE.svm.??.cls > $RESULTBASE.svm.cls"
cat $BASE.svm.??.cls > $RESULTBASE.svm.cls
echo "cat $BASE.svm.??.con > $RESULTBASE.svm.con"
cat $BASE.svm.??.con > $RESULTBASE.svm.con

echo "cat $BASE.??.cls > $RESULTBASE.cls"
cat $BASE.??.cls > $RESULTBASE.out.cls
echo "cat $BASE.??.con > $RESULTBASE.con"
cat $BASE.??.con > $RESULTBASE.con

cat $BASE.tst-??.cls > $RESULTBASE.truth.cls

echo "cls_comp_stats $RESULTBASE.libsvm.cls $RESULTBASE.svm"
cls_comp_stats $RESULTBASE.libsvm.cls $RESULTBASE.svm
agf_correlate $RESULTBASE.libsvm.con $RESULTBASE.svm.con

echo "cls_comp_stats truth.$BASE.cls $RESULTBASE.svm"
cls_comp_stats truth.$BASE.cls $RESULTBASE.svm

echo "cls_comp_stats truth.$BASE.cls $RESULTBASE.out"
cls_comp_stats truth.$BASE.cls $RESULTBASE.out

#clean up, remove temporary files:
if [ -z $KFLAG ]; then
  echo "rm -f $BASE*"
  rm -f $BASE*
fi
  
