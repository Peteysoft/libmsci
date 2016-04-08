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

#output file:
OUTFILE=stats.txt

echo $#

if [ $# -eq 2 ]; then
  CONTROL=$1
  TRAINING_DATA=$2
else
  echo "   N-FOLD FOR ACCELERATED SVM"
  echo
  echo "nfold_accel.sh [-M] [-K] [-s nsamples] [-d ntrial] [-Q method]\\"
  echo "           [svm-options] control train"
  echo
  echo "  svm-options  = options to pass to svm-train"
  echo "  control    = input control file"
  echo "  train     = training data"
  echo
  echo "  K         = keep temporary files"
  echo "  M         = data is in LIBSVM format"
  echo
  echo "  nsamples  = number border samples [$NSAMPLE]"
  echo "  ntrial    = number of folds [$NTRIAL]"
  echo "  method    = method of solving for conditional probabilities [$METHOD]"
  exit 0
fi

BASE=accel$RANDOM.tmp

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

for ((I=0; I<NTRIAL; I++)); do
  ISTRING=$(printf "%2.2d" $I)

  echo "multi_borders -M -- \"svm-train -b 1 $SVMOPT\" $CONTROL $BASE.trn-$ISTRING.svm $BASE $BASE.svc"
  multi_borders -M -- "svm-train -b 1 $SVMOPT" $CONTROL $BASE.trn-$ISTRING.svm $BASE $BASE.svc
  echo "classify_m -M -O \"svm-predict -b 1\" -Q $METHOD $BASE.svc $BASE.tst-$ISTRING.svm $BASE.$ISTRING.svm >> $BASE.txt"
  classify_m -M -O "svm-predict -b 1" -Q $METHOD $BASE.svc $BASE.tst-$ISTRING.svm $BASE.$ISTRING.svmout >> $BASE.txt
  echo "svmout2agf $BASE.$ISTRING.svmout $BASE.svm.$ISTRING >> $BASE.txt"
  svmout2agf $BASE.$ISTRING.svmout $BASE.svm.$ISTRING >> $BASE.txt

  echo "multi_borders -Z -s $NSAMPLE $BASE.svc $BASE.trn-$ISTRING $BASE $BASE.mbc"
  multi_borders -Z -s $NSAMPLE $BASE.svc $BASE.trn-$ISTRING $BASE $BASE.mbc
  echo "classify_m -Q $METHOD $BASE.mbc $BASE.tst-$ISTRING.vec $BASE.$ISTRING >> $BASE.txt"
  classify_m -Q $METHOD $BASE.mbc $BASE.tst-$ISTRING.vec $BASE.$ISTRING >> $BASE.txt
done

echo "cat $BASE.svm.??.cls > $BASE.svm.cls"
cat $BASE.svm.??.cls > $BASE.svm.cls
echo "cat $BASE.svm.??.con > $BASE.svm.con"
cat $BASE.svm.??.con > $BASE.svm.con

echo "cat $BASE.??.cls > $BASE.out.cls"
cat $BASE.??.cls > $BASE.out.cls
echo "cat $BASE.??.con > $BASE.out.con"
cat $BASE.??.con > $BASE.out.con

cat $BASE.tst-??.cls > $BASE.cls

echo "cls_comp_stats $BASE.cls $BASE.svm"
cls_comp_stats $BASE.cls $BASE.svm

echo "cls_comp_stats $BASE.cls $BASE.out"
cls_comp_stats $BASE.cls $BASE.out

#clean up, remove temporary files:
if [ -z $KFLAG ]; then
  echo "rm -f $BASE*"
  #rm -f $BASE*
fi
  
