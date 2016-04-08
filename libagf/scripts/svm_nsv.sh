#!/bin/bash

echo on

#BASE=/home/lenovo/my_software/libmsci/libagf/examples/sample_classes

N="10 20 30 40 50 75 100 200 300 400 500 750 1000 2000 3000 4000 5000 7500 10000"
#N="10 50 100 500 1000 5000 10000"

#$BASE/sample_class -R 1000 2000 sctest > sctest.lvq
#agf2ascii -M sctest sctest.svm

#numbers:
NTRIAL=1
NTEST=10
MAXSAMPLE=100

while getopts 'GKMq:s:c:d:g:h:r:t:' DUM
do
  case $DUM in
    G) LOGFLAG=1
      ;;
    K) KFLAG=1
      ;;
    M) SVMFLAG=1
      ;;
    N) NTRIAL=$OPTARG
      ;;
    q) NTEST=$OPTARG
      ;;
    s) MAXSAMPLE=$OPTARG
      ;;
    c) SVMOPT="$SVMOPT -c $OPTARG"
      ;;
    d) SVMOPT="$SVMOPT -d $OPTARG"
      ;;
    g) SVMOPT="$SVMOPT -g $OPTARG"
      ;;
    h) SVMOPT="$SVMOPT -h $OPTARG"
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
  TRAINING_DATA=$1
  OUTFILE=$2
else
  echo "   SUPPORT-VECTOR-MACHINE NUMBER OF SUPPORT VECTORS"
  echo
  echo "svm_nsv.sh [-g] [-M] [-K] [-N ntrial] [-q ntest]\\"
  echo "        [options] train output"
  echo
  echo "  options   = options to pass to svm-train"
  echo "  train     = training data"
  echo "  output    = output stats"
  echo
  echo "  G         = use logarithmic progression"
  echo "  K         = keep temporary files"
  echo "  M         = data is in LIBSVM format"
  echo
  echo "  ntrial    = number of repeated trials [$NTRIAL]"
  echo "  ntest     = number of sample sizes [$NTEST]"
  echo "  frac      = fraction to reserve for testing [$FRAC]"
  exit 0
fi

BASE=sv$RANDOM.tmp

if test $SVMFLAG; then
  #this is kind of dum:
  echo "svm2agf ${TRAINING_DATA} $BASE.0"
  svm2agf ${TRAINING_DATA} $BASE.0
  NTRAIN=$(wc -l < ${TRAINING_DATA})
  TRAINING_DATA=$BASE.0
else
  echo "agf2ascii -M ${TRAINING_DATA} $BASE.0.svm"
  agf2ascii -M ${TRAINING_DATA} $BASE.0.svm
  NTRAIN=$(wc -l < $BASE.0.svm)
fi


rm -f $OUTFILE

for ((I=0; I<NTEST; I++)); do
  if [ $LOGFLAG ]; then
    n=$(echo "e(($I+1)*l($MAXSAMPLE)/$NTEST)" | bc -l)
  else
    n=$(((I+1)*MAXSAMPLE/NTEST))
  fi

  f=$(echo "$n/$MAXSAMPLE" | bc -l)

  #$BASE/sample_class $n $((n*2)) sctrain > sctrain.lvq
  #agf2ascii -M sctrain sctrain.svm

  echo "agf_preprocess -zf $f ${TRAINING_DATA} $BASE.trn $BASE.tst"
  agf_preprocess -zf $f ${TRAINING_DATA} $BASE.trn $BASE.tst
  echo "agf2ascii -M $BASE.trn $BASE.trn.svm"
  agf2ascii -M $BASE.trn $BASE.trn.svm
  echo "agf2ascii -M $BASE.tst $BASE.tst.svm"
  agf2ascii -M $BASE.tst $BASE.tst.svm

  echo "svm-train $SVMOPT $BASE.tst.svm $BASE.svmmod"
  svm-train $SVMOPT $BASE.tst.svm $BASE.svmmod
  echo "svm-predict $BASE.trn.svm $BASE.svmmod $BASE.svmout"
  svm-predict $BASE.trn.svm $BASE.svmmod $BASE.svmout

  echo -n "$(wc -l < $BASE.tst.svm) " >> $OUTFILE
  grep nr_sv $BASE.svmmod >> $OUTFILE
done

#clean up, remove temporary files:
if [ -z $KFLAG ]; then
  rm -f $BASE.*
fi
 
