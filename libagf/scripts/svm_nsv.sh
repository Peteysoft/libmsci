#!/bin/sh

echo on

SCDIR=/home/lenovo/my_software/libmsci/libagf/examples/sample_classes

N="10 20 30 40 50 75 100 200 300 400 500 750 1000 2000 3000 4000 5000 7500 10000"
#N="10 50 100 500 1000 5000 10000"

#$BASE/sample_class -R 1000 2000 sctest > sctest.lvq
#agf2ascii -M sctest sctest.svm

#numbers:
NTRIAL=10
NTEST=10
MINSAMPLE=10
MAXSAMPLE=100

NSCTEST=5000

while getopts 'GKMq:s:c:d:g:h:r:t:S:N:b:' DUM
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
    s) MINSAMPLE=$OPTARG
      ;;
    S) MAXSAMPLE=$OPTARG
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
    b) BFLAG="-b 1"
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
elif [ $# -eq 1 ]; then
  OUTFILE=$1
else
  echo "   SUPPORT-VECTOR-MACHINE NUMBER OF SUPPORT VECTORS"
  echo
  echo "svm_nsv.sh [-G] [-M] [-K] [-N ntrial] [-q ntest] [-s min] [-S max]\\"
  echo "        [options] [train] output"
  echo
  echo "  options   = options to pass to svm-train"
  echo "  train     = training data--if omitted uses synthetic test classes"
  echo "  output    = output stats"
  echo
  echo "  -G         = use logarithmic progression"
  echo "  -K         = keep temporary files"
  echo "  -M         = data is in LIBSVM format"
  echo
  echo "  ntrial    = number of repeated trials [$NTRIAL]"
  echo "  ntest     = number of sample sizes [$NTEST]"
  echo "  frac      = fraction to reserve for testing [$FRAC]"
  echo "  min       = minimum number of samples [$MINSAMPLE]"
  echo "  max       = maximum number of samples [$MAXSAMPLE]"
  exit 0
fi

BASE=sv$RANDOM.tmp

if [ -z ${TRAINING_DATA} ]; then
  if test $SVMFLAG; then
    #this is kind of dum:
    echo "svm2agf ${TRAINING_DATA} $BASE.0"
    svm2agf ${TRAINING_DATA} $BASE.0
    TRAINING_DATA=$BASE.0
  fi
fi

rm -f $OUTFILE
rm -f test.log

for ((I=0; I<NTEST; I++)); do
  if [ $LOGFLAG ]; then
    n=$(echo "a=e(l($MINSAMPLE)+$I*(l($MAXSAMPLE)-l($MINSAMPLE))/($NTEST-1)); scale=0; (a+1)/1" | bc -l)
  else
    n=$((MINSAMPLE+I*(MAXSAMPLE-MINSAMPLE)/(NTEST-1)))
  fi

  echo $n
  f=$(echo "scale=3; $n/$MAXSAMPLE" | bc -l)

  for ((J=0; J<NTRIAL; J++)); do
    if [ -z ${TRAINING_DATA} ]; then
      $SCDIR/sample_class $((n/3)) $((2*n/3)) $BASE.trn > $BASE.trn.lvq
      $SCDIR/sample_class -R $((NSCTEST/3)) $((2*NSCTEST/3)) $BASE.tst > $BASE.trn.lvq
    else
      echo "agf_preprocess -zf $f ${TRAINING_DATA} $BASE.trn $BASE.tst"
      agf_preprocess -zf $f ${TRAINING_DATA} $BASE.trn $BASE.tst
      echo "agf2ascii -M $BASE.trn $BASE.trn.svm"
      echo "agf2ascii -M $BASE.tst $BASE.tst.svm"
    fi
    agf2ascii -M $BASE.trn $BASE.trn.svm
    agf2ascii -M $BASE.tst $BASE.tst.svm

    echo "svm-train $BFLAG $SVMOPT $BASE.tst.svm $BASE.svmmod"
    svm-train $BFLAG $SVMOPT $BASE.trn.svm $BASE.svmmod
    echo "svm-predict $BFLAG $BASE.tst.svm $BASE.svmmod $BASE.svmout"
    time -o $BASE.tm svm-predict $BFLAG $BASE.tst.svm $BASE.svmmod $BASE.svmout

    echo -n "$(wc -l < $BASE.trn.svm) " >> $OUTFILE
    grep nr_sv $BASE.svmmod | tr -d '\n' >> $OUTFILE
    echo -n " " >> $OUTFILE
    grep -o "[0-9]*\.[0-9]*" $BASE.tm | sed '1q;d' | tr -d '\n' >> $OUTFILE
    echo -n " " >> $OUTFILE
    grep -o "[0-9]*\.[0-9]*" $BASE.tm | sed '2q;d' >> $OUTFILE
  done
done

#clean up, remove temporary files:
if [ -z $KFLAG ]; then
  rm -f $BASE.*
fi
 
