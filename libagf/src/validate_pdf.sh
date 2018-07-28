#!/bin/bash

ALGTYPE=6
NDIV=5
NSAMP=100

#executable names:
SIM_COM=pdf_sim
KNN_COM=pdf_knn
AGF_COM=pdf_agf

TRAIN_BASE=train_pdf
TEST_BASE=test_pdf
SIM_BASE=pdf_simdata

while getopts "c:d:k:K:s:v:V:W:ni:l:I:h" ARG
do
  case $ARG in
    c) ALGTYPE=$OPTARG
       ;;
    d) NDIV=$OPTARG
       ;;
    k) KOPT="-k $OPTARG"
       ;;
    K) AGFOPTS="$AGFOPTS -k $OPTARG"
       ;;
    s) NSAMP=$OPTARG
       ;;
    v) AGFOPTS="$AGFOPTS -v $OPTARG"
       ;;
    V) AGFOPTS="$AGFOPTS -V $OPTARG"
       ;;
    W) AGFOPTS="$AGFOPTS -W $OPTARG"
       ;;
    n) NFLAG="-n"
       ;;
    i) AGFOPTS="$AGFOPTS -i $OPTARG"
       ;;
    I) AGFOPTS="$AGFOPTS -I $OPTARG"
       ;;
    l) AGFOPTS="$AGFOPTS -l $OPTARG"
       ;;
    h) echo "Validates pdf estimates by testing two or more groups of training data"
       echo "at the same test points.  Prints out a cross-correlation matrix."
       echo
       echo "syntax: validate_pdf [options] {train | train1 train2 [train3 [...]] test }"
       echo
       echo "where:"
       echo "  train   is a single file containing vector training data"
       echo "            this data is randomly split up into two or more subsets and"
       echo "            tested with a set of artificially generated test data"
       echo "       or"
       echo "  train1, train2... is a list of files containing training data"
       echo "            representing the same underlying pdf"
       echo "      and"
       echo "  test    is a set of locations at which estimates will be compared"
       echo
       echo "options:"
       echo "  -c      is the type of algorithm to use (5=AGF, 6=KNN)"
       echo "  -d      is the number of divisions into which to split the data"
       echo "  -s      is the number of samples to use in the test data"
       echo "  -K      use -K for k-nearest-neighbours in AGF routines"
       ${AGF_COM} | tail -n 9 
       ;;
  esac
done

shift $(($OPTIND - 1))

NTEST=$#

NTEST=$((NTEST - 1))

if [ $NTEST -eq 0 ]
then
  NTEST=$NDIV
  echo $NTEST
  TRAINBASE=${TRAIN_BASE}-$(date +%N)
  echo "agf_preprocess -C -d $NDIV $1 $TRAINBASE"
  agf_preprocess -C -d $NDIV $1 $TRAINBASE
  TESTFILE=${SIM_BASE}-$(date +%N).vec
  DELTEST=1
  echo "${SIM_COM} -s $NSAMP $AGFOPTS $1 $TESTFILE"
  ${SIM_COM} -s $NSAMP $AGFOPTS $1 $TESTFILE
  for ((I=1; I<=$NTEST; I++))
  do
    TRAINFILE[$I]=$TRAINBASE-$(printf "%2.2d" $((I-1))).vec
    OUTFILE[$I]=${TEST_BASE}-$(date +%N)-$I.dum
    shift
  done
else
  for ((I=1; I<=$NTEST; I++))
  do
    TRAINFILE[$I]=$1
    OUTFILE[$I]=${TEST_BASE}-$(date +%N)-$I.dum
    shift
  done
  TESTFILE=$1
fi

if [ $ALGTYPE -eq 6 ]
then
  COMMAND="${AGF_COM} $NFLAG $AGFOPTS"
elif [ $ALGTYPE -eq 7 ]
then
  COMMAND="${KNN_COM} $NFLAG $KOPT"
else
  COMMAND="${AGF_COM} $NFLAG $AGFOPTS"
fi

for ((I=1; I<=$NTEST; I++))
do
  echo "$COMMAND ${TRAINFILE[$I]} $TESTFILE ${OUTFILE[$I]}"
  $COMMAND ${TRAINFILE[$I]} $TESTFILE ${OUTFILE[$I]}
done

for ((I=1; I<=$NTEST; I++)); do
  for ((J=1; J<=$NTEST; J++)); do
    R=$(agf_correlate ${OUTFILE[$I]} ${OUTFILE[$J]})
    echo -n "$R "
  done
  echo
done

if [ $DELTEST ]
then
  rm $TESTFILE
  for ((I=1; I<=$NTEST; I++))
  do
    rm ${TRAINFILE[$I]}
  done
fi

for ((I=1; I<=$NTEST; I++))
do
  rm ${OUTFILE[$I]}
done

