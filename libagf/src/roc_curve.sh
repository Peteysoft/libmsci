#!/bin/bash

NTRIAL=20

while getopts "0:1:c:d:f:h:k:i:I:l:q:s:t:v:V:W:gHn" ARG
do
  case $ARG in
    0) MIN=$OPTARG
       ;;
    1) MAX=$OPTARG
       ;;
    c) AGFOPTS="$AGFOPTS -c $OPTARG"
       ;;
    d) AGFOPTS="$AGFOPTS -f $((1/$OPTARG))"
       ;;
    f) AGFOPTS="$AGFOPTS -f $OPTARG"
       ;;
    g) GFLAG=1
       ;;
    k) AGFOPTS="$AGFOPTS -k $OPTARG"
       ;;
    q) NTRIAL=$OPTARG
       ;;
    r) RFLAG=1
       ;;
    s) AGFOPTS="$AGFOPTS -s $OPTARG"
       ;;
    t) AGFOPTS="$AGFOPTS -t $OPTARG"
       ;;
    v) AGFOPTS="$AGFOPTS -v $OPTARG"
       ;;
    V) AGFOPTS="$AGFOPTS -V $OPTARG"
       ;;
    W) AGFOPTS="$AGFOPTS -W $OPTARG"
       ;;
    n) AFGOPTS="$AGFOPTS -n"
       ;;
    h) AGFOPTS="$AGFOPTS -h $OPTARG"
       ;;
    i) AGFOPTS="$AGFOPTS -i $OPTARG"
       ;;
    I) AGFOPTS="$AGFOPTS -I $OPTARG"
       ;;
    l) AGFOPTS="$AGFOPTS -l $OPTARG"
       ;;
    H) echo "Generates the ROC-curve in two possible ways:"
       echo "1. by converting floating-point data"
       echo "2. by shifting the discrimination border [-r option]"
       echo
       echo "syntax: roc_curve [-r] [-c type] [-q n] [-g] train"
       echo
       echo "where:"
       echo "  train   is the base-name for the training data files"
       echo
       echo "options:"
       echo "  -c  type = is the type of algorithm to use (0=AGF borders, 1=AGF, 2=KNN)"
       echo "  -q     n = is the number trials"
       echo "  -g       = use a geometric progression"
       nfold | tail -n 15 | head -n 12
       ;;
  esac
done

shift $(($OPTIND - 1))

if [ $RFLAG ]
then
  if [ -z $MIN ] 
  then
    MIN=-1
  fi
  if [ -z $MAX ] 
  then
    MAX=1
  fi
else
  echo float_to_class -M $1.dat
  FDATA=$(float_to_class -M $1.dat)
  I=0
  for VAL in $FDATA
  do
    RESULT[$I]=$VAL
    I=$((I+1))
  done
  if [ -z $MIN ] 
  then
    MIN=${RESULT[0]}
  fi
  if [ -z $MAX ] 
  then
    MAX=${RESULT[1]}
  fi

  echo $MIN
  echo $MAX
fi

if [ $GFLAG ]
then
  for ((I=1; I<$(($NTRIAL-1)); I++)); do
    R[$I]=$(($MIN*($MAX/$MIN)**($I/($NTRIAL-1))))
  done
else
  for ((I=1; I<$(($NTRIAL-1)); I++)); do
    R[$I]=$(($MIN+($MAX-$MIN)*$I/($NTRIAL-1)))
  done
fi

for ((I=0; I<$NTRIAL; I++)); do
  if [ $RFLAG ]
  then
    ROPT="-r ${R[$I]}"
  else
    echo float_to_class $1.dat ${R[$I]} $1.cls
    float_to_class $1.dat ${R[$I]} $1.cls
  fi

  FDATA=$(nfold $AGFOPTS $ROPT $1 | tail -n 23 | head -n 2)
  I=0
  for VAL in FDATA
  do
    RESULT[$I]=$VAL
    I=$((I+1))
  done

  echo $((${RESULT[3]}/${RESULT[2]}+${RESULT[3]})) \
		$((${RESULT[7]}/${RESULT[6]}+${RESULT[7]})) \
		
done


