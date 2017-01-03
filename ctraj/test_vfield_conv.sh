#!/bin/bash

set -e

DATAPATH=.

while getopts 'i:f:n:r:p:B:TK:wH' ARG; do
  case $ARG in
    i) T0=$OPTARG
      ;;
    f) TF=$OPTARG
      ;;
    K) KFLAG=$OPTARG
      ;;
    n) OPTS="$OPTS -n $OPTARG"
      ;;
    r) OPTS="$OPTS -r $OPTARG"
      ;;
    p) DATAPATH=$OPTARG
      ;;
    B) OPTS="$OPTS -B $OPTARG"
      ;;
    T) TFLAG="-T"
      ;;
    w) WFLAG="-w"
      ;;
    H) echo "test_vfield_conv.sh [-p path] [-T] [level -i T0 -f TF | vfieldS vfieldN] ntrial"
      exit 0 
      ;;
  esac
done

shift $(($OPTIND -1))

ID=$RANDOM

if [ $# -eq 2 ]
then
  ZLEV=$1
  NTRIAL=$2
else
  INFO=info.$ID.txt
  NFILE=$1
  SFILE=$2
  vfield_interpolate $SFILE $NFILE > $INFO
  if test -z $T0
  then
    T0=$(head -n 2 $INFO | tail -n 1)
  fi
  if test -z $TF
  then
    TF=$(head -n 3 $INFO | tail -n 1)
  fi
  ZLEV=$(head -n 5 $INFO | tail -n 1)
  NTRIAL=$3
  if test -z $KFLAG
  then
    KFLAG=0
  fi
fi

#ID=0

UWND=uwnd.$ID.txt
VWND=vwnd.$ID.txt
WWND=wwnd.$ID.txt
RESULTS=results.$ID.txt

if test -z $NFILE
then
  NFILE=vfieldN.$ID.ds
  SFILE=vfieldS.$ID.ds
  echo "nc2ds -p $DATAPATH $OPTS -i $T0 -f $TF $TFLAG $ZLEV $SFILE $NFILE"
  nc2ds -p $DATAPATH $OPTS -i $T0 -f $TF $TFLAG $ZLEV $SFILE $NFILE
fi

for ((I=0; I<$NTRIAL; I++))
do
  RANNO=$RANDOM
  echo "DATE=\$(date_calc \"(${TF}_${T0})*$RANNO|32767+$T0\")"
  DATE=$(date_calc "(${TF}_${T0})*$RANNO|32767+$T0")

  echo "extract_ncep_field -L $TFLAG -p $DATAPATH uwnd $ZLEV $DATE > $UWND"
  extract_ncep_field -L $TFLAG -p $DATAPATH uwnd $ZLEV $DATE > $UWND
  echo "extract_ncep_field -L $TFLAG -p $DATAPATH vwnd $ZLEV $DATE > $VWND"
  extract_ncep_field -L $TFLAG -p $DATAPATH vwnd $ZLEV $DATE > $VWND

  echo $DATE >> $RESULTS
  echo "vfield_interpolate u $SFILE $NFILE $UWND | tail -n 1 >> $RESULTS"
  vfield_interpolate -P u $SFILE $NFILE $UWND | tail -n 1 >> $RESULTS
  echo "vfield_interpolate v $SFILE $NFILE $VWND | tail -n 1 >> $RESULTS"
  vfield_interpolate -P v $SFILE $NFILE $VWND | tail -n 1 >> $RESULTS
  if test $WFLAG
  then
    echo "extract_ncep_field -L $TFLAG -p $DATAPATH wwnd $ZLEV $DATE > $WWND"
    extract_ncep_field -L $TFLAG -p $DATAPATH wwnd $ZLEV $DATE > $WWND
    echo "vfield_interpolate w $SFILE $NFILE $WWND | tail -n 1 >> $RESULTS"
    vfield_interpolate w $SFILE $NFILE $WWND | tail -n 1 >> $RESULTS
  fi
done

more $RESULTS

if test -z $KFLAG
then
  rm -f $RESULTS
  rm -f $UWND
  rm -f $VWND
  if test $WFLAG
  then
    rm -f $WWND
  fi
  rm -f $SFILE $NFILE
  if [ $INFO ]; then
    rm -f $INFO
  fi
elif [ $KFLAG -eq 0 ]
then
  rm -f $RESULTS
  rm -f $UWND
  rm -f $VWND
  if test $WFLAG
  then
    rm -f $WWND
  fi
  if [ $INFO ]; then
    rm -f $INFO
  fi
fi

