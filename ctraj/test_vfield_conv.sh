#!/bin/bash

set -e

DATAPATH=.

while getopts 'n:r:p:B:Tw' ARG; do
  case $ARG in
    n)  OPTS="$OPTS -n $OPTARG"
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
    H) echo "test_vfield_conv.sh [-p path] [-T] level t0 tf ntrial"
      ;;
  esac
done

shift $(($OPTIND -1))

ZLEV=$1
T0=$2
TF=$3
NTRIAL=$4

#ID=$RANDOM
ID=0

NFILE=vfieldN.$ID.ds
SFILE=vfieldS.$ID.ds
UWND=uwnd.$ID.txt
VWND=vwnd.$ID.txt
WWND=wwnd.$ID.txt
RESULTS=results.$ID.txt

echo "nc2ds -p $DATAPATH $OPTS -i $T0 -f $TF $TFLAG $ZLEV $SFILE $NFILE"
nc2ds -p $DATAPATH $OPTS -i $T0 -f $TF $TFLAG $ZLEV $SFILE $NFILE

for ((I=0; I<$NTRIAL; I++))
do
  RANNO=$RANDOM
  echo "DATE=\$(date_calc \"(${TF}_${T0})*$RANNO|32767+$T0\")"
  DATE=$(date_calc "(${TF}_${T0})*$RANNO|32767+$T0")

  echo "extract_ncep_field -L $TFLAG -p $DATAPATH uwnd $ZLEV $DATE > $UWND"
  extract_ncep_field -L $TFLAG -p $DATAPATH uwnd $ZLEV $DATE > $UWND
  echo "extract_ncep_field -L $TFLAG -p $DATAPATH vwnd $ZLEV $DATE > $VWND"
  extract_ncep_field -L $TFLAG -p $DATAPATH vwnd $ZLEV $DATE > $VWND

  echo "vfield_interpolate u $SFILE $NFILE $UWND | tail -n 1 >> $RESULTS"
  vfield_interpolate -P u $SFILE $NFILE $UWND | tail -n 1 >> $RESULTS
  echo "vfield_interpolate v $SFILE $NFILE $VWND | tail -n 1 >> $RESULTS"
  vfield_interpolate -P v $SFILE $NFILE $VWND | tail -n 1 >> $RESULTS
  if test $WFLAG
  then
    extract_ncep_field -L $TFLAG -p $DATAPATH wwnd $ZLEV $DATE > $WWND
    vfield_interpolate w $SFILE $NFILE $WWND | tail -n 1 >> $RESULTS
  fi
done

more $RESULTS

rm -f $RESULTS
rm -f $UWND
rm -f $VWND
if test $WFLAG
then
  rm -f $WWND
fi
rm -f $SFILE $NFILE

