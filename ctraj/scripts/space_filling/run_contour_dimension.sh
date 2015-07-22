#!/bin/bash

set -e

ZOPTS=""
STRIDE=1
I0=0
N=-1

while getopts '0:N:O:D:I:F:e:u:U:+-gHZ' ARG; do
  case $ARG in
    Z) BG=1
      ;;
    0) I0=$OPTARG
      ;;
    N) N=$OPTARG
      ;;
    O) STRIDE=$OPTARG
      ;;
    D) ZOPTS="$ZOPTS -D $OPTARG"
      DTYPE=$OPTARG
      ;;
    I) ZOPTS="$ZOPTS -I $OPTARG"
      ;;
    F) ZOPTS="$ZOPTS -F $OPTARG"
      ;;
    e) ZOPTS="$ZOPTS -e $OPTARG"
      ;;
    u) ZOPTS="$ZOPTS -u $OPTARG"
      ;;
    U) ZOPTS="$ZOPTS -U $OPTARG"
      ;;
    +) ZOPTS="$ZOPTS -+"
      ;;
    -) ZOPTS="$ZOPTS --"
      ;;
    g) ZOPTS="$ZOPTS -g"
      ;;
    H) echo "Usage: run_contour_dimension [options] base [outbase]"
       echo "where:"
       echo "  base    is the base file name"
       echo "  outbase is the output base file name"
       echo "options:"
       echo "  -Z   run in background"
       echo "  -0   initial index [0]"
       echo "  -N   number of indices"
       echo "  -O   stride length [1]"
       contour_dimension -?| tail -n +2 | head -n 14
       contour_dimension -?| tail -n 2
       exit 0
       ;;
  esac
done

shift $(($OPTIND - 1))

fbase=$1
if [[ $N -lt 0 ]]; then
  flist=$(ls $fbase*.*.vec)
  for fname in $flist; do
    IND=${fname#$fbase*\.}
    IND=${IND%\.vec}
    #echo $IND
    if [[ $IND -ge $N ]]; then
      N=$[IND+1]
    fi
  done
fi

echo $N
OUTBASE=${2:-$fbase}

for ((INDEX=$I0; INDEX<$[I0+N]; INDEX+=$STRIDE)); do
  echo $INDEX $STRIDE
  flist=$(ls $fbase*.$[INDEX/STRIDE].vec)
  #echo "contour_dimension $ZOPTS $flist > $OUTBASE.$[INDEX/STRIDE].D$DTYPE.txt$BG"
  if [[ $BG ]]; then
    contour_dimension $ZOPTS $flist > $OUTBASE.$[INDEX/STRIDE].D$DTYPE.txt &
  else
    contour_dimension $ZOPTS $flist > $OUTBASE.$[INDEX/STRIDE].D$DTYPE.txt
  fi
done

