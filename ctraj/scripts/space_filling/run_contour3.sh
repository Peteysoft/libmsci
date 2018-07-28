#!/bin/bash

ZOPTS=""

while getopts 'm:s:B:h:k:N:H' ARG; do
  case $ARG in
    m) ZOPTS="$ZOPTS -m $OPTARG"
      ;;
    s) ZOPTS="$ZOPTS -s $OPTARG"
      ;;
    B) ZOPTS="$ZOPTS -B $OPTARG"
      ;;
    h) ZOPTS="$ZOPTS -h $OPTARG"
      ;;
    k) ZOPTS="$ZOPTS -k $OPTARG"
      ;;
    N) ZOPTS="$ZOPTS -N $OPTARG"
      H=$OPTARG
      ;;
    H) echo "Usage: run_contour3 [options] nfile sfile cfile I0 N outbase "
       ;;
  esac
done

shift $(($OPTIND - 1))

nfile=$1
sfile=$2
cfile=$3
I0=$4
N=$5
outbase=$6

echo "contour3 $ZOPTS -0 $I0 $nfile $sfile $cfile $outbase.$I0.vec"
contour3 $ZOPTS -0 $I0 $nfile $sfile $cfile $outbase.$[I0+H].vec
for ((INDEX=$[I0+H]; INDEX<=$[N*H]; INDEX+=$H)); do
  outfile=$outbase.$[INDEX+H].vec
  contour3 $ZOPTS -0 $INDEX $nfile $sfile $cfile $outfile
done

