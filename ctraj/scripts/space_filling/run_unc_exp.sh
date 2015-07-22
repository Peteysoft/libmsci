#!/bin/bash

ZOPTS=""
H=1

while getopts 'N:I:F:e:u:U:+-gH' ARG; do
  case $ARG in
    N) H=$OPTARG
      ;;
    I) ZOPTS="$ZOPTS -I $OPTARG"
      ;;
    F) ZOPTS="$ZOPTS -F $OPTARG"
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
    H) echo "Usage: run_unc_exp [options] base I0 N"
       exit 0
       ;;
  esac
done

shift $(($OPTIND - 1))

outbase=$1
I0=$2
N=$3

for ((INDEX=$[I0]; INDEX<=$[N*H]; INDEX+=$H)); do
  outfile=$outbase.$[INDEX+H].vec
  unc_exp $ZOPTS $outfile
done

