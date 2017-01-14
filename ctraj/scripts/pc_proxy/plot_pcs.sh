#!/bin/bash

VFIELDN=vfield_ncep_500Koneday1997.01.01-1999.12.31_N.ds
VFIELDS=vfield_ncep_500Koneday1997.01.01-1999.12.31_S.ds

NEV=5

while getopts 'Hv:A:lh:k:n:C:G:' ARG; do
  case $ARG in
    v) SVDOPTS="$SVDOPTS -v $OPTARG"
       NEV=$OPTARG
      ;;
    A) SVDOPTS="$SVDOPTS -A $OPTARG"
       WCFLAG=1
      ;;
    l) SVDOPTS="$SVDOPTS -l";
      ;;
    h) TROPTS="$TROPTS -h $OPTARG"
      ;;
    k) TROPTS="$TROPTS -k $OPTARG"
      ;;
    n) TROPTS="$TROPTS -n $OPTARG"
      ;;
    C) TROPTS="$TROPTS -C $OPTARG"
       DSWITCH="-V 3"
      ;;
    G) TROPTS="$TROPTS -G $OPTARG"
       DSWITCH="-V 3"
      ;;
    H) echo "Usage: plot_pcs.sh [-v nv] [-A nA] [-l]"
       echo "            start_date integration_time outbase"
       exit -1
       ;;
  esac
done

shift $((OPTIND-1))

START=$1
INT=$2
BASE=$3

MATFILE=mat$RANDOM.dat
END=$(date_calc $START+0/0/1*$INT)

ctraj_tracer $DSWITCH $TROPTS -i $START -f $END $VFIELDS $VFIELDN $MATFILE

svd_sparse_array $SVDOPTS $MATFILE $BASE.vec

for ((I=0; I<$NEV; I++)); do
  plot_frame.sh -q $BASE.vec $I $BASE$I.ps
done

rm -f $MATFILE
rm -f $BASE.vec

