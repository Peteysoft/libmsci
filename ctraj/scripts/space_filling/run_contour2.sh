#!/bin/bash

set -e

ZOPTS=""
WRITEINT=1
I0=0
N=10

while getopts '0:N:D:M:c:m:s:B:h:k:O:HoZ' ARG; do
  case $ARG in
    Z) BG=1
      ;;
    0) I0=$OPTARG
      ;;
    N) N=$OPTARG
      ;;
    D) CONFIRM=$OPTARG
       DOPT="-D $CONFIRM"
      ;;
    M) ZOPTS="$ZOPTS -M $OPTARG"
        MAXNPT=$OPTARG
      ;;
    c) ZOPTS="$ZOPTS -c $OPTARG"
      ;;
    m) ZOPTS="$ZOPTS -m $OPTARG"
      ;;
    s) ZOPTS="$ZOPTS -s $OPTARG"
      ;;
    B) ZOPTS="$ZOPTS -B $OPTARG"
      ;;
    h) ZOPTS="$ZOPTS -h $OPTARG"
        H=$OPTARG
      ;;
    k) ZOPTS="$ZOPTS -k $OPTARG"
      ;;
    o) OFLAG="-o"
      ;;
    O) WRITEINT=$OPTARG
       ZOPTS="$ZOPTS -O $OPTARG"
      ;;
    H) echo "Usage: run_contour2 [options] nfile sfile cbase"
       echo "where:"
       echo "  nfile   is N. hemi. v-field"
       echo "  sfile   is S. hemi. v-field"
       echo "  cbase   base-name for output contours"
       echo "options:"
       echo "  -Z   run in background (**warning: recursive algorithm)"
       ctraj_contour -? | tail -n +17 | head -n 6
       echo "  -N   number of time steps [$N]"
       ctraj_contour -? | tail -n +26 | head -n 3
       exit 0
       ;;
  esac
done

shift $(($OPTIND - 1))

nfile=$1
sfile=$2
outbase=$3

if [ $[I0 % WRITEINT] -ne 0 ]
then
  echo "Please use a multiple of write interval ($WRITEINT)"
  echo "for the starting index ($I0)"
  exit 1
fi

cp $outbase.bev $outbase.$[I0/WRITEINT].bev

echo "ctraj_contour -N $N $OFLAG $ZOPTS -0 $I0 $nfile $sfile $cfile $outbase.$[I0/WRITEINT].bev"
ctraj_contour -N $N $OFLAG $ZOPTS -0 $I0 $nfile $sfile $outbase.$[I0/WRITEINT].bev

NCOM=$[$(bev2xy $outbase.$[I0/WRITEINT].bev | wc -l)-1]

I0_N=$[I0+NCOM*WRITEINT]
NN=$[N-NCOM*WRITEINT]

if [ $NN -le 0 ]
then
  if test -e $CONFIRM
  then
    if [ $(pgrep -f contour2.*$outbase | wc -l) -eq 1 ]
    then
      echo run_contour2 completed > $CONFIRM
    fi
  fi
  exit 0
fi

newstart=${outbase}_$[I0_N/WRITEINT].txt

bev2xy $outbase.$[I0/WRITEINT].bev $[NCOM] > $newstart

NPT=$(wc -l < $newstart)
echo $NPT

tail -n 1 $newstart > ${outbase}0_0.txt
head -n $[NPT/2+1] ${newstart} >> ${outbase}0_0.txt
xy2bev ${outbase}0.bev < ${outbase}0_0.txt

tail -n $[NPT/2+1] $newstart > ${outbase}1_0.txt
head -n 1 $newstart >> ${outbase}1_0.txt
xy2bev ${outbase}1.bev < ${outbase}1_0.txt

if [[ $BG ]]; then
  echo "./run_contour2.sh $DOPT -o $ZOPTS $nfile $sfile ${outbase}0 $I0_N $NN &"
  echo "./run_contour2.sh $DOPT -o $ZOPTS $nfile $sfile ${outbase}1 $I0_N $NN &"
  ./run_contour2.sh -o $ZOPTS $nfile $sfile ${outbase}0 $I0_N $NN &
  ./run_contour2.sh -o $ZOPTS $nfile $sfile ${outbase}1 $I0_N $NN &
else
  echo "./run_contour2.sh $DOPT -o $ZOPTS $nfile $sfile ${outbase}0 $I0_N $NN"
  echo "./run_contour2.sh $DOPT -o $ZOPTS $nfile $sfile ${outbase}1 $I0_N $NN"
  ./run_contour2.sh -o $ZOPTS $nfile $sfile ${outbase}0 $I0_N $NN
  ./run_contour2.sh -o $ZOPTS $nfile $sfile ${outbase}1 $I0_N $NN
fi

if test -e $CONFIRM
then
  if [ $(pgrep -f contour2.*$outbase | wc -l) -eq 1 ]
  then
    echo run_contour2 completed > $CONFIRM
  fi
fi

