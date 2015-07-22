#!/bin/bash

NGRID=100
NV=5		#number of eigenvectors/which one to plot
N=100
SKIP=1
while getopts 'n:v:S:' ARG; do
  case $ARG in
    n) NGRID=$OPTARG
       ;;
    v) NV=$OPTARG
       ;;
    S) SKIP=$OPTARG
       ;;
  esac
done

shift $((OPTIND-1))


BASE=$1
OUTFILE=$2
N=$3

I0=$SKIP
ZFILE=contours.cpt

for ((INDEX=$I0; INDEX<$N;INDEX+=$SKIP)); do
  FRAME=$BASE$INDEX
  plot_frame.sh -z 100 -n $NGRID $FRAME.dat $[NV-1] $FRAME.ps
  convert $FRAME.ps $FRAME.gif
  LIST="$LIST $FRAME.gif"
done

convert -dispose previous -delay 10 $LIST $OUTFILE

