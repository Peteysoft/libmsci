#!/bin/bash

set -e

ZOPTS=""

while getopts 'H' ARG; do
  case $ARG in
    H) echo "Usage: contour_growth base i0 n"
       exit 0
       ;;
  esac
done

shift $(($OPTIND - 1))

fbase=$1
I0=$2
N=$3
STRIDE=${4:-1}

for ((INDEX=$I0; INDEX<=$N; INDEX+=$STRIDE)); do
  flist=$(ls $fbase*.$[INDEX/STRIDE].vec)
  npt1=($(du -c -B 8 $flist | tail -n 1))
  npt2=${npt1[0]}
  flista=($flist)
  nf=${#flista[@]}
  #header is only 4 bytes, but points on intersections are repeated 
  # (sometimes two points/intersection...)
  # -> npt2 is total size of files in 8 byte blocks,
  #    3 repeated points * 4 bytes + 4 byte header = 16 bytes/8 = 2
  npt=$[npt2-2*nf]
  s=$(contour_length $flist)
  echo $[INDEX/STRIDE] $s $npt
done

