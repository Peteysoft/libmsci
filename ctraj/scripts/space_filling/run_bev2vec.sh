#!/bin/bash

set -e

while getopts 'Zp:H' ARG; do
  case $ARG in
    Z) BG=Z
       ;;
    p) WORKPATH=$OPTARG
       ;;
    H) echo "run_bev2vec [-Z] [-p workpath] [-H] fbase"
       exit
       ;;
  esac
done

shift $(($OPTIND - 1))

fbase=$1

if test -z $WORKPATH; then
  $WORKPATH=${fbase##*/}
else
  if [ ! -d $WORKPATH ]; then
    mkdir $WORKPATH
  fi
fi


#flist="$fbase.0.bev $(ls $fbase*.*.bev)"
flist="$(ls $fbase*.*.bev)"
echo $flist

fout=${fbase##*/}
fout=$WORKPATH/$fout.0

if [[ $BG ]]; then
  echo "bev2xy $fbase.0.bev 0 | lvq2agf -c -h $fout &"
  bev2xy $fbase.0.bev 0 | lvq2agf -c -h $fout &
else
  echo "bev2xy $fbase.0.bev 0 | lvq2agf -c -h $fout"
  bev2xy $fbase.0.bev 0 | lvq2agf -c -h $fout
fi

for fname in $flist
do
  I0=${fname#$fbase*\.}
  I0=${I0%\.bev}

  NC=$(bev2xy $fname | wc -l)
  echo $fname $I0
  for ((I=$[I0+1]; I<$[I0+NC]; I++ ))
  do
    fout=${fname##*/}
    fout=$WORKPATH/${fout%\.*\.bev}.$I
    if [[ $BG ]]; then
      echo "bev2xy $fname $[I-I0] | lvq2agf -c -h $fout &"
      bev2xy $fname $[I-I0] | lvq2agf -c -h $fout &
    else
      echo "bev2xy $fname $[I-I0] | lvq2agf -c -h $fout"
      bev2xy $fname $[I-I0] | lvq2agf -c -h $fout
    fi
  done 
done
