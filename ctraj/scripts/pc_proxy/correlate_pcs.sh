#!/bin/bash

#gets all the pcs for a sparse array...

#set -x

NV=5
#NGRID=200
SKIP=1
N=-1

while getopts 'Hn:v:N:' ARG; do
  case $ARG in
    v) NV=$OPTARG
      ;;
    n) NGRID=$OPTARG
      ;;
    N) N=$OPTARG
      ;;
    S) SKIP=$OPTARG
      ;;
    H) echo "Usage: correlate_pcs.sh [-n ngrid] [-v nv] basename tracer"
       exit -1
       ;;
  esac
done

shift $(($OPTIND - 1))

BASE=$1
TRACER_FILE=$2
NALL=$(extract_field $TRACER_FILE)
if [ $N < 0 || $N > $NALL ]; then
  N=$NALL
fi

echo "$N records found in $TRACER_FILE"
I0=$SKIP

INDEX1=$[NV-1]

#REC=$(extract_field -Rn $NGRID)
for ((INDEX=$I0; INDEX<$N; INDEX+=$SKIP)); do
  OUTFILE=$BASE$INDEX.dat
  correlate_fields -1 $INDEX1 -2 $INDEX $OUTFILE $TRACER_FILE
done

