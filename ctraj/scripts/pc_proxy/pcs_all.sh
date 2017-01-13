#!/bin/bash

#gets all the pcs for a sparse array...

SKIP=1
N=-1
while getopts 'O:N:v:A:lS:' ARG; do
  case $ARG in
    O) OPTS="$OPTS -O $OPTARG"
       WCFLAG=1
      ;;
    N) N=$OPTARG
      ;;
    v) OPTS="$OPTS -v $OPTARG"
       WCFLAG=1
      ;;
    A) OPTS="$OPTS -A $OPTARG"
       WCFLAG=1
      ;;
    l) OPTS="$OPTS -l";
      ;;
    S) SKIP=$OPTARG
      ;;
    h) echo "Usage: pcs_all.sh [-O start] [-N N] [-v nv] [-A nA] [-l]"
       echo "            filename outbase"
       exit -1
       ;;
  esac
done

shift $((OPTIND-1))

MATFILE=$1
BASE=$2

NALL=$(sparse_mat_prod -N $MATFILE)
if [ $N -lt 0 ]; then
  N=$NALL
fi
if [ $N -gt $NALL ]; then
  N=$NALL
fi

I0=$SKIP

for ((INDEX=$I0; INDEX<=$N; INDEX+=$SKIP)); do
  OUTFILE=$BASE$INDEX.dat
  #echo "svd_sparse_array $OPTS -N $INDEX $MATFILE $OUTFILE"
  svd_sparse_array $OPTS -N $INDEX $MATFILE $OUTFILE

  #doesn't work smoothly since I made the forward slash a special character...
  #OUTFILE=$INDEX.dat
  #echo "$OUTFILE=svd($MATFILE($I0:$INDEX)*transpose($MATFILE($I0:$INDEX)))" | sparse_calc -p $BASE
done


