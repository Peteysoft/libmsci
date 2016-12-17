#!/bin/bash

while getopts "a:h:i:I:k:l:N:q:Q:r:s:t:v:V:W:gHKnu" ARG
do
  case $ARG in
    g) F2COPTS="$F2COPTS -g"
       ;;
    q) F2COPTS="$F2COPTS -q $OPTARG"
       N=$OPTARG
       ;;
    Q) QFLAG=$OPTARG
       ;;
    k) AGFOPTS="$AGFOPTS -k $OPTARG"
       ;;
    r) AGFOPTS="$AGFOPTS -r $OPTARG"
       ;;
    s) AGFOPTS="$AGFOPTS -s $OPTARG"
       ;;
    t) AGFOPTS="$AGFOPTS -t $OPTARG"
       ;;
    v) AGFOPTS="$AGFOPTS -v $OPTARG"
       ;;
    V) AGFOPTS="$AGFOPTS -V $OPTARG"
       ;;
    W) AGFOPTS="$AGFOPTS -W $OPTARG"
       ;;
    a) NORMOPTS="$NORMOPTS -a $OPTARG"
       ;;
    K) NORMOPTS="$NORMOPTS -K"
       ;;
    n) NORMOPTS="$NORMOPTS -n"
       ;;
    S) NORMOPTS="$NORMOPTS -S $OPTARG"
       ;;
    u) NORMOPTS="$NORMOPTS -u"
       ;;
    h) AGFOPTS="$AGFOPTS -h $OPTARG"
       ;;
    i) AGFOPTS="$AGFOPTS -i $OPTARG"
       ;;
    I) AGFOPTS="$AGFOPTS -I $OPTARG"
       ;;
    l) AGFOPTS="$AGFOPTS -l $OPTARG"
       ;;
    N) AGFOPTS="$AGFOPTS -N $OPTARG"
       ;;
    H) echo "Training for continuum retrievals using multi-borders"
       echo
       echo "syntax: c_borders.sh [-g] [-q n] [-Q i] [trainopt] [...] \\"
       echo "                 control train border out \\"
       echo "                 {[min max] | [thresh1 [thresh2 [thresh3 ... ]]]}"
       echo
       echo "where:"
       echo "  control        control file"
       echo "  train          base filename for binary training data:"
       echo "                   .vec for vector coordinates"
       echo "                   .dat for floating point ordinates"
       echo "                   .std for transformation/normalization matrix"
       echo "                        (unless -a specified)"
       echo "  border         base name for files containing border data"
       echo "                   (file names are recorded in output control file)"
       echo "  out            output control file to feed to classify_c"
       echo "  min            lowest division in series"
       echo "  max            highest division in series"
       echo "  threshN        Nth threshold in discretization series"
       echo
       echo "  * commands to train the model are written to stdout"
       echo
       echo "options:"
       echo "  -H       = print this help screen"
       echo "  -q     n = is the number divisions"
       echo "  -g       = use a geometric progression"
       echo "  -Q     1 = rewrite control file with strictly hierarchical setup"
       echo "  -Q     2 = rewrite control file with strictly non-hierarchical setup "
       echo
       echo "trainopt:"
       class_borders | tail -n 49 | head -n 27
       exit
       ;;
  esac
done

shift $(($OPTIND - 1))

CONTROL=$1
TRAIN=$2
BORDER=$3
OUT=$4

if [ -z $N ]
then
  N=$(($#-3))
fi

case $QFLAG in
  1) rm -f $CONTROL
     print_control -Q 0 $N
     ;;
  2) rm -f $CONTROL
     print_control -Q 2 $N
     ;;
esac

#echo "multi_borders $NORMOPTS $AGFOPTS $CONTROL $TRAIN $BORDER $OUT"
multi_borders $NORMOPTS $AGFOPTS $CONTROL $TRAIN $BORDER $OUT

for ((I=5; I<=$#; I++)); do
  F2COPTS2="$F2COPTS2 ${!I}"
done

#echo "float_to_class $F2COPTS $TRAIN.dat $TRAIN.cls $F2COPTS2 >> $OUT"
float_to_class $F2COPTS $TRAIN.dat $TRAIN.cls $F2COPTS2 >> $OUT

