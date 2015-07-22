#!/bin/bash

#set -x

#WORKPATH=/work/patmills/ctraj_work
WORKPATH=animation_frames
SKIP=1

while getopts 'aCgHq0:c:F:I:J:N:O:p:r:R:x:y:z:' DUM
do
  case $DUM in
    0) I0=$OPTARG
      ;;
    a) OPTS="$OPTS -a"
      ;;
    c) CFILE=$OPTARG
      ;;
    C) CLEAN_FLAG=1
      ;;
    F) ZOPTS="$ZOPTS -F $OPTARG"
      ;;
    g) OPTS="$ZOPTS -g"
    	QFLAG=1;
      ;;
    I) ZOPTS="$ZOPTS -I $OPTARG"
      ;;
    J) ZOPTS="$OPTS -J $OPTARG"
      ;;
    N) N=$OPTARG
      ;;
    O) SKIP=$OPTARG
      ;;
    p) WORKPATH=$OPTARG
      ;;
    q) ZOPTS="$ZOPTS -q"
    	QFLAG=1;
      ;;
    r) OPTS="$OPTS -r $OPTARG"
      ;;
    R) ZOPTS="$OPTS -R $OPTARG"
      ;;
    x) OPTS="$OPTS -x $OPTARG"
      ;;
    y) OPTS="$OPTS -y $OPTARG"
      ;;
    z) ZOPTS="$ZOPTS -z $OPTARG"
      ;;
    H) echo "Usage: animate_ctraj [-O I0] [-N N] [-a] [-q] [-r slen]" 
       echo "            [-x nlon] [-y nlat] datafile outfile"
       echo "options:"
       echo "  -0 initial index"
       echo "  -c colour palette"
       echo "  -C clean up"
       echo "  -F top-most contour"
       echo "  -g use geometric progression"
       echo "  -I bottom-most contour"
       echo "  -J projection"
       echo "  -N number of frames"
       echo "  -O stride"
       echo "  -p work-path [$WORKPATH]"
       echo "  -R range"
       echo "  -H help"
       exit 0
      ;;
  esac
done

shift $(($OPTIND - 1))
INFILE=$1
OUTFILE=$2

BASE=${INFILE%.*}
BASE=${BASE##*/}

if test -z $I0; then I0=0; fi

if test $N; then IF=$[I0+N-1]; fi

if test -z $QFLAG; then 
	MAX=$(bev2xy $INFILE | wc -l);
fi

if test -z $MAX; then
	QFLAG=1;	
else
  if [ $MAX -eq 0 ]; then
	QFLAG=1;	
  fi;
fi

if test -z $CFILE; then
  CFILE=$BASE.cpt;
fi


if [ $QFLAG ]; then MAX=$(extract_field $OPTS $INFILE); fi

echo $I0

echo $INFILE
echo $MAX
echo $IF

if test -z $IF; then IF=$((MAX-1)); fi

if [ $IF -ge $MAX ]
then 
    IF=$((MAX-1))
fi

echo $IF
echo $((MAX-1))

RANGE=-180/180/0/90
PROJ=S0/90/18

SHIFT=-Yr2

if [ ! -d $WORKPATH ]; then
  mkdir $WORKPATH
  if [ $CLEAN_FLAG ]; then WORK_CLEAN=1; fi
fi

FRAME=$WORKPATH/$BASE.$I0
plot_frame.sh -c $CFILE $ZOPTS $OPTS $INFILE $I0 $FRAME.ps
LIST="$LIST $FRAME.gif"
convert -background white $FRAME.ps $FRAME.gif

for ((INDEX=$[I0+1]; INDEX<=$IF; INDEX+=SKIP)); do
  FRAME=$WORKPATH/$BASE.$INDEX
  plot_frame.sh -c $CFILE $ZOPTS $OPTS $INFILE $INDEX $FRAME.ps
  LIST="$LIST $FRAME.gif"
  convert -background white $FRAME.ps $FRAME.gif
done

convert -background white -loop 0 -dispose previous -delay 10 $LIST $OUTFILE

if [ $CLEAN_FLAG ]; then
  for ((INDEX=$I0; INDEX<=$IF; INDEX+=SKIP)); do
    FRAME1=$WORKPATH/$BASE.$INDEX
    rm $FRAME1.ps $FRAME1.gif
  done
  if [ $WORK_CLEAN ]; then rmdir $WORKPATH; fi
fi

