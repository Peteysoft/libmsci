#!/bin/bash

set -xe

NLON=360
NLAT=181

STRIDE=1

ZGRIDFILE=confidence.cpt

WORKPATH=animation_frames
#WORKPATH=/work/patmills/ctraj_work/animation_frames

PREFIX="gmt "

while getopts '0:I:C:N:O:p:r:F:w:x:y:z:acgH' DUM
do
  case $DUM in
    I) ZOPTS="$ZOPTS -I $OPTARG"
      ;;
    g) ZOPTS="$ZOPTS -g"
      ;;
    p) WORKPATH=$OPTARG
      ;;
    F) ZOPTS="$ZOPTS -F $OPTARG"
      ;;
    z) ZOPTS="$ZOPTS -z $OPTARG"
      ;;
    x) NLON=$OPTARG
      ;;
    y) NLAT=$OPTARG
      ;;
    0) I0=$OPTARG
      ;;
    N) N=$OPTARG
      ;;
    O) STRIDE=$OPTARG
      ;;
    c) PALETTE=$OPTARG
       DDSWTC=1;		
       #dont delete the contour palette before exiting
      ;;
    C) CLEAN_FLAG=1		#clean up afterwards
      ;;
    # these three options are obsolete/unnecessary:
    a) EOPTS="$EOPTS-a";
      ;;
    r) EOPTS="$EOPTS -r $OPTARG"
      ;;
    n) EOPTS="$EOPTS -n $OPTARG"
      ;;
    H) echo "Usage: animate_combined [-0 I0] [-N N] [-a] [-r slen]" 
       echo "            [-x nlon] [-y nlat] tracer contour outfile"
       echo "options:"
       echo "  -0 initial index"
       echo "  -I bottom-most contour"
       echo "  -c colour palette"
       echo "  -C clean up"
       echo "  -g use geometric progression"
       echo "  -J projection"
       echo "  -N number of frames"
       echo "  -O stride"
       echo "  -p work-path [$WORKPATH]"
       echo "  -R range"
       echo "  -F top-most contour"
       echo "  -H help"
       exit 0
      ;;
    w) WORKPATH=$OPTARG
      ;;
  esac
done

shift $(($OPTIND - 1))
TRACER=$1
CONTOUR=$2
OUTFILE=$3

if [ ! -d $WORKPATH ]; then
  mkdir $WORKPATH
  WORK_CLEAN=1
fi

BASE=${OUTFILE%.*}
BASE=$WORKPATH/${BASE##*/}

if test -z $PALETTE; then
  PALETTE=$BASE.cpt;
fi

ZGRIDFILE=$BASE.zgrid;

if test -z $I0; then I0=0; fi

N1=$(bev2xy $CONTOUR | wc -l);

N2=$(extract_field $ROPT $TRACER)

echo $I0

echo $N1 $N2

if test -z $N; then 
  IF=$(($N1-1)); 
else
  IF=$(($N+I0-1));
fi

if test $IF -ge $N1; then IF=$(($N1-1)); fi
if test $IF -ge $N2; then IF=$(($N2-1)); fi

RANGE=0/360/0/90
PROJ=S0/90/18

SHIFT=-Yr2

RANGE2="0/360/-90/90"

#generate colour palette:
extract_field -x $NLON -y $NLAT $EXT_OPT $TRACER 0 | gen_zgrid $ZOPTS > $ZGRIDFILE
${PREFIX}makecpt -T$ZGRIDFILE > $PALETTE

for ((INDEX=$I0; INDEX<=$IF; INDEX+=$STRIDE)); do
  FRAME1=$BASE.$INDEX
  FRAME=$FRAME1.ps

  ${PREFIX}psbasemap -R${RANGE} -J${PROJ} -K -Bg30 > ${FRAME}

  DLON=$((360/$NLON))
  DLAT=$((180/($NLAT-1)))

  GRDFILE=$BASE.$INDEX.grd

  ORIENT=-ZBLx

  extract_field -x $NLON -y $NLAT $EOPTS $TRACER $INDEX | ${PREFIX}xyz2grd -R$RANGE2 -I$DLON/$DLAT $ORIENT -G$GRDFILE
  ${PREFIX}grdimage $GRDFILE -R$RANGE -J$PROJ -C$PALETTE -O -K >> $FRAME;

  bev2xy $CONTOUR $INDEX | ${PREFIX}psxy -R$RANGE -J$PROJ -W2,black -O -K >> $FRAME;

  #${PREFIX}pscoast -R$RANGE -J$PROJ -Dl -W -O >> $FRAME

  LIST="$LIST $FRAME1.gif"
  convert -background white $FRAME $FRAME1.gif

  rm $GRDFILE
done

rm $ZGRIDFILE
if test -z $DDSWTC; then
  rm $PALETTE;
fi

convert -loop 0 -dispose previous -delay 10 $LIST $OUTFILE

if [ $CLEAN_FLAG ]; then
  for ((INDEX=$I0; INDEX<=$IF; INDEX+=$STRIDE)); do
    FRAME1=$WORKPATH/$BASE.$INDEX
    rm $FRAME1.ps $FRAME1.gif
  done
  if [ $CLEAN_WORK ] 
  then
    rm $WORKPATH
  fi
fi
  


