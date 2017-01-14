#!/bin/bash
#set -x
set -e

#default values for GMT projection:
RANGE="0/360/0/90"
RANGE2="0/360/-90/90"
PROJ="S0/90/18"

#number of longitude and latitude grids:
NLON=360
NLAT=181

#dimension of frame:
XLEN=360
YLEN=180

#add outlines of landforms:
VTYPE=0

#vertical grid:
Z0=-1
ZF=1
NZ=21

#line thickness for line plots:
THICK=2

#command for contour plotting:
CONTOUR=grdimage

#plot continents using lines?
LFLAG=1
CUSTOMR=0

while getopts 'c:F:g:I:J:R:T:V:x:X:y:Y:z:HLq' ARG; do
  case $ARG in
    c) PALETTE=$OPTARG
       DDSWTC=1;		
       #dont delete the contour palette before exiting
      ;;
    F) ZF=$OPTARG
       WCFLAG=1			#auto-generate colour palette
      ;;
    g) ZOPTS="$ZOPTS -g"
       WCFLAG=1
      ;;
    I) Z0=$OPTARG
       WCFLAG=1
      ;;
    J) PROJ=$OPTARG
      ;;
    L) CONTOUR=grdcontour
       LFLAG=0
      ;;
    q) QFLAG=1
      ;;
    R) RANGE=$OPTARG
       CUSTOMR=1
      ;;
    T) THICK=$OPTARG
      ;;
    V) VTYPE=$OPTARG
      ;;
    x) NLON=$OPTARG
      ;;
    X) XLEN=$OPTARG
       CUSTOMRANGE=1
      ;;
    y) NLAT=$OPTARG
      ;;
    Y) YLEN=$OPTARG
       CUSTOMRANGE=1
      ;;
    z) NZ=$OPTARG
       WCFLAG=1
      ;;
    H) echo "Usage: quick_plot [-q] [-L] [-V 1]"
       echo "            [-z nz] [-I bottom] [-F top] [-g]"
       echo "            [-x nlon] [-y nlat] [outfile]"
       echo "options:"
       echo "  -J projection"
       echo "  -R range"
       echo "  -L plot lines instead of solid contours"
       echo "  -H help"
       exit 0
       ;;
  esac
done

shift $(($OPTIND - 1))

OUTFILE=$1

ZOPTS="$ZOPTS -z $NZ -I $Z0 -F $ZF"

# for non-global and non-projection plots:
# (this will need to be sorted...)
if test $CUSTOMRANGE; then
  RANGE2=0/$XLEN/0/$YLEN
  if [[ $VTYPE -ne 0 && $CUSTOMR -eq 0 ]]; then
    RANGE=$RANGE2
  fi
fi

NOW=$(date +"%Y.%m.%d_%H-%M-%S")
BASE=${OUTFILE%.*}.$NOW
#BASE=${INFILE%.*}
#BASE=${BASE##*/}

#colour palette stuff:
if test -z $DDSWTC; then
  WCFLAG=1;
fi

if test -z $PALETTE; then
  PALETTE=$BASE.cpt;
fi

if [ ! -f $PALETTE ]; then
  WCFLAG=1;
fi

ZGRIDFILE=$BASE.zgrid;

if test -z $OUTFILE; then OUTFILE="$BASE.$INDEX.tmp.ps"; fi

echo $QFLAG
echo "psbasemap -R${RANGE} -J${PROJ} -Bg30 -K > ${OUTFILE}"
psbasemap -R${RANGE} -J${PROJ} -Bg30 -K > ${OUTFILE}

echo $QFLAG

#add coastlines:
if test -z $QFLAG; then
  LFLAG=0
fi

if [[ $VTYPE -eq 0 && $LFLAG -eq 0 ]]
then
  echo "pscoast -R${RANGE} -J${PROJ} -Dl -K -O >> ${OUTFILE}"
  pscoast -R${RANGE} -J${PROJ} -G220 -Dl -K -O >> ${OUTFILE}
fi

echo "psbasemap -R${RANGE} -J${PROJ} -K -O -Bg30 >> ${OUTFILE}"
psbasemap -R${RANGE} -J${PROJ} -K -O -Bg30 >> ${OUTFILE}

if test $QFLAG; then
  #plot field:
  DLON=$(date_calc "($XLEN|$NLON)")
  DLAT=$(date_calc "($YLEN|(${NLAT}_1))")

  GRDFILE=$BASE.$INDEX.grd

  ORIENT=-ZBLx

  #auto-generate colour palette/vertical levels:
  if test $WCFLAG; then
    echo "gen_zgrid $ZOPTS > $ZGRIDFILE"
    gen_zgrid $ZOPTS > $ZGRIDFILE
    echo "makecpt -T$ZGRIDFILE > $PALETTE"
    makecpt -T$ZGRIDFILE > $PALETTE
  fi

  echo "xyz2grd -R$RANGE2 -I$DLON/$DLAT $ORIENT -G$GRDFILE"
  xyz2grd -R$RANGE2 -I$DLON/$DLAT $ORIENT -G$GRDFILE
  echo "$CONTOUR $GRDFILE -R$RANGE -J$PROJ -C$PALETTE -O -K >> $OUTFILE;"
  $CONTOUR $GRDFILE -R$RANGE -J$PROJ -C$PALETTE -O -K >> $OUTFILE;
else
  #plot single contour:
  echo "psxy -R${RANGE} -J${PROJ} -W$THICK,black -O -K >> ${OUTFILE};"
  psxy -R${RANGE} -J${PROJ} -W$THICK,black -O -K >> ${OUTFILE};
fi

#add coastlines:
if [[ $VTYPE -eq 0 && $LFLAG -eq 1 ]]
then
  echo "pscoast -R${RANGE} -J${PROJ} -Dl -W -O >> ${OUTFILE}"
  pscoast -R${RANGE} -J${PROJ} -Dl -W -O >> ${OUTFILE}
fi

#clean up:
rm -f $ZGRIDFILE

rm -f $GRDFILE

if test -z $1; then 
	gv ${OUTFILE}
	rm ${OUTFILE};
fi

if test -z $DDSWTC; then 
  if test $QFLAG; then
	rm $PALETTE;
  fi;
fi

