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

#go through gmt "master" command:
PREFIX="gmt "

#command for contour plotting:
CONTOUR=${PREFIX}grdimage

#plot continents using lines?
LFLAG=1
CUSTOMR=0

while getopts 'c:F:g:h:I:J:R:T:V:W:x:X:y:Y:z:HLqS' ARG; do
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
    h) if [[ $OPTARG -lt 0 ]]
       then
         RANGE="0/360/-90/0"
         PROJ="S0/-90/18"
       elif [[ $OPTARG -eq 0 ]]
       then
         RANGE="0/360/-90/90"
         PROJ="Q0/0/9i"
       fi
      ;;
    I) Z0=$OPTARG
       WCFLAG=1
      ;;
    J) PROJ=$OPTARG
      ;;
    L) CONTOUR=${PREFIX}grdcontour
       LFLAG=0
      ;;
    q) QFLAG=1
      ;;
    R) RANGE=$OPTARG
       CUSTOMR=1
      ;;
    S) SFLAG=1
      ;;
    T) TITLE=:.\"$OPTARG\":
      ;;
    V) VTYPE=$OPTARG
      ;;
    W) THICK=$OPTARG
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
    H) echo "Usage: quick_plot [-q] [-L] [-V 1] [-h hemi]"
       echo "            [-z nz] [-I bottom] [-F top] [-g]"
       echo "            [-x nlon] [-y nlat] [outfile]"
       echo "options:"
       echo "  -J projection"
       echo "  -R range"
       echo "  -L plot lines instead of solid contours"
       echo "  -W width of line/size of symbol"
       echo "  -S use symbols instead of lines"
       echo "  -q plot field instead of lines/sybols"
       echo "  -I bottom contour level"
       echo "  -F top contour level"
       echo "  -z number of contour levels"
       echo "  -g levels have logarithmic progression"
       echo "  -h hemisphere: 1=N (default); 0=globe; -1=S"
       echo "  -T title"
       echo "  -H help"
       echo "  -c colour palette file"
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
echo "${PREFIX}psbasemap -R${RANGE} -J${PROJ} -Bg30$TITLE -K > ${OUTFILE}"
${PREFIX}psbasemap -R${RANGE} -J${PROJ} -Bg30 -K > ${OUTFILE}

echo $QFLAG

#add coastlines:
if test -z $QFLAG; then
  LFLAG=0
fi

if [[ $VTYPE -eq 0 && $LFLAG -eq 0 ]]
then
  echo "---${PREFIX}pscoast -R${RANGE} -J${PROJ} -G220 -Dl -K -O >> ${OUTFILE}---"
  #${PREFIX}pscoast -R${RANGE} -J${PROJ} -G220 -Dl -K -O >> ${OUTFILE}
fi

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
    ${PREFIX}makecpt -T$ZGRIDFILE > $PALETTE
  fi

  echo "xyz2grd -R$RANGE2 -I$DLON/$DLAT $ORIENT -G$GRDFILE"
  ${PREFIX}xyz2grd -R$RANGE2 -I$DLON/$DLAT $ORIENT -G$GRDFILE
  echo "$CONTOUR $GRDFILE -R$RANGE -J$PROJ -C$PALETTE -O -K >> $OUTFILE;"
  $CONTOUR $GRDFILE -R$RANGE -J$PROJ -C$PALETTE -O -K >> $OUTFILE;
else
  #plot single contour:
  if test $SFLAG; then
    echo "${PREFIX}psxy -R${RANGE} -J${PROJ} -Sc${THICK}p -O -K >> ${OUTFILE};"
    ${PREFIX}psxy -R${RANGE} -J${PROJ} -Sc${THICK}p -O -K >> ${OUTFILE};
  else
    echo "${PREFIX}psxy -R${RANGE} -J${PROJ} -W$THICK,black -O -K >> ${OUTFILE};"
    ${PREFIX}psxy -R${RANGE} -J${PROJ} -W$THICK,black -O -K >> ${OUTFILE};
  fi
fi


#add coastlines:
if [[ $VTYPE -eq 0 && $LFLAG -eq 1 ]]
then
  echo "${PREFIX}psbasemap -R${RANGE} -J${PROJ} -O -Bg30 >> ${OUTFILE}"
  ${PREFIX}psbasemap -R${RANGE} -J${PROJ} -O -Bg30 >> ${OUTFILE}
  #${PREFIX}psbasemap -R${RANGE} -J${PROJ} -K -O -Bg30 >> ${OUTFILE}
  echo "----${PREFIX}pscoast -R${RANGE} -J${PROJ} -Dl -W -O >> ${OUTFILE}----"
  #${PREFIX}pscoast -R${RANGE} -J${PROJ} -Dl -W -O >> ${OUTFILE}
else
  # reinforce grid lines:
  echo "${PREFIX}psbasemap -R${RANGE} -J${PROJ} -O -Bg30 >> ${OUTFILE}"
  ${PREFIX}psbasemap -R${RANGE} -J${PROJ} -O -Bg30 >> ${OUTFILE}
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

