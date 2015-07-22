#!/bin/bash
set -e

#defaults:
RANGE="0/360/0/90"
RANGE2="0/360/-90/90"
PROJ="S0/90/18"
ORIENT=-ZBLx

NLON=144
NLAT=73

Z0=-1
ZF=1
NZ=21

while getopts 'qc:g:z:I:J:F:Har:x:y:Z:R:2:' ARG; do
  case $ARG in
    Z) ORIENT=-Z$OPTARG
      ;;
    z) ZOPTS="$ZOPTS -z $OPTARG"
       NZ=$OPTARG
       WCFLAG=1
      ;;
    I) ZOPTS="$ZOPTS -I $OPTARG"
       Z0=$OPTARG
       WCFLAG=1
      ;;
    F) ZOPTS="$ZOPTS -F $OPTARG"
       ZF=$OPTARG
       WCFLAG=1
      ;;
    c) PALETTE=$OPTARG
       DDSWTC=1;		
       #dont delete the contour palette before exiting
      ;;
    g) ZOPTS="$ZOPTS -g"
       WCFLAG=1
      ;;
    a) AOPT="-a";
      ;;
    r) EXT_OPT="$EXT_OPT -r $OPTARG"
      ;;
    x) NLON=$OPTARG
      ;;
    y) NLAT=$OPTARG
      ;;
    e) CFLAG=1
      ;;
    q) QFLAG=1
      ;;
    R) RANGE=$OPTARG
      ;;
    2) RANGE2=$OPTARG
      ;;
    J) PROJ=$OPTARG
      ;;
    H) echo "Usage: quick_plot [-a] [-q] [-r slen] "
       echo "            [-z nz] [-I bottom] [-F top] [-g]"
       echo "            [-x nlon] [-y nlat] [outfile]"
       echo "options:"
       echo "  -J projection"
       echo "  -R range"
       echo "  -H help"
       exit 0
       ;;
  esac
done

ZOPTS="-z $NZ -I $Z0 -F $ZF"

shift $(($OPTIND - 1))

INFILE=/dev/stdin
#INFILE=dum.in
#INDEX=$2
#OUTFILE=dum.ps

NOW=$(date +"%Y.%m.%d_%H-%M-%S")
BASE=dum.$NOW
#BASE=${INFILE%.*}
#BASE=${BASE##*/}
OUTFILE=${1:-$BASE.ps}

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

echo "psbasemap -R${RANGE} -J${PROJ} -K -Bg30 > ${OUTFILE}"
psbasemap -R${RANGE} -J${PROJ} -K -Bg30 > ${OUTFILE}

  DLON=$(echo "scale=2; 360./$NLON" | bc)
  DLAT=$(echo "scale=2; 180./($NLAT-1)" | bc)
  echo $NLON $NLAT $DLON $DLAT

  GRDFILE=$BASE.grd

  echo "gen_zgrid $ZOPTS > $ZGRIDFILE"
  gen_zgrid $ZOPTS > $ZGRIDFILE
  echo "makecpt -T$ZGRIDFILE > $PALETTE"
  makecpt -T$ZGRIDFILE > $PALETTE

  echo "xyz2grd -R$RANGE2 -I$DLON/$DLAT $ORIENT -G$GRDFILE < $INFILE"
  xyz2grd -R$RANGE2 -I$DLON/$DLAT $ORIENT -G$GRDFILE < $INFILE
  echo "grdimage $GRDFILE -R$RANGE -J$PROJ -C$PALETTE -O -K >> $OUTFILE"
  grdimage $GRDFILE -R$RANGE -J$PROJ -C$PALETTE -O -K >> $OUTFILE

echo "pscoast -R${RANGE} -J${PROJ} -Dl -W -O >> ${OUTFILE}"
pscoast -R${RANGE} -J${PROJ} -Dl -W -O >> ${OUTFILE}

rm -f $ZGRIDFILE

rm -f $GRDFILE

if test -z $1; then 
	gv ${OUTFILE}
	rm ${OUTFILE};
fi

if test -z $DDSWTC; then 
  rm $PALETTE;
fi

