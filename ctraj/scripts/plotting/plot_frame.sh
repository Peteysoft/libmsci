#!/bin/bash
#set -x
set -e

RANGE="0/360/0/90"
RANGE2="0/360/-90/90"
PROJ="S0/90/18"

NLON=360
NLAT=181

XLEN=360
YLEN=180

VTYPE=2

while getopts 'qc:g:z:I:J:F:Har:x:y:V:n:qX:Y:R:' ARG; do
  case $ARG in
    z) ZOPTS="$ZOPTS -z $OPTARG"
       WCFLAG=1
      ;;
    I) ZOPTS="$ZOPTS -I $OPTARG"
       WCFLAG=1
      ;;
    F) ZOPTS="$ZOPTS -F $OPTARG"
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
    V) EXT_OPT="$EXT_OPT -V $OPTARG"
       VTYPE=$OPTARG
      ;;
    n) EXT_OPT="$EXT_OPT -n $OPTARG"
      ;;
    X) EXT_OPT="$EXT_OPT -X $OPTARG"
       XLEN=$OPTARG
      ;;
    Y) EXT_OPT="$EXT_OPT -Y $OPTARG"
       YLEN=$OPTARG
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
    J) PROJ=$OPTARG
      ;;
    H) echo "Usage: plot_frame [-a] [-q] [-r slen] "
       echo "            [-z nz] [-I bottom] [-F top] [-g]"
       echo "            [-x nlon] [-y nlat] datafile [index [outfile]]"
       echo "options:"
       echo "  -J projection"
       echo "  -R range"
       echo "  -H help"
       exit 0
       ;;
  esac
done

shift $(($OPTIND - 1))

INFILE=$1
INDEX=$2
OUTFILE=$3

NOW=$(date +"%Y.%m.%d_%H-%M-%S")
BASE=${OUTFILE%.*}.$NOW
#BASE=${INFILE%.*}
#BASE=${BASE##*/}

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

if test -z $INDEX; then INDEX=$((0)); fi

if test -z $CFLAG; then
  N=0
  if test -z $QFLAG; then
    N=$(bev2xy $INFILE | wc -l);
  fi

  if test $[N] -eq 0; then
    QFLAG=1; 
  fi

  if test $QFLAG; then
    N=$(extract_field $EXT_OPT $INFILE);
  fi

  echo "plot frame: n=$N"

  if test $[INDEX] -ge $[N]; then
	echo "Index out of bounds: $INDEX; $N records";
  fi
fi

if test $[INDEX] -lt 0; then
	echo "Index must be greater than 0: $INDEX; $N records";
fi

if test -z $OUTFILE; then OUTFILE="$BASE.$INDEX.tmp.ps"; fi

echo $QFLAG
echo "psbasemap -R${RANGE} -J${PROJ} -K -Bg30 > ${OUTFILE}"
psbasemap -R${RANGE} -J${PROJ} -K -Bg30 > ${OUTFILE}

echo $QFLAG

if test $QFLAG; then
  DLON=$(date_calc "($XLEN|$NLON)")
  DLAT=$(date_calc "($YLEN|($NLAT_1))")

  GRDFILE=$BASE.$INDEX.grd

  ORIENT=-ZBLx

  if test $WCFLAG; then
    echo "extract_field -x $NLON -y $NLAT $EXT_OPT $INFILE $INDEX | gen_zgrid $ZOPTS > $ZGRIDFILE"
    extract_field -x $NLON -y $NLAT $EXT_OPT $INFILE $INDEX | gen_zgrid $ZOPTS > $ZGRIDFILE
    echo "makecpt -T$ZGRIDFILE > $PALETTE"
    makecpt -T$ZGRIDFILE > $PALETTE
  fi

  NLON2=$(( $NLON + 1 ))

  echo "extract_field -x $NLON -y $NLAT $EXT_OPT $INFILE $INDEX | xyz2grd -R$RANGE2 -I$NLON+/$NLAT+ $ORIENT -G$GRDFILE"
  extract_field -x $NLON -y $NLAT $EXT_OPT $INFILE $INDEX | xyz2grd -R$RANGE2 -I$DLON/$DLAT $ORIENT -G$GRDFILE
  echo "grdimage $GRDFILE -R$RANGE -J$PROJ -C$PALETTE -O -K >> $OUTFILE;"
  grdimage $GRDFILE -R$RANGE -J$PROJ -C$PALETTE -O -K >> $OUTFILE;
else
  echo "bev2xy ${INFILE} ${INDEX} | psxy -R${RANGE} -J${PROJ} -W5,black -O -K >> ${OUTFILE};"
  bev2xy ${INFILE} ${INDEX} | psxy -R${RANGE} -J${PROJ} -W5,black -O -K >> ${OUTFILE};
fi

if [ $VTYPE == 2 ]
then
  echo "pscoast -R${RANGE} -J${PROJ} -Dl -W -O >> ${OUTFILE}"
  pscoast -R${RANGE} -J${PROJ} -Dl -W -O >> ${OUTFILE}
fi

rm -f $ZGRIDFILE

rm -f $GRDFILE

if test -z $3; then 
	gv ${OUTFILE}
	rm ${OUTFILE};
fi

if test -z $DDSWTC; then 
  if test $QFLAG; then
	rm $PALETTE;
  fi;
fi

