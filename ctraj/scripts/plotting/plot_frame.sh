#!/bin/bash
#set -x
set -e

NLON=360
NLAT=181

XLEN=360
YLEN=180

VTYPE=0

CONTOUR=grdimage

while getopts 'c:F:g:h:I:J:n:r:R:T:V:W:x:X:y:Y:z:HLq' ARG; do
  case $ARG in
    c) PALETTE=$OPTARG
       DDSWTC=1;		
       #dont delete the contour palette before exiting
      ;;
    F) ZOPTS="$ZOPTS -F $OPTARG"
       WCFLAG=1			#auto-generate colour palette
      ;;
    g) ZOPTS="$ZOPTS -g"
       WCFLAG=1
      ;;
    h) PFOPTS="$PFOPTS -h $OPTARG"
      ;;
    I) ZOPTS="$ZOPTS -I $OPTARG"
       WCFLAG=1
      ;;
    # options to pass to "quick_plot":
    J) PFOPTS="$PFOPTS -J $OPTARG"
      ;;
    L) CONTOUR=grdcontour
       PFOPTS="$PFOPTS -L"
      ;;
    n) EXT_OPT="$EXT_OPT -n $OPTARG"
      ;;
    q) QFLAG=1
      ;;
    r) EXT_OPT="$EXT_OPT -r $OPTARG"
      ;;
    R) PFOPTS="$PFOPTS -R $OPTARG"
      ;;
    T) PFOPTS="$PFOPTS -T \"$OPTARG\""
      ;;
    V) EXT_OPT="$EXT_OPT -V $OPTARG"
       VTYPE=$OPTARG
      ;;
    W) PFOPTS="$PFOPTS -W $OPTARG"
      ;;
    x) NLON=$OPTARG
       PFOPTS="$PFOPTS -x $OPTARG"
      ;;
    X) EXT_OPT="$EXT_OPT -X $OPTARG"
       XLEN=$OPTARG
       PFOPTS="$PFOPTS -X $OPTARG"
      ;;
    y) NLAT=$OPTARG
       PFOPTS="$PFOPTS -y $OPTARG"
      ;;
    Y) EXT_OPT="$EXT_OPT -Y $OPTARG"
       YLEN=$OPTARG
       PFOPTS="$PFOPTS -X $OPTARG"
      ;;
    z) ZOPTS="$ZOPTS -z $OPTARG"
       WCFLAG=1
      ;;
    H) echo "Usage: plot_frame [-C] [-q] [-r slen] [-h hemi] "
       echo "            [-z nz] [-I bottom] [-F top] [-g]"
       echo "            [-x nlon] [-y nlat] datafile [index [outfile]]"
       echo "options:"
       echo "  -J projection"
       echo "  -R range"
       echo "  -C plot lines instead of solid contours"
       echo "  -H help"
       exit 0
       ;;
  esac
done

PFOPTS="$PFOPTS -x $NLON -y $NLAT -V $VTYPE"

shift $(($OPTIND - 1))

INFILE=$1
INDEX=$2
OUTFILE=$3

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

#indexing into data file:
if test -z $INDEX; then INDEX=$((0)); fi

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

if test $[INDEX] -lt 0; then
	echo "Index must be greater than 0: $INDEX; $N records";
fi

if test -z $OUTFILE; then OUTFILE="$BASE.$INDEX.tmp.ps"; fi

echo $QFLAG

if test $QFLAG; then
  #plot field:
  DLON=$(date_calc "($XLEN|$NLON)")
  DLAT=$(date_calc "($YLEN|(${NLAT}_1))")

  GRDFILE=$BASE.$INDEX.grd

  ORIENT=-ZBLx

  #auto-generate colour palette/vertical levels:
  if test $WCFLAG; then
    echo "extract_field -x $NLON -y $NLAT $EXT_OPT $INFILE $INDEX | gen_zgrid $ZOPTS > $ZGRIDFILE"
    extract_field -x $NLON -y $NLAT $EXT_OPT $INFILE $INDEX | gen_zgrid $ZOPTS > $ZGRIDFILE
    echo "makecpt -T$ZGRIDFILE > $PALETTE"
    makecpt -T$ZGRIDFILE > $PALETTE
  fi

  echo "extract_field -x $NLON -y $NLAT $EXT_OPT $INFILE $INDEX | quick_plot.sh -q -c $PALETTE $PFOPTS $OUTFILE"
  extract_field -x $NLON -y $NLAT $EXT_OPT $INFILE $INDEX | quick_plot.sh -q -c $PALETTE $PFOPTS $OUTFILE
else
  #plot single contour:
  echo "bev2xy ${INFILE} ${INDEX} | quick_plot.sh $PFOPTS $OUTFILE"
  bev2xy ${INFILE} ${INDEX} | quick_plot.sh $PFOPTS $OUTFILE
fi

#clean up:
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

