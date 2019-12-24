#!/bin/bash

set -e
set -e

LEVUNIT="mb"
DT="24:0"
#OPTS2="-R"

function extract_ncep_uv {
  extract_ncep_field -p $DATAPATH $TFLAG $UVAR $1 $2
  extract_ncep_field -p $DATAPATH $TFLAG $VVAR $1 $2
}

function extract_pp_uv {
  FIELDS=($(date_calc -f "%y %m %d %H" $DATE))
  YEAR=$(date_calc -f "%Y" $DATE)
  YR=${FIELDS[0]}
  MONTH=${FIELDS[1]}
  DAY=${FIELDS[2]}
  HOUR=${FIELDS[3]}
  FNAME=$(ls $DATAPATH/$YEAR/ppassm_oper?_y${YR}_m${MONTH}_d${DAY}_h${HOUR}.pp)
  read_pp $TFLAG "$FNAME" $1 $2
}

while getopts 'V:0:N:i:f:n:r:p:B:TH' ARG; do
  case $ARG in
    V) SOURCE=$OPTARG
       case $SOURCE in
         0) DATAPATH="/home/peter/data/badc_data"
            DT="24:0"
		 ;;
         1) DATAPATH="/home/peter/data/ncep_data"
            DT="6:0"
		 ;;
        esac
      ;;
    h)  DT=$OPTARG
      ;;
    N)  OPTS="$OPTS -N $OPTARG"
        N=$OPTARG
      ;;
    i)  OPTS="$OPTS -i $OPTARG"
        T0=$OPTARG
      ;;
    f)  TF=$OPTARG
      ;;
    n)  OPTS="$OPTS -n $OPTARG"
      ;;
    r)  OPTS="$OPTS -r $OPTARG"
      ;;
    x) NX=$OPTARG
      ;;
    y) NY=$OPTARG
      ;;
    p) DATAPATH="$OPTARG"
      ;;
    B)  OPTS2="$OPTS2 -B $OPTARG"
      ;;
    T)  TFLAG="-T"
        LEVUNIT="K"
      ;;
    H) echo "import_vfield.sh [-V source]"
       echo "            [-i t0] [-f tf] [-n ngrid] [-r sidelength/2] [-N nt]"
       echo "            [-h dt] [-p datapath] [-T] [-B pagesize] [-I initfile]"
       echo "            level Sfile Nfile"
       echo
       echo "  options:"
       echo "    -V: 0 = UKMO data in 'pp' format from the BADC"
       echo "        1 = NCEP data in NetCDF format"
       echo "          or specify a script or command"
       exit 0
      ;;

  esac
done

shift $(($OPTIND -1))

ZLEV=$1
OPTS="$OPTS -z $ZLEV"
NFILE=$2
SFILE=$3

if test -n $TF
then
  if test -z $N
  then
    echo "date_calc \"(${TF}_${T0})|$DT\""
    N=$(date_calc "(${TF}_${T0})|$DT")
  fi
fi

if test -z $TF
then
  TF=$(date_calc "$T0+$N*$DT")
fi
OPTS="$OPTS -f $TF"

echo DATAPATH=$DATAPATH

DATE=$(date_calc $T0)

case $SOURCE in
  0) FIELDS=($(date_calc -f "%y %m %d %H" $DATE))
    YEAR=$(date_calc -f "%Y" $DATE)
    YR=${FIELDS[0]}
    MONTH=${FIELDS[1]}
    DAY=${FIELDS[2]}
    HOUR=${FIELDS[3]}
    FNAME=$(ls $DATAPATH/$YEAR/ppassm_oper?_y${YR}_m${MONTH}_d${DAY}_h${HOUR}.pp)
    HEADER_FILE=pp${YEAR}$MONTH$DAY${HOUR}_headers.txt
    read_pp -Q $TFLAG "$FNAME" > $HEADER_FILE
    FIELDS=($(grep -m 1 "[ \t]* [0-9]*\/[0-9]*\/[0-9]*.*[ \t]1[ \t]" $HEADER_FILE))
    NX=${FIELDS[5]}
    NY=${FIELDS[4]}
    COMMAND=extract_pp_uv
  ;;
  1) extract_ncep_field -H $TFLAG -p $DATAPATH $UVAR $ZLEV $DATE > $UBASE.txt
    UVAR=uwnd
    VVAR=vwnd
    HEADER=$(head -n 1 $UBASE.txt)
    NX=${HEADER%[ \t]*[^ \t]}
    NY=${HEADER#[^ \t]*[ \t]}
    echo "$NX $NY"
    COMMAND=extract_ncep_uv
  ;;
  *) COMMAND=$SOURCE
  ;;
esac

OPTS="$OPTS -x $NX -y $NY"
  
echo "$COMMAND $ZLEV $T0 | vdata2ds $OPTS $OPTS2 0 $SFILE $NFILE"
$COMMAND $ZLEV $T0 | vdata2ds $OPTS $OPTS2 0 $SFILE $NFILE

for ((I=1; I<=$N; I++))
do
  DATE=$(date_calc "$T0+$I*$DT")
  echo "$COMMAND $ZLEV $DATE | vdata2ds $OPTS2 $I $SFILE $NFILE"
  $COMMAND $ZLEV $DATE | vdata2ds $OPTS2 $I $SFILE $NFILE
done

