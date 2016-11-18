#!/bin/bash

set -e

LEVUNIT="mb"
DT="6:0"
#OPTS2="-R"

UVAR=uwnd
VVAR=vwnd

VPATH=.

while getopts '0:N:i:f:n:r:p:B:T' ARG; do
  case $ARG in
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
    p)  VPATH=$OPTARG
      ;;
    B)  OPTS="$OPTS -B $OPTARG"
      ;;
    T)  EXTOPTS="$EXTOPTS -T"
        LEVUNIT="K"
      ;;
    H) echo "test_vfield_conv.sh [-i t0] [-f tf] [-p path] level date
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

if test -z $DATAPATH
then
  DFLAG=1
  DATAPATH=.
fi
echo DFLAG=$DFLAG

#this one used for filenames:
echo "date_calc -f \"%Y%m%d%H\" $T0"
DATE1=$(date_calc -f "%Y%m%d%H" $T0)
#this one used to call extract_ncep_field:
DATE2=$(date_calc $T0)
  
UBASE=$DATAPATH/$DATE1.$UVAR.$ZLEV$LEVUNIT
VBASE=$DATAPATH/$DATE1.$VVAR.$ZLEV$LEVUNIT

extract_ncep_field -H $EXTOPTS $UVAR $ZLEV $T0 > $UBASE.txt
HEADER=$(head -n 1 $UBASE.txt)
NX=${HEADER%[ \t]*[^ \t]}
NY=${HEADER#[^ \t]*[ \t]}
echo "$NX $NY"
OPTS="$OPTS -x $NX -y $NY"

echo "extract_ncep_field $EXTOPTS $UVAR $ZLEV $T0 > $UBASE.txt"
echo "extract_ncep_field $EXTOPTS $VVAR $ZLEV $T0 > $VBASE.txt"

extract_ncep_field $EXTOPTS $UVAR $ZLEV $T0 > $UBASE.txt
extract_ncep_field $EXTOPTS $VVAR $ZLEV $T0 > $VBASE.txt

echo "cat $UBASE.txt $VBASE.txt | vdata2ds $OPTS $OPTS2 0 $SFILE $NFILE"
cat $UBASE.txt $VBASE.txt | vdata2ds $OPTS $OPTS2 0 $SFILE $NFILE

rm $UBASE.txt
rm $VBASE.txt
    
for ((I=1; I<=$N; I++))
do
  #this one used for filenames:
  DATE1=$(date_calc -f "%Y%m%d%H" "$T0+$I*$DT")
  #this one used to call wgrib: 
  DATE2=$(date_calc "$T0+$I*$DT")
  
  UBASE=$DATAPATH/$DATE1.$UVAR.$ZLEV$LEVUNIT
  VBASE=$DATAPATH/$DATE1.$VVAR.$ZLEV$LEVUNIT

  echo "extract_ncep_field $EXTOPTS $UVAR $ZLEV $DATE2 > $UBASE.grib"
  echo "extract_ncep_field $EXTOPTS $VVAR $ZLEV $DATE2 > $VBASE.grib"

  extract_ncep_field $EXTOPTS $UVAR $ZLEV $DATE2 > $UBASE.txt
  extract_ncep_field $EXTOPTS $VVAR $ZLEV $DATE2 > $VBASE.txt

  echo "cat $UBASE.txt $VBASE.txt | vdata2ds $OPTS2 $I $SFILE $NFILE"
  cat $UBASE.txt $VBASE.txt | vdata2ds $OPTS2 $I $SFILE $NFILE

  rm $UBASE.txt
  rm $VBASE.txt

done

