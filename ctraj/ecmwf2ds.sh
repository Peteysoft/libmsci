#!/bin/bash

set -e

LEVTYPE="pl"
DT="6:0"
OPTS2="-R"

UVAR=131.128
VVAR=132.128

while getopts '0:N:i:f:n:r:p:B:TI:OH' ARG; do
  case $ARG in
    h)  DT=$OPTARG
      ;;
    N)  OPTS="$OPTS -N $OPTARG"
        N=$OPTARG
      ;;
    i)  OPTS="$OPTS -i $OPTARG"
        T0=$OPTARG
      ;;
    I)  ECMWF_OPTS="$ECWMF_OPTS --init $OPTARG"
      ;;
    f)  OPTS="$OPTS -f $OPTARG"
        TF=$OPTARG
      ;;
    n)  OPTS="$OPTS -n $OPTARG"
      ;;
    r)  OPTS="$OPTS -r $OPTARG"
      ;;
    p)  DATAPATH=$OPTARG
      ;;
    B)  OPTS2="$OPTS2 -B $OPTARG"
      ;;
    T)  LEVTYPE="pt"
      ;;
    O)  ECMWF_OPTS="$ECMWF_OPTS --old"
      ;;
    H) echo "ecmwf2ds.sh [-i t0] [-f tf] [-n ngrid] [-r sidelength/2] [-N nt] [-O]"
       echo "            [-h dt] [-p datapath] [-T] [-B pagesize] [-I initfile]"
       echo "            level Sfile Nfile"
       echo
       echo "- requires wgrib and an account on the ECMWF server"
       echo "  (data is downloaded directly from the ECMWF)"
       echo "- if datapath is specified, saves grib files for later use"
       echo "- calls vdata2ds (most of the parameters are merely passed to it)"
       echo "- initfile contains account information for ECWMF data servers"
       echo "- -O=use old data server/account"
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
DATE1=$(date_calc -f "%Y%m%d%H" $T0)
#this one used to call wgrib: 
DATE2=$(date_calc -f "%Y%m%d %H" $T0)
  
UBASE=$DATAPATH/$DATE1.$UVAR.$LEVTYPE.$ZLEV
VBASE=$DATAPATH/$DATE1.$VVAR.$LEVTYPE.$ZLEV

echo "get_ecmwf.py $ECMWF_OPTS $DATE2 $LEVTYPE $ZLEV $UVAR $UBASE.grib"
echo "get_ecmwf.py $ECMWF_OPTS $DATE2 $LEVTYPE $ZLEV $VVAR $VBASE.grib"

if test -z $DFLAG
then
  if test -f $UBASE.grib
  then
    echo
    echo fuck2
  else
    echo fuck3
    get_ecmwf.py $ECMWF_OPTS $DATE2 $LEVTYPE $ZLEV $UVAR $UBASE.grib
  fi
  if test -f $VBASE.grib
  then
    echo
  else
    get_ecmwf.py $ECMWF_OPTS $DATE2 $LEVTYPE $ZLEV $VVAR $VBASE.grib
  fi
else
  echo fuck1
  get_ecmwf.py $ECMWF_OPTS $DATE2 $LEVTYPE $ZLEV $UVAR $UBASE.grib
  get_ecmwf.py $ECMWF_OPTS $DATE2 $LEVTYPE $ZLEV $VVAR $VBASE.grib
fi

wgrib -d all -h -text -o $UBASE.txt $UBASE.grib
HEADER=$(head -n 1 $UBASE.txt)
NX=${HEADER%[ \t]*[^ \t]}
NY=${HEADER#[^ \t]*[ \t]}
echo "$NX $NY"
OPTS="$OPTS -x $NX -y $NY"

echo "wgrib -d all -nh -text -o $UBASE.txt $UBASE.grib"
echo "wgrib -d all -nh -text -o $VBASE.txt $VBASE.grib"

wgrib -d all -nh -text -o $UBASE.txt $UBASE.grib
wgrib -d all -nh -text -o $VBASE.txt $VBASE.grib

echo "cat $UBASE.txt $VBASE.txt | vdata2ds $OPTS $OPTS2 $I $SFILE $NFILE"
cat $UBASE.txt $VBASE.txt | vdata2ds $OPTS $OPTS2 0 $SFILE $NFILE

#I know there is a way to do this sensibly
#bash syntax for conditional expressions suck so badly I can't figure
#out what it is tho...
if test -z $DFLAG
then
  echo
else
  rm $UBASE.grib
  rm $VBASE.grib
fi

rm $UBASE.txt
rm $VBASE.txt
    
for ((I=1; I<=$N; I++))
do
  #this one used for filenames:
  DATE1=$(date_calc -f "%Y%m%d%H" "$T0+$I*$DT")
  #this one used to call wgrib: 
  DATE2=$(date_calc -f "%Y%m%d %H" "$T0+$I*$DT")
  
  UBASE=$DATAPATH/$DATE1.$UVAR.$LEVTYPE.$ZLEV
  VBASE=$DATAPATH/$DATE1.$VVAR.$LEVTYPE.$ZLEV

  echo "get_ecmwf.py $ECMWF_OPTS $DATE2 $LEVTYPE $ZLEV $UVAR $UBASE.grib"
  echo "get_ecmwf.py $ECMWF_OPTS $DATE2 $LEVTYPE $ZLEV $VVAR $VBASE.grib"

  if test -z $DFLAG
  then
    if test -f $UBASE.grib
    then
      echo #if only bash had a sensible syntax for conditional expressions...
    else
      get_ecmwf.py $ECMWF_OPTS $DATE2 $LEVTYPE $ZLEV $UVAR $UBASE.grib
    fi
    if test -f $VBASE.grib
    then
      echo #if only bash had a sensible syntax for conditional expressions...
    else
      get_ecmwf.py $ECMWF_OPTS $DATE2 $LEVTYPE $ZLEV $VVAR $VBASE.grib
    fi
  else
    get_ecmwf.py $ECMWF_OPTS $DATE2 $LEVTYPE $ZLEV $UVAR $UBASE.grib
    get_ecmwf.py $ECMWF_OPTS $DATE2 $LEVTYPE $ZLEV $VVAR $VBASE.grib
  fi

  echo "wgrib -d all -nh -text -o $UBASE.txt $UBASE.grib"
  echo "wgrib -d all -nh -text -o $VBASE.txt $VBASE.grib"

  wgrib -d all -nh -text -o $UBASE.txt $UBASE.grib
  wgrib -d all -nh -text -o $VBASE.txt $VBASE.grib

  echo "cat $UBASE.txt $VBASE.txt | vdata2ds $OPTS2 $I $SFILE $NFILE"
  cat $UBASE.txt $VBASE.txt | vdata2ds $OPTS2 $I $SFILE $NFILE

  if test -z $DFLAG
  then
    echo
  else
    rm $UBASE.grib
    rm $VBASE.grib
  fi

  rm $UBASE.txt
  rm $VBASE.txt

done

