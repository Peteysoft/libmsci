#!/bin/bash

DATASET="shuttle segment sat pendigits usps vehicle"

#PATH=/cygdrive/f/peteysoft/libmsci/libagf/examples/orthogonal
PATH=.

SUFFIX=.txt

BINPATH=../statlog

# col  value
# 0    1 versus 1 test time user (s)
# 1    1 versus 1 test time system
# 2    ECC test time user
# 3    ECC test time system
# 4    orthogonal test time user
# 5    orthogonal test time system
# 6    1 versus 1 accuracy
# 7    1 versus 1 U.C.
# 8    ECC accuracy
# 9    ECC U.C.
# 10   orthogonal accuracy
# 11   orthogonal U.C.
# 12   1 versus 1 Brier score
# 13   ECC Brier score
# 14   orthogonal Brier score

for D in $DATASET; do
  #printf "%10s" $D
   #echo "$BINPATH/sum_col 0 1 / 7 / 12 / 2 3 / 9 / 13 / 4 / 5 / 11 / 14 < $PATH/$D$SUFFIX | $BINPATH/calc_stats 9 | ./write_results 9 $D"
   #$BINPATH/sum_col 0 1 / 7 / 12 / 2 3 / 9 / 13 / 4 5 / 11 / 14 < $PATH/$D$SUFFIX | $BINPATH/calc_stats 9 | ./write_results 9 $D
   $BINPATH/sum_col 0 1 / 5 / 8 / 2 3 / 7 / 9 < $PATH/$D$SUFFIX | $BINPATH/calc_stats 6 | ./write_results 6 $D
done

