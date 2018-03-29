#!/bin/bash

DATASET="shuttle segment sat pendigits usps vehicle"

PATH=/cygdrive/f/peteysoft/libmsci/libagf/examples/orthogonal

SUFFIX=.txt

BINPATH=../statlog

# col  value
# 0    1 versus 1 test time user (s)
# 1    1 versus 1 test time system
# 2    orthogonal test time user
# 3    orthogonal test time system
# 4    1 versus 1 accuracy
# 5    1 versus 1 U.C.
# 6    orthogonal accuracy
# 7    orthogonal U.C.
# 8    1 versus 1 Brier score
# 9    orthogonal Brier score

for D in $DATASET; do
  #printf "%10s" $D
   echo "$BINPATH/sum_col 0 1 / 5 / 8 / 2 3 / 7 / 9 < $PATH/$D$SUFFIX | $BINPATH/calc_stats 6 | ./write_results 6 $D"
   $BINPATH/sum_col 0 1 / 5 / 8 / 2 3 / 7 / 9 < $PATH/$D$SUFFIX | $BINPATH/calc_stats 6 | ./write_results 6 $D
done

