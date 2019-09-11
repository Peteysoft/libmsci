#!/bin/bash

DATASET="pendigits sat segment shuttle urban usps vehicle"

#METHOD=lin
#METHOD=svm
METHOD=acc

#PATH=/cygdrive/f/peteysoft/libmsci/libagf/examples/orthogonal
PATH=.

SUFFIX=_$METHOD.txt

BINPATH=../statlog

# col  value
# 0    1 versus 1 test time user (s)
# 1    1 versus 1 test time system
# 2    1 versus 1 accuracy
# 3    1 versus 1 U.C.
# 4    1 versus 1 Brier score
# 5    1 versus rest time user (s)
# 6    1 versus rest time system
# 7    1 versus rest accuracy
# 8    1 versus rest U.C.
# 9    1 versus rest Brier score
# 10   ECC test time user
# 11   ECC test time system
# 12   ECC accuracy
# 13   ECC U.C.
# 14   ECC Brier score
# 15   orthogonal test time user
# 16   orthogonal test time system
# 17   orthogonal accuracy
# 18   orthogonal U.C.
# 19   orthogonal Brier score
# 20   orthogonal "non-strict" test time user
# 21   orthogonal "non-strict" test time system
# 22   orthogonal "non-strict" accuracy
# 23   orthogonal "non-strict" U.C.
# 24   orthogonal "non-strict" Brier score

for D in $DATASET; do
  #printf "%10s" $D
   #echo "$BINPATH/sum_col 0 1 / 7 / 12 / 2 3 / 9 / 13 / 4 / 5 / 11 / 14 < $PATH/$D$SUFFIX | $BINPATH/calc_stats 9 | ./write_results 9 $D"
   #$BINPATH/sum_col 0 1 / 7 / 12 / 2 3 / 9 / 13 / 4 5 / 11 / 14 < $PATH/$D$SUFFIX | $BINPATH/calc_stats 9 | ./write_results 9 $D
   $BINPATH/sum_col 0 1 / 3 / 4 / 5 6 / 8 / 9 / 10 11 / 13 / 14 / 15 16 / 18 / 19 / 20 21 / 23 / 24 < $PATH/$D$SUFFIX | $BINPATH/calc_stats | ./write_results 15 3 $D
done

