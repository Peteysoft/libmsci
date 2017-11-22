DATASET="shuttle segment sat pendigits usps vehicle"

BINPATH=../statlog

# col  value
# 1    1 versus 1 test time user (s)
# 2    1 versus 1 test time system
# 3    orthogonal test time user
# 4    orthogonal test time system
# 5    1 versus 1 accuracy
# 6    1 versus 1 U.C.
# 7    orthogonal accuracy
# 8    orthogonal U.C.
# 9    1 versus 1 prob. corr.
# 10   1 versus 1 prob. slope
# 11   orthogonal prob. corr.
# 12   orthogonal prob. slope

for D in $DATASET; do
  #printf "%10s" $D
   $BINPATH/sum_col 0 1 / 4 / 5 / 8 / 9 / 2 3 / 6 / 7 / 10 / 11 < $D.txt | $BINPATH/calc_stats 10 | ./write_results $D
done

