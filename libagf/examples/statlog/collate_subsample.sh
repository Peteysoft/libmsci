DATASET="shuttle segment sat codrna pendigits usps ijcnn1 mushrooms phishing humidity"

# col  value
# 0    training time user
# 1    training time system
# 2    test time user
# 3    test time system
# 4    accuracy
# 5    U.C.
# 6    training samples
# 7    support vectors

for D in $DATASET; do
  #printf "%10s" $D
   ./sum_col.exe 6 / 7 / 0 1 / 2 3 / 4 / 5  < $D.sub.txt | ./calc_stats 10 | ./write_results2 $D
   echo "\\hline"
done

