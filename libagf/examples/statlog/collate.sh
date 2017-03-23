#DATASET="heart shuttle segment sat letter codrna pendigits dna seismic usps"
DATASET="ijcnn1 leu mushrooms madelon phishing splice humidity mnist"

# col  value
# 1    KNN time user (s)
# 2    KNN time system
# 3    AGF train time user
# 4    AGF train time system
# 5    AGF test time user
# 6    AGF test time system
# 7    SVM train time user
# 8    SVM train time system
# 9    SVM test time user
# 10   SVM test time system
# 11   ACC train time user
# 12   ACC train time system
# 13   ACC test time user
# 14   ACC test time system
# 15   KNN accuracy
# 16   KNN U.C.
# 17   AGF accuracy
# 18   AGF U.C.
# 19   SVM accuracy
# 20   SVM U.C.
# 21   ACC accuracy
# 22   ACC U.C.
# 23   number of suport vectors

for D in $DATASET; do
  #printf "%10s" $D
   ./sum_col.exe 0 1 / 2 3 / 4 5 / 6 7 / 8 9 / 10 11 / 12 13 / 14 / 15 / 16 / 17 / 18 / 19 / 20 / 21 / 22 / 23 / 24 < $D.txt | ./calc_stats 10 | ./write_results $D
   echo "\\hline"
done

