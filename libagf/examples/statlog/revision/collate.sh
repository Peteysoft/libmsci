#DATASET="heart shuttle sat segment dna splice codrna letter pendigits"
DATASET="usps mnist ijcnn1 madelon seismic mushrooms phishing humidity"
#DATASET="seismic"

# col  value
# 0    Lin train time user (s)
# 1    Lin train time system
# 2    Lin test time user
# 3    Lin test time system
# 4    SVM train time user
# 5    SVM train time system
# 6    SVM test time user
# 7    SVM test time system
# 8    ACC train time user
# 9    ACC train time system
# 10   ACC test time user
# 11   ACC test time system
# 12   Lin accuracy
# 13   Lin U.C.
# 14   SVM accuracy
# 15   SVM U.C.
# 16   ACC accuracy
# 17   ACC U.C.
# 18   number of suport vectors

for D in $DATASET; do
  #printf "%10s" $D
   ../sum_col 0 1 / 2 3 / 4 5 / 6 7 / 8 9 / 10 11 / 12 / 13 / 14 / 15 / 16 / 17 / 18 < $D.txt | ../calc_stats 10 | ./write_results $D
   echo "\\hline"
done

