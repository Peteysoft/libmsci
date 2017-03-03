FILES="heart.txt shuttle.txt segment.txt sat.txt"

for F in $FILES; do
   ./sum_col.exe 0 1 / 2 3 / 4 5 / 6 7 / 8 9 / 10 11 / 12 13 / 14 / 15 / 16 / 17 / 18 / 19 / 20 / 21 < $F | ./calc_stats
done

