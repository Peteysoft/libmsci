BASE=shuttle
#BASE=usps

MODEL=${BASE}_svm/${BASE}.ortho.mbc
CLASS=${BASE}_svm/${BASE}.tst
PROB=${BASE}_svm/${BASE}.ortho.txt

SUB="10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000"

FRAC="0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1"

for F in $FRAC; do
	validate_probabilities -H bar.cls foobar.txt
	echo "$F " | tr -d '\n'
	agf_preprocess -z -f $F $CLASS foo bar
	classify_m -Q 8 $MODEL bar.vec foobar > foobar.txt
	validate_probabilities -H bar.cls foobar.txt | grep -E -o "[0-9]+\.?[0-9]*" | tr '\n' ' '
	echo
done

