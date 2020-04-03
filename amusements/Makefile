all: quine count
	./quine > quine.copy
	diff quine.c quine.copy

clean: 
	rm -f quine count quine.copy

quine: quine.c
	gcc -o quine quine.c

count: count.cc
	g++ -o count count.cc
