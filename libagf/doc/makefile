
all:
	make -C agf2
	make -C multiborders
	make -C directions
	make -C multi2
	make -C svm_accel

clean:
	make -C agf2 clean
	make -C multiborders clean
	make -C directions clean
	make -C multi2 clean
	make -C svm_accel clean

libagf_guide.pdf: libagf_guide.tex agf_bib.bib
	latex libagf_guide.tex
	bibtex libagf_guide
	latex libagf_guide.tex
	latex libagf_guide.tex
	dvipdf libagf_guide
	
