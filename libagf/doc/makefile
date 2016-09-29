
all:
	make -C agf2
	make -C multiborders
	make -C directions
	make -C multi2

clean:
	make -C agf2 clean
	make -C multiborders clean
	make -C directions clean

all:
	make -C agf2
	make -C multiborders
	make -C directions
	make -C multi2

libagf_guide.pdf: libagf_guide.tex agf_bib.bib
	latex libagf_guide.tex
	bibtex libagf_guide
	latex libagf_guide.tex
	latex libagf_guide.tex
	dvipdf libagf_guide
	
