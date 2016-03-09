all: agf2/agf_paper.pdf multiborders/multiborders_arxiv.pdf directions/agf_directions.pdf libagf_guide.pdf

clean:
	make -C agf2 clean
	make -C multiborders clean
	make -C directions clean

agf2/agf_paper.pdf:
	make -C agf2

multiborders/multiborders_arxiv.pdf:
	make -C multiborders

directions/agf_directions.pdf:
	make -C directions

libagf_guide.pdf: libagf_guide.tex agf_bib.bib
	latex libagf_guide.tex
	bibtex libagf_guide
	latex libagf_guide.tex
	latex libagf_guide.tex
	dvipdf libagf_guide
	
