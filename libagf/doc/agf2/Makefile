BASE_NAME = agf_paper
SOURCE = $(BASE_NAME).tex agf_class_comp_body.tex adgaf_abstract.tex
FIGURES = sample_class_bw.ps pdf1_agf.eps pdf2_agf.eps con_knn.eps con_agf.eps con_agf_bord.eps con_svm.eps lvq_vs_brd2.eps
BIBFILE = ../agf_bib.bib

TP = latex

all: $(BASE_NAME).pdf 

$(BASE_NAME).pdf: $(BASE_NAME).dvi $(FIGURES)
	dvipdf $(BASE_NAME).dvi

view: $(BASE_NAME).pdf
	gv $(BASE_NAME).pdf &

print: $(BASE_NAME).ps
	lpr $(BASE_NAME).ps

$(BASE_NAME).ps: $(BASE_NAME).dvi
	dvips $(BASE_NAME).dvi

clean:
	rm -f $(BASE_NAME).pdf $(BASE_NAME).dvi $(BASE_NAME).aux $(BASE_NAME).bbl $(BASE_NAME).blg $(BASE_NAME).toc $(BASE_NAME).log

$(BASE_NAME).dvi: $(BASE_NAME).aux $(SOURCE)
	$(TP) $(BASE_NAME).tex
	$(TP) $(BASE_NAME).tex

$(BASE_NAME).aux: $(SOURCE) $(BIBFILE)
	$(TP) $(BASE_NAME).tex
	bibtex $(BASE_NAME)

