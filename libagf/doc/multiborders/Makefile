BASE_NAME = multiborders

JOURNAL = arxiv

WRAPPER = $(BASE_NAME)_$(JOURNAL)

all: $(WRAPPER).pdf

clean:
	rm -f $(WRAPPER).aux
	rm -f $(WRAPPER).log
	rm -f $(WRAPPER).bbl
	rm -f $(WRAPPER).blg
	rm -f $(WRAPPER).pdf

multiborders_arxiv.pdf: $(WRAPPER).tex $(BASE_NAME)_body.tex $(BASE_NAME)_abstract.tex ../agf_bib.bib multiclass_figure2.pdf
	pdflatex $(WRAPPER)
	bibtex $(WRAPPER)
	pdflatex $(WRAPPER)
	pdflatex $(WRAPPER)

