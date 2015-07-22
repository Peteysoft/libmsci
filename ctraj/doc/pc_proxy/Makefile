
NAME = PC_proxy
JOURNAL = arXiv

ABSTRACT = $(NAME)_abstract.txt
BIBFILE = $(NAME).bib
ACK = $(NAME)_ack.txt
BODY = $(NAME).tex

WRAPPER = $(NAME)_$(JOURNAL)

INCLUDE = $(ABSTRACT) $(ACK) $(BODY)

TP = latex

BIB = bibtex

FIGURES = lyap_spec.eps meas_loc.eps pc1.eps POAM_proxyscat.eps proxycorr.eps proxyfields.eps \
		proxyO3.0017.stn.eps proxyO3.0017.stn_legend.eps proxyvssonde.eps tcorr.eps tsample.eps

$(WRAPPER).pdf: $(WRAPPER).dvi
	dvipdf $(WRAPPER)

clean:
	rm -f $(WRAPPER).pdf $(WRAPPER).dvi $(WRAPPER).aux $(WRAPPER).bbl \
		$(WRAPPER).blg $(WRAPPER).log

$(WRAPPER).dvi: $(WRAPPER).tex $(INCLUDE) $(FIGURES) $(WRAPPER).aux $(WRAPPER).bbl
	$(TP) $(WRAPPER).tex
	$(TP) $(WRAPPER).tex

$(WRAPPER).bbl: $(WRAPPER).aux $(BIBFILE)
	$(BIB) $(WRAPPER)

$(WRAPPER).aux: $(WRAPPER).tex $(INCLUDE) $(FIGURES) 
	$(TP) $(WRAPPER).tex

