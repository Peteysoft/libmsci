agf_directions.pdf: agf_directions.dvi
	dvipdf agf_directions

agf_directions.dvi: agf_directions.tex optimal_agf.tex continuum_classification.tex ../agf_bib.bib
	latex agf_directions.tex
	bibtex agf_directions
	latex agf_directions.tex
	latex agf_directions.tex

clean:
	rm -f agf_directions.aux
	rm -f agf_directions.blg
	rm -f agf_directions.bbl
	rm -f agf_directions.log
	rm -f agf_directions.dvi
	rm -f agf_directions.pdf
