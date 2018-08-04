#JOURNAL=_elsevier

all: tracer_transport$(JOURNAL).pdf manual.pdf
	make -C proxy2
	make -C pc_proxy

clean:
	make -C proxy2 clean
	make -C pc_proxy clean
	make TARGET=tracer_transport clean_target
	make TARGET=manual clean_target

manual.pdf: description.tex manual.tex grid_example.ps semi_Lagrangian.ps tracer_coord.ps tracer.bib data_flow_general.eps
	latex manual
	latex manual
	dvipdf manual


tracer_transport$(JOURNAL).pdf: tracer_transport$(JOURNAL).tex tracer.bib tt_cover_graphic.ps tt_symbol_list.tex tracer_transport_body.tex
	make TARGET=tracer_transport$(JOURNAL)

$(TARGET).pdf: $(TARGET).tex
	latex $(TARGET)
	bibtex $(TARGET)
	latex $(TARGET)
	latex $(TARGET)
	dvipdf $(TARGET)

clean_target:
	rm -f $(TARGET).aux
	rm -f $(TARGET).bbl
	rm -f $(TARGET).blg
	rm -f $(TARGET).log
	rm -f $(TARGET).pdf
	rm -f $(TARGET).dvi
	rm -f $(TARGET).toc

