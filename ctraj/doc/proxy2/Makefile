TARGET = proxy2_springer

DEPENDENTS = proxy2_abstract.txt proxy2.bib proxy2_body.tex proxy2_ack.txt proxy2_apdx.tex proxy2_linearity2.tex

FIGURES =    classic_O3.1998.10.01.ps \
	     firstpc1997.12.01.ps \
	     global_classic_sonde_error.eps \
	     global_classic_sonde_Nhemi.eps \
	     global_classic_sonde_scatter.eps \
	     global_classic_xval_error.eps \
	     global_classic_xval_scatter.eps \
	     global_PC_sonde_error.eps \
	     global_PC_sonde_Nhemi.eps \
	     global_PC_sonde_scatter.eps \
	     global_PC_xval_error.eps \
	     global_PC_xval_scatter.eps \
	     lyap_spec.eps \
	     Nhemi_classic_sonde.eps \
	     Nhemi_PC_sonde.eps \
	     passive1997_12_01.eps \
	     pc_O3_1998_10_01.eps \
	     ozone_sonde_stations.ps \
	     sample_sv.ps \
	     tracer_correlation.eps \
	     POAM_loc.eps
#	     svl1.ps svl2.ps svl3.ps svl4.ps \
#	     svr1.ps svr2.ps svr3.ps svr4.ps 

$(TARGET).pdf: $(TARGET).tex $(DEPENDENTS) $(FIGURES)
	latex $(TARGET)
	bibtex $(TARGET)
	latex $(TARGET)
	latex $(TARGET)
	dvipdf $(TARGET)

clean:
	rm -f $(TARGET).aux
	rm -f $(TARGET).bbl
	rm -f $(TARGET).blg
	rm -f $(TARGET).log
	rm -f $(TARGET).pdf
	rm -f $(TARGET).dvi
	rm -f $(TARGET).toc
