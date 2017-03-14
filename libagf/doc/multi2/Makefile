BIBFILE = ../agf_bib.bib

SOURCE1 = binary_class_map.eps config1.eps rhist.eps sample_roc.eps skillvsr0.eps

SOURCE2 = landclasstree.eps

all:
	make TARGET=constrained SOURCE="$(SOURCE1)" constrained.pdf
	make TARGET=multiclass_review SOURCE="$(SOURCE2)" multiclass_review.pdf
	make TARGET=calibration SOURCE="$(SOURCE1)" calibration.pdf

clean:
	make TARGET=constrained clean_target
	make TARGET=multiclass_review clean_target
	make TARGET=calibration clean_target

clean_target:
	rm -f $(TARGET).dvi
	rm -f $(TARGET).aux
	rm -f $(TARGET).bbl
	rm -f $(TARGET).pdf

$(TARGET).pdf: $(TARGET).dvi
	dvipdf $(TARGET)

$(TARGET).dvi: $(TARGET).tex $(SOURCE) $(BIBFILE)
	latex $(TARGET)
	bibtex $(TARGET)
	latex $(TARGET)
	latex $(TARGET)
