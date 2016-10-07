BIBFILE = ../agf_bib.bib

TARGET = svm_accel

SOURCE = border_vectors.eps skill_v_nb.eps skill_v_nt.eps support_vectors.eps nsv.eps

all: $(TARGET).pdf

clean:
	make TARGET=svm_accel clean_target

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

