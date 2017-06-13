BIBFILE = ../agf_bib.bib

TARGET = svm_accel_jmlr

FIG = border_vectors.eps \
      support_vectors.eps \
      skill_v_nb.eps \
      skill_v_nt.eps \
      nsv.eps \
	svm_time.eps \
	border_time.eps \
	break_even.eps

SOURCE = summary.tex \
	param.tex \
	results.1.tex \
	results2.1.tex \
	sum_nsv.tex \
	subsample.1.tex \
	shuttle_subsampling.tex \
	verdict.1.tex \
	svm_accel_body.tex \
	svm_accel_abstract.txt \
	svm_accel.bib \
	../agf_bib.bib \
	../pwl.bib

all: $(TARGET).pdf

clean:
	make TARGET=svm_accel_jmlr clean_target

clean_target:
	rm -f $(TARGET).dvi
	rm -f $(TARGET).aux
	rm -f $(TARGET).bbl
	rm -f $(TARGET).pdf

$(TARGET).pdf: $(TARGET).dvi
	dvipdf $(TARGET)

$(TARGET).dvi: $(TARGET).tex $(SOURCE) $(BIBFILE) $(FIG)
	latex $(TARGET)
	bibtex $(TARGET)
	latex $(TARGET)
	latex $(TARGET)

