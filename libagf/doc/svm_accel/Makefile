#JOURNAL = _jmlr
#JOURNAL = _springer

TARGET = svm_accel$(JOURNAL)

FIG = border_time.eps \
	border_vectors.eps \
	break_even.eps \
      nsv.eps \
      skill_v_nb.eps \
      skill_v_nt.eps \
      support_vectors.eps \
	svm_time.eps

SOURCE = param.tex \
	results.1.tex \
	results2.1.tex \
	shuttle_subsampling.tex \
	subsample.1.tex \
	sum_nsv.tex \
	summary.tex \
	svm_accel_abstract.txt \
	svm_accel.bib \
	svm_accel_body.tex \
	svm_accel_symbol_def.tex \
	svm_accel_symbol_table.tex \
	verdict.1.tex \
	../agf_bib.bib \
	../pwl.bib

all: $(TARGET).pdf

clean:
	rm -f $(TARGET).dvi
	rm -f $(TARGET).aux
	rm -f $(TARGET).bbl
	rm -f $(TARGET).pdf

$(TARGET).pdf: $(TARGET).dvi
	dvipdf $(TARGET)

$(TARGET).dvi: $(TARGET).tex $(SOURCE) $(FIG)
	latex $(TARGET)
	bibtex $(TARGET)
	latex $(TARGET)
	latex $(TARGET)

