#JOURNAL = _jmlr
JOURNAL = _springer

TARGET = svm_accel$(JOURNAL)

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
	svm_accel_symbol_def.tex \
	svm_accel_symbol_table.tex \
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

