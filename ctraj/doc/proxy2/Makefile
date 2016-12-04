TARGET = proxy2

$(TARGET).pdf: $(TARGET).tex
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
