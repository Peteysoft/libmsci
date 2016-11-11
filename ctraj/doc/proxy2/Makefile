TARGET = proxy2

$(TARGET).pdf: $(TARGET).tex
	latex $(TARGET)
	bibtex $(TARGET)
	latex $(TARGET)
	latex $(TARGET)
	dvipdf $(TARGET)

