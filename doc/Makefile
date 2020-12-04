TARGET := doc

all:
	@pdflatex $(TARGET)
	@bibtex $(TARGET)
	@pdflatex $(TARGET)
	@pdflatex $(TARGET)

clean:
	$(RM) $(TARGET).{aux,bbl,blg,dvi,log,out,pdf}

.PHONY: all clean
