OBJ = $(wildcard *.tex)

doc: ${OBJ}
	@xelatex -synctex=1 $^
	make clean

all: ${OBJ}
	@xelatex -synctex=1 $^
	@bibtex main.aux
	@xelatex -synctex=1 $^
	@xelatex -synctex=1 $^
	make clean

.PHONY: clean
clean:
	rm -f *.aux *.log *.out *.toc *.bbl *.run.xml *.blg 
