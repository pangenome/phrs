PAPER = paper
TEX   = pdflatex
BIB   = bibtex

.PHONY: all clean

all: $(PAPER).pdf

$(PAPER).pdf: $(PAPER).tex bibliography.bib jnl.cls mathphys.bst
	bash compile.sh

clean:
	rm -f $(PAPER).pdf $(PAPER).aux $(PAPER).bbl $(PAPER).blg \
	      $(PAPER).log $(PAPER).out $(PAPER).toc \
	      Meth.aux Meth.bbl Meth.blg \
	      Supp.aux Supp.bbl Supp.blg
