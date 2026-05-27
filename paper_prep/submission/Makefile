.PHONY: all clean

all: paper.pdf

paper.pdf: paper.tex bibliography.bib figures/fig1.pdf figures/fig2.pdf \
           figures/fig3.pdf figures/fig4.pdf \
           figures/ed1.pdf figures/ed2.pdf figures/ed3.pdf \
           figures/ed4.pdf figures/ed5.pdf figures/ed8.pdf \
           figures/nj_tree.pdf
	pdflatex -interaction=nonstopmode -file-line-error paper.tex
	bibtex paper
	pdflatex -interaction=nonstopmode -file-line-error paper.tex
	pdflatex -interaction=nonstopmode -file-line-error paper.tex

clean:
	rm -f paper.aux paper.bbl paper.bcf paper.blg paper.log \
	      paper.out paper.run.xml paper.toc paper.fls paper.fdb_latexmk
