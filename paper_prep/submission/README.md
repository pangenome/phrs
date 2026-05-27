# paper_prep/submission

LaTeX submission package for "Concerted evolution and unorthodox recombination of human subtelomeres".

## How to build

```
cd paper_prep/submission
make
```

Or equivalently:

```
bash compile.sh
```

This runs `pdflatex` + `bibtex` + `pdflatex` + `pdflatex` and produces `paper.pdf`.

## Dependencies

- `texlive-latex-extra` (or equivalent TeX distribution with `natbib`)
- `bibtex` (bibliography processor, part of standard TeX distributions)
- `pdflatex` (from `texlive-latex-base`)
- The `unsrtnat` BibTeX style (shipped with `natbib`; `kpsewhich unsrtnat.bst` should return a path)

No Springer Nature class file (`sn-jnl.cls`) is required: the document falls back to the standard `article` class with `natbib` numerical superscript style (`unsrtnat.bst`), which approximates Nature's reference formatting. For final submission, replace `\documentclass[12pt]{article}` with `\documentclass{sn-jnl}` and install the Springer Nature LaTeX template (available at https://www.springernature.com/gp/authors/campaigns/latex-author-support).

## Contents

| File | Description |
|---|---|
| `paper.tex` | Main LaTeX document |
| `bibliography.bib` | Filtered BibTeX (76 cited entries only) |
| `figures/fig1.pdf` | Main Figure 1 |
| `figures/fig2.pdf` | Main Figure 2 |
| `figures/fig3.pdf` | Main Figure 3 |
| `figures/fig4.pdf` | Main Figure 4 |
| `figures/ed1.pdf` | Extended Data Figure 1 |
| `figures/ed2.pdf` | Extended Data Figure 2 |
| `figures/ed3.pdf` | Extended Data Figure 3 |
| `figures/ed4.pdf` | Extended Data Figure 4 |
| `figures/ed5.pdf` | Extended Data Figure 5 |
| `figures/ed8.pdf` | Extended Data Figure 8 |
| `figures/nj_tree.pdf` | Extended Data NJ tree |
| `Makefile` | Build recipe |
| `compile.sh` | Shell build script |
| `BUILD_LOG.md` | Compile log and validation results |
