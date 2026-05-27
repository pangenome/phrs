# paper_prep/submission — Build Instructions

## Quick build

```bash
cd paper_prep/submission
bash compile.sh
```

Or via Make:

```bash
cd paper_prep/submission
make
```

## Output

`paper.pdf` — compiled Springer Nature jnl-class submission.

## Dependencies

| Tool | Purpose |
|------|---------|
| `pdflatex` | LaTeX engine (TeX Live 2020+ recommended) |
| `bibtex` | Bibliography processing |
| `jnl.cls` | Springer Nature consolidated class (bundled) |
| `mathphys.bst` | Math/physical sciences bibliography style (bundled) |

`jnl.cls` and `mathphys.bst` are copied from the PGGB template at
`/home/guarracino/Downloads/_PGGB__Building_pangenome_graphs/` and are
bundled here so the build is self-contained.

## Directory layout

```
paper_prep/submission/
  paper.tex               — single self-contained LaTeX source
  jnl.cls                 — Springer Nature class (bundled)
  mathphys.bst            — bibliography style (bundled)
  bibliography.bib        — 76 cited references (filtered from REFERENCES_v6.bib)
  Makefile                — make target: `make` builds paper.pdf
  compile.sh              — pdflatex -> bibtex (main+Meth+Supp) -> pdflatex x 2
  BUILD_LOG.md            — compile log and validation evidence
  README.md               — this file
  fig/
    MainFigures/          — Figure1.pdf ... Figure4.pdf
    ExtendedDataFigures/  — ED_Fig1.pdf ... ED_Fig7.pdf
    SupplementaryFigures/ — (empty; reserved)
```

## Figure provenance

| File | Source |
|------|--------|
| `fig/MainFigures/Figure1.pdf` | `paper_prep/figures/fig1/figure_fig1.pdf` |
| `fig/MainFigures/Figure2.pdf` | `paper_prep/figures/fig2/figure_fig2.pdf` |
| `fig/MainFigures/Figure3.pdf` | `paper_prep/figures/fig3/figure_fig3.pdf` |
| `fig/MainFigures/Figure4.pdf` | `paper_prep/figures/fig4/figure_fig4.pdf` |
| `fig/ExtendedDataFigures/ED_Fig1.pdf` | `paper_prep/figures/ed1/figure_ed1.pdf` |
| `fig/ExtendedDataFigures/ED_Fig2.pdf` | `paper_prep/figures/ed2/figure_ed2.pdf` |
| `fig/ExtendedDataFigures/ED_Fig3.pdf` | `paper_prep/figures/ed3/figure_ed3.pdf` |
| `fig/ExtendedDataFigures/ED_Fig4.pdf` | `paper_prep/figures/ed4/figure_ed4.pdf` |
| `fig/ExtendedDataFigures/ED_Fig5.pdf` | `paper_prep/figures/ed5/figure_ed5.pdf` |
| `fig/ExtendedDataFigures/ED_Fig6.pdf` | `paper_prep/figures/ed8/figure_ed8.pdf` (renumbered) |
| `fig/ExtendedDataFigures/ED_Fig7.pdf` | `paper_prep/figures/nj_tree_arms/nj_tree_annotated.pdf` |

## Bibliography

`bibliography.bib` contains exactly 76 entries filtered from
`paper_prep/synthesis/REFERENCES_v6.bib` to match the 76 bibkeys listed in
`paper_prep/synthesis/RENDERED_REFERENCES_v6.md`.

The main bibliography uses the `mathphys` style (Springer Nature numerical).
Methods-only references are processed by `multibib` into a separate Methods
References section using the `unsrt` style.
