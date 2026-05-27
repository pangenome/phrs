# paper_prep/submission — Build Instructions

Self-contained LaTeX submission for *Concerted evolution and unorthodox recombination of human subtelomeres*, in the Springer Nature `jnl.cls` template.

## TL;DR

```bash
cd paper_prep/submission
make            # or: bash compile.sh
xdg-open paper.pdf
```

Produces `paper.pdf` (~37 pages, ~2 MB).

## Dependencies

| Tool          | Purpose                                    |
|---------------|--------------------------------------------|
| `pdflatex`    | LaTeX engine (TeX Live 2020 or newer)      |
| `bibtex`      | Bibliography processing                    |
| `jnl.cls`     | Springer Nature class (bundled here)       |
| `mathphys.bst`| Springer Nature math/phys bib style (bundled) |

`jnl.cls` and `mathphys.bst` are bundled in this directory, so the build is self-contained. No system-wide install of the Springer Nature template is needed.

### Install pdflatex + bibtex

Ubuntu / Debian:

```bash
sudo apt-get install texlive-latex-base texlive-latex-extra texlive-fonts-recommended texlive-bibtex-extra
```

Fedora / RHEL:

```bash
sudo dnf install texlive-scheme-medium
```

macOS (Homebrew):

```bash
brew install --cask mactex-no-gui
```

Verify:

```bash
pdflatex --version
bibtex --version
```

## Build

### Option 1 — `make`

```bash
cd paper_prep/submission
make
```

`make` runs `compile.sh`, which executes:

```
pdflatex paper      # pass 1: write aux files
bibtex   paper      # main bibliography
bibtex   Meth       # Methods References (multibib)
bibtex   Supp       # Supplementary References (multibib, may be empty)
pdflatex paper      # pass 2: resolve citations
pdflatex paper      # pass 3: resolve cross-references
```

### Option 2 — direct script

```bash
cd paper_prep/submission
bash compile.sh
```

Identical to `make`.

### Clean

```bash
make clean
```

Removes `paper.pdf` and all `.aux`, `.bbl`, `.blg`, `.log`, `.out`, `.toc` files (including the `Meth.*` and `Supp.*` multibib auxiliaries).

## Output

| File         | Description                                       |
|--------------|---------------------------------------------------|
| `paper.pdf`  | Final compiled submission                         |
| `paper.log`  | Full LaTeX log (read this on build errors)        |
| `paper.blg`  | BibTeX log for the main bibliography              |
| `Meth.blg`   | BibTeX log for the Methods References bibliography|

## Directory layout

```
paper_prep/submission/
  paper.tex               — single self-contained LaTeX source
  jnl.cls                 — Springer Nature class (bundled)
  mathphys.bst            — bibliography style (bundled)
  bibliography.bib        — 76 cited references (filtered from REFERENCES_v6.bib)
  Makefile                — `make` builds paper.pdf
  compile.sh              — pdflatex -> bibtex (main + Meth + Supp) -> pdflatex x 2
  BUILD_LOG.md            — compile log + 10-check validation evidence
  README.md               — this file
  fig/
    MainFigures/          — Figure1.pdf ... Figure4.pdf
    ExtendedDataFigures/  — ED_Fig1.pdf ... ED_Fig7.pdf
```

## Figure provenance

| Submission file                          | Source                                                |
|------------------------------------------|-------------------------------------------------------|
| `fig/MainFigures/Figure1.pdf`            | `paper_prep/figures/fig1/figure_fig1.pdf`             |
| `fig/MainFigures/Figure2.pdf`            | `paper_prep/figures/fig2/figure_fig2.pdf`             |
| `fig/MainFigures/Figure3.pdf`            | `paper_prep/figures/fig3/figure_fig3.pdf`             |
| `fig/MainFigures/Figure4.pdf`            | `paper_prep/figures/fig4/figure_fig4.pdf`             |
| `fig/ExtendedDataFigures/ED_Fig1.pdf`    | `paper_prep/figures/ed1/figure_ed1.pdf`               |
| `fig/ExtendedDataFigures/ED_Fig2.pdf`    | `paper_prep/figures/ed2/figure_ed2.pdf`               |
| `fig/ExtendedDataFigures/ED_Fig3.pdf`    | `paper_prep/figures/ed3/figure_ed3.pdf`               |
| `fig/ExtendedDataFigures/ED_Fig4.pdf`    | `paper_prep/figures/ed4/figure_ed4.pdf`               |
| `fig/ExtendedDataFigures/ED_Fig5.pdf`    | `paper_prep/figures/ed5/figure_ed5.pdf`               |
| `fig/ExtendedDataFigures/ED_Fig6.pdf`    | `paper_prep/figures/ed8/figure_ed8.pdf` (renumbered)  |
| `fig/ExtendedDataFigures/ED_Fig7.pdf`    | `paper_prep/figures/nj_tree_arms/nj_tree_annotated.pdf` |

## Bibliography

`bibliography.bib` contains exactly **76 entries**, filtered from `paper_prep/synthesis/REFERENCES_v6.bib` (374 entries total) to match the 76 cited bibkeys listed in `paper_prep/synthesis/RENDERED_REFERENCES_v6.md`.

- Main bibliography: `mathphys` style (Springer Nature numerical superscripts).
- Methods References: emitted by `multibib` into a separate section, `unsrt` style. Methods-only citations use `\citeMeth{...}` (16 keys); main-text citations use `\cite{...}` (60 keys).

## Source of truth

Content is a faithful transcription of `paper_prep/synthesis/NATURE_DRAFT_v6.md` (200-word abstract + 3299-word main + 1591-word Methods). Figure captions are taken from the NATURE_DRAFT_v6.md `## Figure list` section.

To re-sync content after editing the markdown source, dispatch a workgraph task pointing at the updated `NATURE_DRAFT_v6.md`; do **not** hand-edit `paper.tex` unless the change is purely typographic.

## Troubleshooting

### `! LaTeX Error: File 'jnl.cls' not found.`

Verify you are running `make` / `compile.sh` from inside `paper_prep/submission/`. Both `jnl.cls` and `mathphys.bst` must sit beside `paper.tex`.

### `Citation 'foo' on page N undefined`

A bibkey is referenced in `paper.tex` but missing from `bibliography.bib`. Cross-check against `paper_prep/synthesis/RENDERED_REFERENCES_v6.md`. The 76 cited keys are listed there.

### `Package multibib Error`

`multibib` requires running `bibtex paper`, `bibtex Meth`, and (optionally) `bibtex Supp` in addition to the main bibtex pass. `compile.sh` handles this; if you run bibtex manually, do all three.

### Build hangs at `pdflatex`

A LaTeX error is waiting for input. The compile script uses `-interaction=nonstopmode`, so hangs should not occur in normal use. If they do, read the last 30 lines of `paper.log`.

### `\S` produces "Incomplete \iffalse" error

Known interaction between `jnl.cls`, `multibib`, and `natbib`. The source uses `\textsection` instead of `\S` (renders identically as §). Do not revert this.

## Author block status

The current `\author*[...]` block lists Andrea Guarracino (TGen) and Erik Garrison (UTHSC) as co-corresponding authors. A `% TODO: complete co-author list and affiliations` comment marks where the full author and affiliation list should be filled in before submission.

The `\subsection*{Acknowledgments}`, `\subsection*{Author Contributions Statement}`, and `\subsection*{Competing Interests Statement}` blocks are TODO stubs.

## Journal-portfolio switching

`jnl.cls` ships several Springer Nature reference styles. To switch (e.g., for a non-Nature SN journal), change line 8 of `paper.tex`:

```latex
\documentclass[pdflatex,mathphys]{jnl}       % current (Math & Phys Sci)
% \documentclass[pdflatex,standardnature]{jnl}  % Standard Nature Portfolio
% \documentclass[pdflatex,vancouver]{jnl}       % Vancouver
% \documentclass[pdflatex,apa]{jnl}             % APA
% \documentclass[pdflatex,chicago]{jnl}         % Chicago
% \documentclass[pdflatex,aps]{jnl}             % American Physical Society
% \documentclass[pdflatex,basic]{jnl}           % Basic SN / Chemistry
```

See the PGGB template (`/home/guarracino/Downloads/_PGGB__Building_pangenome_graphs/main.tex` lines 18–34) for the full menu.

## Validation evidence

`BUILD_LOG.md` records 10 validation checks: exit-0 compile, no undefined refs, jnl documentclass, jnl.cls and mathphys.bst bundled, 76 bib entries, all 76 bibkeys cited, 11 figure PDFs present, Methods References section emitted, and detex word count check. Open it after every rebuild that changes structure.
