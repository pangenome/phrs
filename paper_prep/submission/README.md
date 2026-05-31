# paper_prep/submission — Build Instructions

Self-contained LaTeX submission for *Concerted evolution and unorthodox recombination of human subtelomeres*, in the Springer Nature `jnl.cls` template.

## TL;DR

```bash
cd paper_prep/submission
make            # or: bash compile.sh
xdg-open paper.pdf
```

Produces `paper.pdf`.

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

`compile.sh` uses `set -euo pipefail`, so the build **aborts on the first
error**. A consequence worth knowing: if pass 1 dies, bibtex and passes 2–3
never run, so the build looks like it has hundreds of "undefined reference"
errors when the real problem is a single fatal error earlier (see
Troubleshooting → Unicode). Always judge ref/citation resolution from the
**final** pass, i.e. `paper.log`, not from the multi-pass `compile.sh` stdout:

```bash
grep -c 'undefined' paper.log     # expect 0
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
  bibliography.bib        — entries cited in main + in Methods
  Makefile                — `make` builds paper.pdf
  compile.sh              — pdflatex -> bibtex (main + Meth + Supp) -> pdflatex x 2
  BUILD_LOG.md            — compile log (predates the 5-figure restructure; stale)
  README.md               — this file
  fig/
    MainFigures/          — Fig1*.png … Fig5*.png (slide-derived, multi-panel)
    ExtendedDataFigures/  — ED_Fig1_mouse_zygotene.png (the only ED figure used);
                            ED_Fig1.pdf … ED_Fig7.pdf are leftover and unused
```

## Figure structure and provenance

The manuscript narrative follows the BoG 2026 talk. The **source of truth** for
the figures and the story is the talk deck and its transcript:

- `paper_prep/paper figures for Concerted evolution and unorthodox recombination of human subtelomeres.pptx`
- `paper_prep/Session7-PopulationGenomics.en.srt`

Five main figures plus one Extended Data figure, each built directly from
specific deck slides:

| Figure | Panels (files in `fig/`)                                                            | Deck slides |
|--------|-------------------------------------------------------------------------------------|-------------|
| **Fig 1** | `Fig1a_genomewide.png`, `Fig1b_lengths.png`                                      | 3, 4        |
| **Fig 2** | `Fig2a_pggb_layout.png` (top), `Fig2b_tree_jaccard.png` + `Fig2c_community_jaccard.png` (below) | 6, 10, 12 |
| **Fig 3** | `Fig3a_C1_chr4q.png`+`Fig3a_C1_chr10q.png`; `Fig3b_C2_chr10p.png`+`Fig3b_C2_chr18p.png`; `Fig3c_C11_chr5q.png`+`Fig3c_C11_chr6q.png` | 14, 15, 16 |
| **Fig 4** | `Fig4a_human_scatter.png`, `Fig4b_porec_community.png`                            | 21, 22      |
| **Fig 5** | `Fig5_pedigree_untangle.png`                                                      | 24          |
| **ED 1**  | `ExtendedDataFigures/ED_Fig1_mouse_zygotene.png`                                  | 20          |

> **ASSET STATUS — placeholders.** All figure files above are PNGs extracted
> from the deck slides (some still carry the speaker-note annotations baked into
> the slide). They make the PDF compile and review correctly, but
> **publication-quality vector versions must be regenerated** from the per-figure
> scripts under `paper_prep/figures/` (`fig1`…`fig4`, etc.).

The earlier reviewer-driven analyses (within-community heterogeneity,
population/$F_{\mathrm{ST}}$ structure, the full 3D forest plot and controls,
CEPH1463 cross-assembler validation, RPE-1, gene enrichment) are **retained as
text** in the body and Methods but no longer have figures. Only `fig:fig1`–
`fig:fig5` and `fig:ed1` (mouse) exist — do **not** re-introduce `\ref{fig:ed2}`
and friends. The original `ExtendedDataFigures/ED_Fig1.pdf`…`ED_Fig7.pdf` remain
on disk but are unreferenced and can be deleted.

## Bibliography

`bibliography.bib` contains **76 entries**. After the restructure the paper
cites **keys in the main bibliography** (`\cite{...}`) and **in Methods
References** (`\citeMeth{...}`, emitted as a separate `multibib` section in
`unsrt` style); some keys are cited in both. Entries in the `.bib` that are no
longer cited (e.g. `hao2024snul`, `bouquet_ChenCEP164Cilia2025`, which were only
in dropped Extended Data captions) are harmlessly ignored by bibtex.

## Source of truth and editing

The narrative is driven by the deck + transcript listed under *Figure structure
and provenance* above. `paper.tex` is now hand-edited directly to match that
narrative — there is no markdown→LaTeX re-sync step and no workgraph dispatch.
(The previous `paper_prep/synthesis/NATURE_DRAFT_v6.md` content was reorganised
into this 5-figure structure and is superseded as the layout authority.)

Writing style follows Guarracino et al. 2023 (*Nature*, recombination between
heterologous human acrocentric chromosomes; `~/Downloads/papers/RecombAcroChro_*.pdf`):
short bold results-section headers, a measured "we find / we observe" cadence,
and a single-paragraph abstract that lands on the cytogenetic-confirmation note.
To stay consistent with the talk, reviewer-era analyses that no longer have a
figure (within-community heterogeneity, population/$F_{\mathrm{ST}}$ structure,
the 14-test 3D forest and its controls, CEPH1463 cross-assembler validation,
RPE-1) were cut from the body and Methods; recover them from git history if a
reviewer asks.

## Troubleshooting

### `! LaTeX Error: Unicode character ... not set up for use with LaTeX.`

A literal non-ASCII character (commonly `≥`, `×`, `→`) reached `pdflatex`. The
usual source is a `note`/title field in `bibliography.bib` — `mathphys.bst`
prints `note` fields, so any raw Unicode there breaks the build. Replace with a
LaTeX macro: `$\geq$`, `$\times$`, `$\rightarrow$`. Find offenders with:

```bash
grep -nP '[^\x00-\x7F]' bibliography.bib
```

Because of `set -e`, this fatal error in pass 1 makes every later cross-ref look
undefined; fix the Unicode and rebuild before chasing "undefined reference"
warnings.

### `! LaTeX Error: File 'jnl.cls' not found.`

Verify you are running `make` / `compile.sh` from inside `paper_prep/submission/`. Both `jnl.cls` and `mathphys.bst` must sit beside `paper.tex`.

### `Citation 'foo' on page N undefined`

A bibkey is referenced in `paper.tex` but missing from `bibliography.bib`. Add the entry, or fix the typo'd key.

### `Reference 'fig:edN' undefined`

Only `fig:ed1` (the mouse figure) exists; `fig:ed2`…`fig:ed11` were removed when
the Extended Data was reduced to a single figure. Do not reference them.

### `Package multibib Error`

`multibib` requires running `bibtex paper`, `bibtex Meth`, and (optionally) `bibtex Supp` in addition to the main bibtex pass. `compile.sh` handles this; if you run bibtex manually, do all three.

### Build hangs at `pdflatex`

A LaTeX error is waiting for input. The compile script uses `-interaction=nonstopmode`, so hangs should not occur in normal use. If they do, read the last 30 lines of `paper.log`.

### `\S` produces "Incomplete \iffalse" error

Known interaction between `jnl.cls`, `multibib`, and `natbib`. The source uses `\textsection` instead of `\S` (renders identically as §). Do not revert this.

## Author block status

The `\author*[...]` block lists three authors: **Andrea Guarracino** (TGen,
co-corresponding), **Angela Gyamfi** (St. Jude), and **Erik Garrison** (UTHSC,
co-corresponding) — matching the BoG talk credits.

- `\subsection*{Acknowledgments}` — **filled** (Heather Mefford; Rob Williams,
  Pjotr Prins, Vincenza Colonna; HPRC Pangenomes Working Group).
- `\subsection*{Author Contributions}` and `\subsection*{Competing Interests}`
  remain TODO stubs (not fabricated; confirm CRediT roles with the authors).

## Journal-portfolio switching

`jnl.cls` ships several Springer Nature reference styles. To switch (e.g., for a non-Nature SN journal), change the `\documentclass` line of `paper.tex`:

```latex
\documentclass[pdflatex,mathphys]{jnl}       % current (Math & Phys Sci)
% \documentclass[pdflatex,standardnature]{jnl}  % Standard Nature Portfolio
% \documentclass[pdflatex,vancouver]{jnl}       % Vancouver
% \documentclass[pdflatex,apa]{jnl}             % APA
% \documentclass[pdflatex,chicago]{jnl}         % Chicago
% \documentclass[pdflatex,aps]{jnl}             % American Physical Society
% \documentclass[pdflatex,basic]{jnl}           % Basic SN / Chemistry
```
