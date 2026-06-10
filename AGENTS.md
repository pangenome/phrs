# CLAUDE.md

## Project: subtelomeric PHR analysis

Research repo, not a software product. Inter-chromosomal subtelomeric sequence
sharing across 233 HPRCv2 samples (465 haplotypes, 15,668 PHRs, 41 signal-bearing
arms, 15 arm-level / 50 sequence-level communities), plus 3D-genome validation
(Hi-C, Pore-C, CiFi, Dip-C, sperm sc, RPE-1, mouse meiosis). Outputs are a
14-section analysis report, paper figures, and slide decks.

GitHub: `https://github.com/ekg/phrs`.

## Where data lives

Source-of-truth bioinformatics data is NOT in the repo. Scripts read absolute
paths:

- `/moosefs/guarracino/HPRCv2/PHR_III/` — similarity matrices, Leiden/UPGMA
  community assignments, 3D-validation tables, enrichment outputs.
- `/moosefs/erikg/phrs/` — mirror of this repo.

A worktree without moosefs access cannot rebuild figures or run figure scripts.
Pre-rendered upstream artifacts that figures vendor are committed in dedicated
folders: input data (BEDs, GFF3, `*.tsv.gz`, PHR gene lists) in `data/`;
genome-wide and inter-chromosomal coverage plots (`p_genome_wide_*`,
`p_combined_*`, `p_num_chromosomes*`) in `inter-chr-plots/`; per-chromosome
identity heatmaps (`identity_heatmap_chr*.pdf`) in `identity_heatmaps/`.

## Layout

- `end-to-end-report/report/01_pipeline.md` … `14_pedigree_recombination.md` —
  14-section analysis report; `end-to-end-report/README.md` is the TOC.
- `paper_prep/figures/{fig1..fig4, ed1..ed5, ed8, nj_tree_arms}/` — one directory
  per paper figure. `ed6` and `ed7` slots are intentionally empty; do not
  create stubs. Each `fig*/`+`ed*/` has a single `figure_<id>.R` (`fig2/`
  uses `figure_fig2.py`), `sources.tsv` (input/output manifest), `caption.md`,
  and rendered PDF+PNG. `nj_tree_arms/` is an exception: script is
  `nj_tree.R`, outputs `nj_tree_annotated.{pdf,png}`, no `sources.tsv` /
  `caption.md`.
- `paper_prep/synthesis/` — paper-level artifacts (`ABSTRACT_BoG.md`,
  `ABSTRACT_nature.md`, `REFERENCES_v3.bib`, `CROSSWALK.md`).
- `paper_prep/submission/` — the active LaTeX manuscript (Springer Nature
  `jnl.cls`): `paper.tex` (single file, no `\input`), `bibliography.bib`,
  bundled `jnl.cls`+`mathphys.bst`, `fig/{MainFigures,ExtendedDataFigures}/`.
  See "Paper submission (LaTeX)" below.
- `paper_prep/lit_review/topic_NN_*.{md,bib}` — 14 per-topic literature
  reviews; `SYNTHESIS_v2.typ` renders the unified review.
- `paper_prep/surveys/SURVEY_NN_*.md` — section-by-section evidence
  inventories. Figure `sources.tsv` lines cite them as `SURVEY_NN §X.Y`.
- `paper_prep/_brainstorming/` — scratch (analysis experiments, validation
  reports); authoritative content is upstream in `paper_prep/figures/` or
  `synthesis/`.
- `slides/v2-review-zoom/` — active deck (BoG 2026). Other `slides/`
  subdirs are earlier rounds (`v2`, `v2-zoom`, `v2-review`) or background
  decks (`chm13-phr-ucsc-browser/`, `meiosis-prophase-background/`).
- `scripts/plot-impg-coverage.R` — standalone R script (run from repo root)
  that reads `data/` inputs (`chm13-annotations.bed`, `hprc25272.CHM13.*.tsv.gz`)
  and emits the coverage plots into `inter-chr-plots/` and the per-chromosome
  identity heatmaps into `identity_heatmaps/`.
- `scripts/{ci,cladistics,hic,mouse,pedigree,popgen}/` — standalone
  reviewer-response statistics (bootstrap CIs, FST jackknife, Mantel tests,
  Monte-Carlo nulls, crossover rates). Filenames carry reviewer-comment tags
  (`*_d_m12`, `*_d_peerq3`, `*_f34`); paired `.py`/`.R` reimplementations and
  `.tsv`/`.json` results sit beside the scripts.
- `scripts/merge_archive_into_graph.py` — wg utility (not analysis code).

## Slide decks (Typst)

Active deck: `slides/v2-review-zoom/_typst/zoom_review_deck.typ` ->
`slides/v2-review-zoom/BoG_2026_review_zoom_v9.pdf`.

Build from `slides/v2-review-zoom/_typst/`:

```
typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v9.pdf
typst compile --root .. --ppi 144 zoom_review_deck.typ page-{0p}.png
```

`typst` is at `~/.cargo/bin/typst`. Each round (`v2` … `v9`) is a frozen
artifact: own `REVISION_NOTES_V<N>.md` (v2 is the unsuffixed
`REVISION_NOTES.md`; v3+ follow the `_V<N>` convention), plus per-slide
experiments under `_revision_assets/v<N>/<slide-id>/` (`README.md`,
`SLIDE_PATCH.md` handoff, `make_*.R` source, candidate PDF/PNG). Bump to
`v<N+1>` for a new revision round; do not edit a frozen `v<N>` in place.

## Paper submission (LaTeX)

Build the manuscript from `paper_prep/submission/` (self-contained; needs only
`pdflatex` + `bibtex`, no moosefs):

```
cd paper_prep/submission
make            # or: bash compile.sh
```

`compile.sh` runs `pdflatex -> bibtex (main) -> bibtex Meth -> bibtex Supp ->
pdflatex -> pdflatex` and must exit 0; it produces `paper.pdf`. Three separate
bibliographies via `multibib`/`newcites` (main, `Meth`, `Supp`) — each needs
its own `bibtex` pass, hence the `.aux`/`.bbl` for `Meth` and `Supp`. `Supp`
may be empty (compile.sh tolerates that bibtex failure). Keep `paper.tex` a
single file — the Springer Nature template forbids `\input{...}`. `make clean`
removes all build artifacts. `BUILD_LOG.md` records the last validated compile.
If a build dies on `Unicode character ... not set up for use with LaTeX`, the
literal Unicode is almost always in a `bibliography.bib` `note`/title field
(the bst prints notes); replace it with a LaTeX macro (`$\geq$`, `$\times$`,
`$\rightarrow$`). `set -e` in `compile.sh` then aborts before bibtex/passes
2-3, so cross-refs come back undefined until the real error is fixed; verify
the FINAL pass via `paper.log` (`grep -c 'undefined' paper.log`), not the
multi-pass `compile.sh` stdout.

The manuscript narrative follows the BoG-2026 talk: deck `paper_prep/paper
figures for Concerted evolution ... .pptx` + transcript
`paper_prep/Session7-PopulationGenomics.en.srt` are the source of truth for the
story (extract slide images with python-pptx; recurse `GROUP` shapes or you
miss grouped figures, e.g. slide 10 and the mouse panel on slide 20). Five main
figures map to deck slides, each panel stacked or side-by-side per author
direction: Fig1 genome-wide homology + PHR lengths (slides 3,4); Fig2 PGGB
hairball on top, tree-ordered + community-ordered Jaccard heatmaps side-by-side
below (6,10,12); Fig3 browser views, each pair stacked — C1 4q/10q DUX4
(14), C2 10p/18p TUBB8B (15), C11 5q/6q OR4F (16); Fig4 sequence-vs-3D scatter +
Pore-C community heatmap (21,22); Fig5 WashU pedigree untangle (24/26).
Extended Data is a SINGLE figure: ED1 = mouse meiotic Hi-C (slide 20),
`fig/ExtendedDataFigures/ED_Fig1_mouse_zygotene.png`. Only `fig:ed1` (mouse) and
`fig:fig1..fig5` exist — do not re-add `\ref{fig:ed2..}`.
Writing follows Guarracino et al. 2023 (the Nature acrocentric-recombination
paper, `~/Downloads/papers/RecombAcroChro_*.pdf`): short bold results-section
headers (`\subsection*{...}` inside the body), a measured "we find / we
observe" cadence, single-paragraph abstract. To keep the draft consistent with
the talk, reviewer-era analyses without a figure (within-community
heterogeneity, popgen/FST, the 14-test 3D forest + controls, CEPH1463, RPE-1,
gene enrichment) were CUT from the body AND the corresponding Methods
subsections; recover from git history if a reviewer asks.
ASSET STATUS: `fig/MainFigures/Fig*.png` (and the ED1 mouse png) are
slide-derived placeholders, some with speaker-note annotations baked in;
publication-quality vector versions must be regenerated from
`paper_prep/figures/`. The original `fig/ExtendedDataFigures/ED_Fig1-7.pdf` are
no longer referenced by `paper.tex`.

## Figures (R / Python)

Each `paper_prep/figures/<id>/` is self-contained: one script reads the
`sources.tsv` inputs and writes `figure_<id>.{pdf,png}` next to it. R scripts
run under guix:

```
guix shell -m /tmp/manifest_fig<id>.scm -- Rscript paper_prep/figures/<id>/figure_<id>.R
```

`guix` lives on the moosefs deploy environment, not on every dev machine. The
exact manifest path is in each script's header comment.

`sources.tsv` schemas are not standardized across figures (columns differ:
`panel/status/role/path/note`, `panel/content/status/source_path/survey_section`,
or `panel/source_path/role/headline_metric`). Common statuses: `READY` (already
on disk, vendored), `GENERATE` (script will produce it). `/moosefs/...` paths
are external inputs.

## Naming and units

- Arms: `chr1_p`, …, `chr22_q`, `chrX_p`, `chrX_q`, `chrY_p`, `chrY_q`
  (48 total, 41 with signal). Some upstream tables use `_parm`/`_qarm`; figure
  scripts normalise before joining.
- Communities: arm-level Leiden k=15 (cross-checked by UPGMA k=14);
  sequence-level Leiden k=50.
- PHR (pseudo-homolog region; formerly "pseudohomologous region"): telomere-anchored 500 kb flank from a contig
  >= 1 Mb that aligned at >= 95% identity to another arm. 15,668 PHRs across
  18,827 flanks.
- CHM13-coordinate PHR BEDs in `data/`: `data/chm13.phrs.bed`,
  `data/chm13.phrs.no_acro.bed`, `data/CHM13-HG002.sub-telo-phrs.bed` (all
  repo-root input data now lives under `data/`).

## Commit convention

Commits originate from workgraph worker agents. Messages look like
`feat: <task-slug> (agent-NNN)`. The `(agent-NNN)` suffix is workgraph
provenance for the graph; keep it.

This guide is written to both `CLAUDE.md` and `AGENTS.md` and kept in
lock-step. The two files exist because Claude Code and Codex CLI look for
different filenames, but they should never drift in content. Any divergence is
a bug. Update both together.
