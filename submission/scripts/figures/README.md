# submission/scripts/figures — main-figure generators

Reproduces every manuscript figure from inputs already in the repo (`data/`,
base R, no moosefs). Run from the repo root unless noted, then rebuild:
`( cd submission && make )`.

## What is computed vs copied

**9 figures are COMPUTED from `data/`** (a script plots them from numbers).
**7 image files are VENDORED — copied, NOT computed here** (Fig3's 6 browser
panels + Fig5): they are screenshots / pre-rendered PDFs whose generators need
the live UCSC browser or off-repo moosefs data, so only the rendered images
live in the repo.

| Figure | How | Generator / source |
|---|---|---|
| Fig1a | computed | HPRCv2 repo `plot-impg-coverage.inter-chr-map.R` (external clone) |
| Fig1b | computed | `make_fig1b_lengths.R` (here) — per-end PHR length heatstrip (terminal 500 kbp, 10 kbp columns, q arms flipped); input `data/all-vs-all.1Mb.p95.id95.len.tsv` |
| Fig2a | computed | `make_fig2a_pggb_layout.R` (here) — odgi layout main component, nodes colored by arm-level Leiden community (same palette as Fig2b); inputs `data/fig2a_pggb_layout.og.lay.tsv.gz` + `data/fig2a_node_community.tsv.gz` |
| Fig2b+2c | computed | `make_fig2bc_jaccard_heatmaps.R` (here) — one combined `Fig2bc_jaccard` file, both panels + one shared legend |
| Fig4a | computed | `make_fig4a_human_scatter.R` (here) — HG002 Pore-C all-points (single panel); full 2-panel HG002+CHM13 kept in `_backup/make_fig4a_human_scatter_full.R` |
| Fig4b | computed | `make_fig4b_porec_community.R` (here) — Pore-C matrix + community colour bands/labels |
| Fig4c | computed | `make_fig4c_mouse_zygotene.R` (here) — mouse all-points 20 kb: zygotene scatter + per-stage per-pair Spearman trajectory (peaks at bouquet); was ED1 |
| ED1 | computed | `make_ed1_human_contacts.R` (here) — 3-panel all-points replicate of Fig4a across genomes/assays (CHM13 Hi-C, HG002 Hi-C, HG002 CiFi), 50 kbp |
| **Fig3a/b/c (6 PNG)** | **VENDORED (not computed)** | UCSC hs1 browser screenshots in `slides/chm13-phr-ucsc-browser/_assets/ucsc/panels/` |
| **Fig5 (PDF)** | **VENDORED (not computed)** | `end-to-end-report/pedigree-plots/washu/...untangle.pdf` |

`paper.tex` includes vectors (`.pdf`) for every computed panel except Fig2a
(`.png`); Fig2b/2c are a single combined `Fig2bc_jaccard.pdf`; the vendored
Fig3 panels are `.png`.

## Computed — scripts here (location-aware, run from anywhere)

```bash
Rscript submission/scripts/figures/make_fig1b_lengths.R           # Fig1b  per-end PHR length heatstrip (terminal 500 kbp, 10 kbp columns, q arms flipped)
Rscript submission/scripts/figures/make_fig2a_pggb_layout.R        # Fig2a  odgi layout main component, nodes colored by arm-level Leiden community
Rscript submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R  # Fig2bc combined tree- + community-ordered Jaccard heatmaps
Rscript submission/scripts/figures/make_fig4a_human_scatter.R      # Fig4a  HG002 Pore-C all-points rho=0.381 n=2830 (single panel)
Rscript submission/scripts/figures/make_ed1_human_contacts.R     # ED1    3-panel replicate: CHM13 Hi-C 0.716/652, HG002 Hi-C 0.662/2544, HG002 CiFi 0.191/2757 (50 kbp)
Rscript submission/scripts/figures/make_fig4b_porec_community.R    # Fig4b  Pore-C contact matrix ordered by sequence community
Rscript submission/scripts/figures/make_fig4c_mouse_zygotene.R     # Fig4c  mouse 20kb per-pair: zygo rho=0.614 n=1135; trajectory lepto/zygo/pachy/diplo 0.419/0.614/0.576/0.245 (needs the 4 mouse_*_phr_20000bp_seqlevel.tsv)
```

These resolve `data/` and `submission/fig/` from their own path (override with
`DATA_DIR=` / `OUT_DIR=`). Fig4a (HG002 Pore-C) and ED1 (CHM13 Hi-C) are each one
dot per inter-chromosomal PHR sequence pair, single-sample (a sample's own PHRs
vs its own length-normalised contact). The full 2-panel HG002+CHM13 version of
the old Fig4a is preserved in `_backup/make_fig4a_human_scatter_full.R`.

## Computed — deck generators (migrated)

Fig1b and Fig2a were previously rendered by frozen slide generators under
`slides/v2-review-zoom/_revision_assets/` and copied in. Both are now built
directly by `make_fig1b_lengths.R` and `make_fig2a_pggb_layout.R` above (the
frozen slide assets are left untouched for provenance;
`paper_prep/20260521_pptx-slide-image-provenance.md`).

## Computed — Fig1a (external `pangenome/HPRCv2` repo)

Cloned at `/home/guarracino/Dropbox/git/HPRCv2`; script + inputs
(`data/hprc25272-wf.CHM13.100kb-xm5-id098-l50000.tsv.gz`,
`data/chm13-annotations.bed`) are all there.

```bash
( cd /home/guarracino/Dropbox/git/HPRCv2 && Rscript scripts/plot-impg-coverage.inter-chr-map.R )
cp /home/guarracino/Dropbox/git/HPRCv2/p_interchrom_karyogram_count_rainbow_inset.100kb.pdf \
   submission/fig/MainFigures/Fig1a_genomewide.pdf
```

## VENDORED — copied, not computed

```bash
# Fig3 — UCSC hs1 browser screenshots (rebuild needs internet + a UCSC session:
#   python3 slides/chm13-phr-ucsc-browser/_scripts/render_ucsc_browser_suite.py --force)
P=slides/chm13-phr-ucsc-browser/_assets/ucsc/panels; M=submission/fig/MainFigures
cp $P/chr4q_chr4_193304946_193574945.png   $M/Fig3a_C1_chr4q.png
cp $P/chr10q_chr10_134488135_134758134.png $M/Fig3a_C1_chr10q.png
cp $P/chr10p_chr10_1_105002.png            $M/Fig3b_C2_chr10p.png
cp $P/chr18p_chr18_1_397502.png            $M/Fig3b_C2_chr18p.png
cp $P/chr5q_chr5_181752940_182045439.png   $M/Fig3c_C11_chr5q.png
cp $P/chr6q_chr6_171909129_172126628.png   $M/Fig3c_C11_chr6q.png

# Fig5 — pedigree untangle (generator is off-repo on moosefs; only the PDF is here)
cp end-to-end-report/pedigree-plots/washu/PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf \
   $M/Fig5_pedigree_untangle.pdf
```
