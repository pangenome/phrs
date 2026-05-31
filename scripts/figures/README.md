# scripts/figures — main-figure generators

Regenerate manuscript figures from files already in this repo (`data/` inputs,
base R, no moosefs). **Generators write straight into the submission**
(`paper_prep/submission/fig/...`), so there is no copy step — just run from the
repo root and rebuild the PDF:

```bash
Rscript scripts/figures/make_<fig>.R          # writes into paper_prep/submission/fig/...
( cd paper_prep/submission && make )          # rebuild paper.pdf
```

Override paths with `DATA_DIR=...` / `OUT_DIR=...` if needed.

## Generators (runnable now)

```bash
Rscript scripts/figures/make_fig4a_human_scatter.R   # -> fig/MainFigures/Fig4a_human_scatter.{png,pdf}
Rscript scripts/figures/make_ed1_mouse_zygotene.R    # -> fig/ExtendedDataFigures/ED_Fig1_mouse_zygotene.{png,pdf}
Rscript scripts/figures/make_fig4_allpoints.R        # -> fig/ExtendedDataFigures/ED_Fig2_allpoints.{png,pdf}
```

| Figure | Inputs (`data/`) | Verified |
|--------|------------------|----------|
| Fig 4a | `hg002_porec_phr_pair_correlation.tsv`, `chm13_phr_pair_correlation.tsv` | HG002 ρ=0.485 (n=803); CHM13 ρ=0.674 (n=688) |
| ED 1   | `zuo2021_zygotene_phr_pair_correlation.tsv` | mouse zygotene ρ=0.715 (n=344); stage peak zygo 0.718 |
| ED 2 (all-points) | `human_HG002_porec_50000bp_seqlevel.tsv`, `human_CHM13_hic_50000bp_seqlevel.tsv` | HG002 Pore-C ρ=0.381 (n=2830); CHM13 Hi-C ρ=0.716 (n=652) |

ED 2 plots **one dot per inter-chromosomal PHR sequence pair** (single-sample: each sample's own PHRs vs its own contact; x = `jaccard`, y = `hic_contact_norm`), the assumption-free companion to Fig 4a, which plots **one dot per arm pair** (x = population mean Jaccard over the 233 samples). It is referenced in `paper.tex` as Extended Data Fig. 2.

## Generator pending (inputs in `data/`, script not written yet)

| Figure | Inputs (`data/`) |
|--------|------------------|
| Fig 1b lengths | `all-vs-all.1Mb.p95.id95.len.tsv` |
| Fig 2a PGGB layout | `hprcv2.1Mb.telo_trimmed.p95.id95...smooth.final.og.lay.tsv` |
| Fig 2b/2c Jaccard heatmaps | `hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`, `hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`, `chm13.chrom.sizes`, `all-vs-all.1Mb.p95.id95.len.tsv`, `chm13.phrs.bed` (+ `paper_prep/figures/fig1/architecture_per_arm.tsv`) |
| Fig 4b Pore-C community heatmap | `hg002_porec_contact_matrix.tsv`, `hg002_porec_hic.arm-leiden.communities.tsv`, `hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`, `hg002_porec_global_test.tsv` |

## Vendored (real in-repo assets — copy, not generated here)

```bash
# Fig 3 — UCSC hs1 browser panels
P=slides/chm13-phr-ucsc-browser/_assets/ucsc/panels
M=paper_prep/submission/fig/MainFigures
cp $P/chr4q_chr4_193304946_193574945.png   $M/Fig3a_C1_chr4q.png
cp $P/chr10q_chr10_134488135_134758134.png $M/Fig3a_C1_chr10q.png
cp $P/chr10p_chr10_1_105002.png            $M/Fig3b_C2_chr10p.png
cp $P/chr18p_chr18_1_397502.png            $M/Fig3b_C2_chr18p.png
cp $P/chr5q_chr5_181752940_182045439.png   $M/Fig3c_C11_chr5q.png
cp $P/chr6q_chr6_171909129_172126628.png   $M/Fig3c_C11_chr6q.png
# (panels themselves are rebuilt by, needs internet + UCSC session:
#  python3 slides/chm13-phr-ucsc-browser/_scripts/render_ucsc_browser_suite.py --force )

# Fig 5 — pedigree untangle
cp end-to-end-report/pedigree-plots/washu/PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf \
   paper_prep/submission/fig/MainFigures/Fig5_pedigree_untangle.pdf
```

## Deferred

- Fig 1a (slide-3 ideogram): needs `plot-impg-coverage.inter-chr-map.R` + its
  `hprc25272-wf.CHM13.100kb-xm5-id098-l50000.tsv.gz`, neither in the repo.
