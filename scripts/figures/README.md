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

## Generators in the deck revision assets (run against `data/`, then copy)

These four panels are NOT pending — their generators already exist as the
BoG-2026 deck revision-asset scripts under `slides/v2-review-zoom/_revision_assets/`
(the submission `.png` placeholders were extracted from those very slides; see
`paper_prep/20260521_pptx-slide-image-provenance.md`). The scripts default to
`/moosefs` paths but every input is vendored in `data/`, so override the paths
and copy the output into the submission. Run from the repo root.

**Fig 1b** (slide 4, PHR length heatstrip) — `PHR_LENGTH_TSV` env override:

```bash
S=slides/v2-review-zoom/_revision_assets/v9/06a_q_axis_kbp
PHR_LENGTH_TSV=data/all-vs-all.1Mb.p95.id95.len.tsv Rscript $S/make_06a_q_axis_kbp.R
cp $S/phr_length_arm_heatstrip_10kbp.pdf paper_prep/submission/fig/MainFigures/Fig1b_lengths.pdf
```

**Fig 2a** (slide 6, PGGB ODGI layout component 8) — `--layout-tsv` flag (PNG only):

```bash
S=slides/v2-review-zoom/_revision_assets/v6/pggb_graph_black
Rscript $S/render_pggb_layout_component8_black.R \
  --layout-tsv data/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.tsv \
  --out $S/pggb_graph_2d_black.png
cp $S/pggb_graph_2d_black.png paper_prep/submission/fig/MainFigures/Fig2a_pggb_layout.png
```

**Fig 2b + 2c** (slides 10 & 12, tree- and community-ordered Jaccard) — one run
produces both; six positional path args + out dir:

```bash
S=slides/v2-review-zoom/_revision_assets/v5/07a_tree_then_community_heatmap
Rscript $S/make_07a_tree_then_community_heatmap.R \
  data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv \
  data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv \
  paper_prep/figures/fig1/architecture_per_arm.tsv \
  data/chm13.phrs.bed \
  data/all-vs-all.1Mb.p95.id95.len.tsv \
  data/chm13.chrom.sizes \
  $S
cp $S/07a_tree_ordered_heatmap.pdf      paper_prep/submission/fig/MainFigures/Fig2b_tree_jaccard.pdf
cp $S/07a_community_ordered_heatmap.pdf paper_prep/submission/fig/MainFigures/Fig2c_community_jaccard.pdf
```

**Fig 4b** (slide 22, Pore-C contacts by community) — paths are hardcoded at
lines 6-9 with NO override hook, so edit those four to the `data/` copies
(`hg002_porec_contact_matrix.tsv`, `hg002_porec_hic.arm-leiden.communities.tsv`,
`hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`, `hg002_porec_global_test.tsv`)
before running:

```bash
S=slides/v2-review-zoom/_revision_assets/v4/10a_xaxis_orientation
# edit S/make_10a_xaxis_orientation.R lines 6-9 -> data/<basename>, then:
Rscript $S/make_10a_xaxis_orientation.R
cp $S/candidate_10a_xaxis_orientation.pdf paper_prep/submission/fig/MainFigures/Fig4b_porec_community.pdf
```

Switching any panel to its `.pdf` means updating the matching
`\includegraphics{...png}` line in `paper.tex` to `.pdf`.

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

## Fig 1a — external repo generator (not in this repo, but reproducible locally)

Fig 1a (genome-wide karyogram, # other chromosomes per 100 kb window) is
generated by the `pangenome/HPRCv2` repo, cloned locally at
`/home/guarracino/Dropbox/git/HPRCv2`. Both the script and its inputs are
present there:

- script: `scripts/plot-impg-coverage.inter-chr-map.R`
- inputs (absolute paths hardcoded in the script):
  `data/hprc25272-wf.CHM13.100kb-xm5-id098-l50000.tsv.gz`,
  `data/chm13-annotations.bed`
- outputs (written to the HPRCv2 repo root):
  `p_interchrom_karyogram_count_rainbow_inset.100kb.{png,pdf}`

```bash
# regenerate, then copy the vector PDF into the submission
( cd /home/guarracino/Dropbox/git/HPRCv2 && Rscript scripts/plot-impg-coverage.inter-chr-map.R )
cp /home/guarracino/Dropbox/git/HPRCv2/p_interchrom_karyogram_count_rainbow_inset.100kb.pdf \
   paper_prep/submission/fig/MainFigures/Fig1a_genomewide.pdf
# paper.tex currently \includegraphics{Fig1a_genomewide.png}; switch to the .pdf for vector
```

A pre-rendered `p_interchrom_karyogram_count_rainbow_inset.100kb.{png,pdf}`
already sits in the HPRCv2 clone, so the copy step works even without re-running.
