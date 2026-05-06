# figure_manifest.md — BoG v2 deck

Every figure / image asset referenced across `slides/v2/slide_01_*.md` … `slide_15_*.md`. Synthesised by `bog-v2-slides` (agent-809) on 2026-05-06 from the 15 fanout slide files.

Path conventions:
- relative paths (e.g. `slides/v2/...`, `paper_prep/figures/...`) are anchored to the repo root `/moosefs/erikg/phrs/` (this branch is `wg/agent-809/bog-v2-slides`)
- absolute `/moosefs/guarracino/HPRCv2/...` paths live in Andrea's pipeline tree (read-only from this worktree, listing was confirmed)
- `Exists?` was checked on 2026-05-06 with `ls -e`/`Bash` from the agent worktree
- `Status` enumerates whether the asset is finished, ready-to-render, missing, or optional polish

## 1. Primary figures (one per slide)

| Slide | Source path | Exists? | Status | Notes |
|---|---|:-:|---|---|
| 01 title | (none — title slide, task spec) | n/a | none required | Slide spec explicitly says "no figure". |
| 02 implicit interval tree | `slides/20260204_Subtelomics_overview_EG.pdf` page 2 | **N** | **MISSING — v1 deck PDF not in worktree** | The v1 deck file is referenced by slides 02, 03, and 09 but is **not** present in `slides/`. Decision needed: re-import the v1 PDF, render a substitute interval-tree figure, or use the v1 page 2 image from another tree. |
| 03 IMPG workflow | `slides/20260204_Subtelomics_overview_EG.pdf` page 3 + new ER inset (`slides/v2/slide_03_er_callout.{pdf,png}`) | **N** for v1 page; ER inset **not yet rendered** | **MISSING + READY-TO-RENDER** | v1 deck PDF missing (same root cause as slide 02). The ER-callout R script (n=18,827, p\*=log(n)/n≈5.21e-4, 230× threshold) is embedded in `slide_03_impg_workflow.md`; runs in seconds locally but no PDF/PNG yet exists in `slides/v2/`. |
| 04 genome-wide identity | `paper_prep/figures/fig1/figure_fig1.pdf` panel (a) | Y | ready | Polished publication-quality version of v1 slide 4; chr18q inset already embedded. Optional PAR2 overlay R script provided in `slide_04_genome_wide_identity.md`; not yet rendered. |
| 04 alt | `identity_heatmap_chr18.zoom_last1mb.pdf` | unverified | optional | Single-chromosome zoom alternative; not checked into the worktree (likely lives in `/moosefs/guarracino/...`). Verify before swapping. |
| 05 interchrom similarities | `p_num_chromosomes_wide.pdf` (worktree root) | Y | ready | The v1-slide-5 figure preserved verbatim. Cross-ref `paper_prep/figures/fig1/figure_fig1.pdf` panel **(b)** as manuscript-quality companion (do **not** swap). Optional PHR-scale annotation R snippet inline. |
| 06 length distributions | `/moosefs/guarracino/HPRCv2/PHR_III/plots/all-vs-all.1Mb.p95.id95.len_length_dist_by_chr_arm.pdf` | Y | ready (off-tree) | Faceted histogram per arm, p blue / q orange, pink for missing introvert arms. v1 figure preserved verbatim. Read-only path. |
| 06 callout asset | `slides/v2/slide_06_clade_callouts.{pdf,png}` (R script in slide file) | **N** | READY-TO-RENDER | Pure ggplot2; runs in seconds; no SBATCH; renders the C1/C2/C7/C14/C15 callout legend. |
| 07 heatmap | `paper_prep/figures/fig1/figure_fig1.pdf` panel (c) | Y | ready | 41×41 arm-level Jaccard distance heatmap with UPGMA dendrogram + 15 Leiden community boxes. |
| 07 NJ tree | `paper_prep/figures/nj_tree_arms/nj_tree_annotated.pdf` (PNG companion `.png`) | Y | ready | Produced by upstream `nj-tree-from` task (commit `602a9d3` per slide 07 notes). All six abstract clades labelled in bold. |
| 07 composite | `slides/v2/slide_07_heatmap_nj_combined.png` (R script in slide file) | **N** | READY-TO-RENDER (optional) | Light `magick` stitch of the two PDFs into one landscape composite. Optional polish. |
| 08 PCA-by-chrom | `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-chromosome.pdf` (PNG companion) | Y | ready (off-tree) | Existing artifact. Note: this is **MDS / PCoA**, not strict PCA — see `coherence_check.md` for the labelling decision. |
| 08 PCA-by-superpop | `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-superpopulation.pdf` (PNG companion) | Y | ready (off-tree) | Existing artifact, paired with the by-chrom panel. |
| 08 composite | `slides/v2/slide_08_pca_combined.png` (R script in slide file) | **N** | READY-TO-RENDER (optional) | Light `magick` stitch into a single side-by-side. Optional polish. |
| 09 PCA-by-community | `slides/20260204_Subtelomics_overview_EG.pdf` page 10 | **N** | **MISSING** | The v1-deck PCA scatter is the recommended primary; the v1 deck PDF is not in this worktree (same gap as slides 02 / 03). Alt: render a fresh community-coloured version of the MDS scatter from `/moosefs/guarracino/...arm-leiden-k15.communities.pdf` (which **does** exist on disk and appears suitable — see "Available alternates" below). |
| 09 callout asset | `slides/v2/slide_09_clade_legend.{pdf,png}` (R script in slide file) | **N** | READY-TO-RENDER | Pure ggplot2; per slide-09 logs the R script was verified renders cleanly in this worktree's R install — but the **rendered files are not committed**. Re-run before deck assembly. |
| 10 Hi-C bulk | `paper_prep/figures/fig3/figure_fig3.pdf` panel (a) (PNG companion `.png`) | Y | ready | HG002 Pore-C contact matrix, 50 kb, 77 arm-haplotypes ordered by sequence community. |
| 10 Mantel exclusions | `paper_prep/figures/ed5/figure_ed5.pdf` panel (b) (PNG companion `.png`) | Y | ready | Mantel ρ before vs after acrocentric+sex exclusion. |
| 10 composite | `slides/v2/slide_10_bulk_mantel_composite.pdf` (R script in slide file) | **N** | READY-TO-RENDER (optional) | `cowplot` composite of the two panels with the headline-numbers banner. |
| 11 single-cell 3D | `paper_prep/figures/fig3/figure_fig3.pdf` panel (c) | Y | ready | Per-cell W/B for both GM12878 (Tan 2018) and sperm (Xu 2025), with S_all negative control. |
| 11 alt: GM12878 | `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_16cells_phr_mantel_scatter.pdf` | Y | optional (off-tree) | Stand-alone GM12878 Mantel scatter. |
| 11 alt: sperm | `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_mantel_scatter.pdf` | Y | optional (off-tree) | Stand-alone sperm Mantel scatter. |
| 11 composite | `slides/v2/slide_11_dipc_sperm_pair.pdf` (R script in slide file) | **N** | READY-TO-RENDER (optional) | `cowplot`/`magick` two-panel composite of the alts above. |
| 12 mouse meiotic | `paper_prep/figures/fig4/figure_fig4.pdf` panel (d) | Y | ready | Mouse zygotene per-PHR-pair scatter (B6+CAST T2T, Zuo 2021). |
| 12 alt | `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp/zuo2021_zygotene_phr_pair_scatter.pdf` | unverified | optional (off-tree) | Standalone zygotene scatter, larger. Existence not yet checked. |
| 12 trajectory inset | `slides/v2/slide_12_stage_trajectory.{pdf,png}` (R script in slide file) | **N** | READY-TO-RENDER | 4 hard-coded Mantel ρ values (lepto/zygo/pachy/diplo); ggplot in seconds. **Load-bearing inset** per slide-12 notes — without it, the "zygotene peak" claim is a single number on faith. |
| 13 pedigree lead | `end-to-end-report/pedigree-plots/washu/PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf` | Y | ready | The *odgi untangle* ribbon for PAN027's maternal haplotype 1 against mother PAN010. The lead figure for slide 13. |
| 13 triptych part 2 | `end-to-end-report/pedigree-plots/washu/PAN027.paternal_hap2_from_PAN011_father.untangle.pdf` | Y | ready (optional) | Paternal-haplotype companion. |
| 13 triptych part 3 | `end-to-end-report/pedigree-plots/washu/PAN028.maternal_hap1_from_PAN027_mother.untangle.pdf` | Y | ready (optional) | Three-generation transmission. |
| 13 triptych composite | `slides/v2/slide_13_washu_triptych.pdf` (R script in slide file) | **N** | READY-TO-RENDER (optional) | `cowplot` three-panel landscape; drop the optional asset entirely if 90 s budget gets cut. |
| 14 gene biology | `slides/v2/slide_14_gene_biology.pdf` (and `.png`) | Y | **rendered & committed** | Three-panel figure (DUX4 / OR4F / TAR1) generated by `slides/v2/slide_14_gene_biology.R` from existing cluster TSVs. Already on disk in this worktree. |
| 14 alt | `paper_prep/figures/ed4/figure_ed4.pdf` panels (c)+(d) | Y | reference only | Manuscript ED4 — do **not** substitute for the talk version per slide-14 notes (ED4 has 4 panels including a GO:BP and a copy-weighted panel that are off-topic for the talk). |
| 15 ed8(a) feedback loop | `paper_prep/figures/ed8/figure_ed8.pdf` (PNG `.png`) | Y | ready | Closer schematic — sequence sharing → 3D proximity → ectopic exchange → propagation. Built by Andrea as the discussion synthesis figure. |
| 15 callout overlay | `slides/v2/slide_15_loop_with_callouts.pdf` (R script in slide file) | **N** | READY-TO-RENDER (optional) | Slide-number tag overlay on ed8(a). Optional polish; bullets already do the slide-tagging in words. |

## 2. Source data tables referenced by slide R scripts

(For lineage / methods reproducibility — the talk does not show these, but the slide R scripts read them.)

| Path | Used by | Exists? |
|---|---|:-:|
| `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv` | slide 06 (PHR length distribution data) | Y (off-tree) |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/d4z4_dux4l_by_community.tsv` | slide 14 panel a (DUX4) | Y (off-tree) |
| `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv` | slide 14 panel b (OR4F) | Y (off-tree) |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_by_arm.tsv` | slide 14 panel c (TAR1) | Y (off-tree) |
| `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/all_pedigrees_patches.tsv` | slide 13 numbers (5,984 HQ patches) | unverified (off-tree) |
| `community_based/50000bp/hg002_porec_global_test.tsv` | slide 10 HG002 Pore-C B/W (3.9e-85) | unverified |

## 3. Available alternates (verified on 2026-05-06)

These exist in `/moosefs/guarracino/HPRCv2/PHR_III/similarity/` and are useful fallbacks if the v1 deck PDF stays missing:

- `hprcv2.1Mb.subtelo.mds.arm-leiden-k15.communities.pdf` — community-colored arm-level MDS (good slide-09 substitute)
- `hprcv2.1Mb.subtelo.mds.arm-leiden-k15.superpop-hulls.pdf` — same with superpop hull overlay
- `hprcv2.1Mb.subtelo.mds.arm-upgma-k14.communities.pdf` — UPGMA-coloured equivalent
- `hprcv2.1Mb.subtelo.mds.color-by-arm.pdf` — by-arm coloring
- `hprcv2.1Mb.subtelo.mds.seq-upgma-k150.communities.pdf` — sequence-level UPGMA

## 4. Generation-status summary

- **15 ready (committed/published assets, including 1 in this worktree at slide_14_gene_biology.{pdf,png})**
- **8 ready-to-render** (R scripts inline in slide files; need a single `Rscript` pass before deck assembly): slide 03 ER inset, slide 06 callout, slide 07 stitch, slide 08 stitch, slide 09 legend, slide 10 composite, slide 11 composite, slide 12 trajectory, slide 13 triptych, slide 15 overlay (8 of these are flagged optional / polish; 2 are load-bearing — slide 09 legend and slide 12 trajectory inset).
- **3 missing** — all are pages of the v1 deck PDF (`slides/20260204_Subtelomics_overview_EG.pdf`, used by slides 02, 03, 09). The v1 deck PDF is **not in this worktree**. See `coherence_check.md` §1 for the decision flag.

## 5. Critical-path actions (lead-author / synthesizer call)

1. **Resolve the missing v1 deck PDF** — slides 02, 03, 09 all reference it. Either (a) restore the file to `slides/`, (b) re-render the page-2 / page-3 / page-10 visuals from current data, or (c) substitute with an alt (slide-09 has a clean alt at `arm-leiden-k15.communities.pdf`; slides 02 and 03 are harder — the implicit-interval-tree diagram and the wave→dotplot→forest visual would each need to be rebuilt or extracted).
2. **Render the load-bearing R scripts** — slide_09_clade_legend, slide_12_stage_trajectory. Without these the named-clades legend and the zygotene-peak trajectory rely on hand-drawn Keynote substitutes.
3. **Decide on optional composites** — eight optional R scripts produce side-by-side composites. They are pure `magick`/`cowplot` wrappers (no SBATCH); a 5-minute `Rscript` pass would render all of them. The synthesizer's call is whether the deck is built around composite PDFs (cleaner exports) or live-arranged in Keynote (more flexible).
