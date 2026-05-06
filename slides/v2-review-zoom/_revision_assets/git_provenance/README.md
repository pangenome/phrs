# Review Zoom Deck Git Provenance Audit

Task: `review-zoom-git-provenance-audit`  
Date: 2026-05-06 UTC  
Scope: `slides/v2-review-zoom` only. No deck source was edited by this audit.

This note reconstructs the current review-zoom deck lineage from git history and
blob identity. It is intended as a routing document for the downstream slide
revision tasks, not as an instruction to keep the current assets.

## Executive Summary

- The current review zoom deck is the agent-951 render, committed as
  `10bee88` on branch `wg/agent-951/render-bog-annotated-zoom` and
  squash-merged to `main` as `4862ec7` (`feat: render-bog-annotated-zoom
  (agent-951)`) at 2026-05-06 19:24 UTC.
- The previous zoom review deck is agent-878, committed as `cb973ab` on
  `wg/agent-878/build-bog-v2-3` and squash-merged as `fd9a250`. Its Typst
  source used many absolute `.wg-worktrees/agent-878` paths.
- The v2 synthesis layer is agent-809, committed as `fa8d410` on
  `wg/agent-809/bog-v2-slides` and squash-merged as `8e44af8`. It created
  `SLIDES_v2_PLAN.md`, `figure_manifest.md`, and `coherence_check.md`.
- The current `slides/v2-review-zoom/_typst/zoom_review_deck.typ` has local
  relative asset paths only. `git blame` attributes the whole file to the
  squash merge `4862ec7`; the branch commit `10bee88` is the direct pre-merge
  implementation.
- Blob-hash comparison shows that most current assets are not newly generated
  in a reproducible way by the review-zoom task. They are exact copies of
  earlier committed PNGs from one of these sources:
  `slides/v2`, `paper_prep/figures`, `slides/v2-zoom/_typst/assets`, or the
  unmerged/side branch `wg/agent-944/fix-bog-review` (`759b673`) review assets.
- Highest provenance risks: the agent-878 absolute paths in
  `slides/v2-zoom/_typst/zoom_deck.typ`, crop PNGs with no committed crop
  recipe, off-tree `/moosefs/guarracino/HPRCv2/...` products, and current assets
  that cite "canonical review assets" even though `slides/v2-review/_typst/assets`
  is not present on `main` at `4862ec7`.

## Evidence Commands Run

Primary commands used:

```bash
git log --oneline --decorate --all -- slides/v2-review-zoom slides/v2-zoom slides/v2 .wg/output/render-bog-annotated-zoom .wg/output/build-bog-v2-3 .wg/output/bog-v2-slides
git show --stat --oneline --decorate 10bee88
git show --stat --oneline --decorate cb973ab
git show --stat --oneline --decorate fa8d410
git show --stat --oneline --decorate 4862ec7
git log --oneline --decorate --all -- slides/v2/slide_03* slides/v2/slide_04* slides/v2/slide_06* slides/v2/slide_07* slides/v2/slide_08* slides/v2/slide_09* slides/v2/slide_10* slides/v2/slide_11* slides/v2/slide_12* slides/v2/slide_14*
git ls-tree -r HEAD slides/v2-review-zoom/_typst/assets
git ls-tree -r HEAD slides/v2-zoom/_typst/assets
git ls-tree -r 759b673 slides/v2-review/_typst/assets
git blame --line-porcelain -- slides/v2-review-zoom/_typst/zoom_review_deck.typ
rg -n "agent-878|agent-813|/moosefs/guarracino|github.com|pangenome/HPRCv2|erdos_renyi" slides/v2-review-zoom slides/v2-zoom slides/v2 .wg/output
```

External hooks checked:

```bash
git ls-remote https://github.com/pangenome/HPRCv2.git HEAD refs/heads/*
git ls-remote https://github.com/ekg/erdos_renyi.git HEAD refs/heads/*
git clone --filter=blob:none --no-checkout https://github.com/pangenome/HPRCv2.git /tmp/HPRCv2-provenance
git clone --filter=blob:none --no-checkout https://github.com/ekg/erdos_renyi.git /tmp/erdos-renyi-provenance
```

External heads observed on 2026-05-06:

- `github.com/pangenome/HPRCv2`: `main` and `HEAD` at
  `d14883c314e683063abe8b461992f12825ccd5ed`; relevant branches include
  `inter-chr-plot` (`ab5c9bed...`) and `num-chromosomes` (`253c675...`).
- `github.com:ekg/erdos_renyi`: `main` and `HEAD` at
  `d9ec48f1945d14f38d0131f56ed288cfc7883e73`; repo contains
  `erdos_renyi_viz.R` and `viz.html`.

## Timeline

| UTC time / order | Commit | Agent/task | What it introduced | Provenance impact |
|---|---:|---|---|---|
| Per-slide fanout | `d18c247` / `31ac5f4` | `bog-v2-slide-03` / agent-759 | `slides/v2/slide_03_impg_workflow.md` | Defines IMPG + ER threshold narrative; rendered ER callout comes later. |
| Per-slide fanout | `6192d63` / `661562f` | `bog-v2-slide-04` / agent-760 | `slides/v2/slide_04_genome_wide_identity.md` | Names `paper_prep/figures/fig1/figure_fig1` and optional identity zooms. |
| Per-slide fanout | `5695f4d` / `72c1256` | `bog-v2-slide-06` / agent-780 | `slides/v2/slide_06_length_distributions.md` | Documents off-tree HPRCv2 length-distribution source and callout plan. |
| Per-slide fanout | `a6ea1da` / `3845658` | `bog-v2-slide-07` / agent-801 | `slides/v2/slide_07_allvsall_heatmap_nj_clades.md` | Points to Fig 1 panel c and NJ tree; NJ tree upstream is `602a9d3`. |
| Per-slide fanout | `f1d2478` / `468582d` | `bog-v2-slide-08` / agent-767 | `slides/v2/slide_08_pca_chromosome_superpop.md` | Records that assets are MDS/PCoA, not strict PCA, from HPRCv2 similarity path. |
| Per-slide fanout | `6312a50` / `1628919` | `bog-v2-slide-09` / agent-768 | `slides/v2/slide_09_pca_communities_clades.md` | Reuses v1 page 10 and warns PCA/MDS terminology can drift. |
| Per-slide fanout | `bca2237` / `0e76f59` | `bog-v2-slide-10` / agent-769 | `slides/v2/slide_10_hic_bulk_mantel_exclusions.md` | Defines Fig 3 panel a + ED5 panel source pair. |
| Per-slide fanout | `1255c69` / `b784e3e` | `bog-v2-slide-11` / agent-770 | `slides/v2/slide_11_single_cell_3d.md` | Defines Fig 3 panel c / single-cell source. |
| Per-slide fanout | `964ce4c` / `0772fbe` | `bog-v2-slide-12` / agent-772 | `slides/v2/slide_12_mouse_meiotic_zygotene_bouquet.md` | Defines Fig 4 panel d and trajectory inset need. |
| Per-slide fanout | `863e756` / `6cf052c` | `bog-v2-slide-14` / agent-774 | `slide_14_gene_biology.R/.md/.pdf/.png` | Only slide-specific commit that also rendered its figure. |
| Synthesis | `fa8d410` / `8e44af8` | `bog-v2-slides` / agent-809 | v2 plan, manifest, coherence check | Identified v1 deck gaps, off-tree sources, PCA/MDS risk. |
| v2 render | `beb5036` / `28c2337` | `build-bog-v2-2` / agent-813 | `slides/v2/BoG_2026.pdf`, rendered v1 page images, R-derived callouts | Provides exact source blobs for current 02/03/09, 03b/06b/09b/12b. |
| zoom render | `cb973ab` / `fd9a250` | `build-bog-v2-3` / agent-878 | `slides/v2-zoom`, 28-page A4 zoom deck and crop PNGs | Creates many focused crop PNGs now reused byte-for-byte by current deck; Typst embeds stale absolute worktree paths. |
| annotated review | `f5f9495` / `8413de2`, then `759b673` | `render-bog-annotated`, `fix-bog-review` / agents 941, 944 | `slides/v2-review` review deck and assets on side branch | `759b673` contains "canonical review assets" that match current review-zoom blobs but are not present under `slides/v2-review/_typst/assets` on current `main`. |
| current zoom | `10bee88` / `4862ec7` | `render-bog-annotated-zoom` / agent-951 | `slides/v2-review-zoom`, 26 pages, local assets | Current deliverable. Localizes assets and removes `.wg-worktrees/agent-878` references from this deck. |

## Current Asset Lineage Table

The "source blob" column uses the strongest evidence found: exact git blob
matches when available, otherwise source-path references in committed scripts or
manifests. "Generated script" is intentionally conservative; crop operations are
marked unknown unless a committed script records the crop.

| Slide | Current asset(s) | Likely source path / commit | Generated script if known | Open questions / risk | Recommendation |
|---|---|---|---|---|---|
| 01 | Text-only page; no asset | `zoom_review_deck.typ` in `10bee88` / `4862ec7` | Typst only | No lineage risk. | Keep as deck source text; no asset work needed. |
| 02 | `assets/s02_interval_tree.png` | Exact blob match to `slides/v2/_typst/v1_page_02-02.png` from `28c2337` and to `slides/v2-review/_typst/assets/s02_interval_tree.png` on `759b673` (`5669db27...`). | Extracted/rasterized by agent-813 from v1 deck page 2; extraction recipe not committed. | Source v1 PDF `slides/20260204_Subtelomics_overview_EG.pdf` is still missing from current repo; current asset is a copied raster. | Use current PNG if acceptable; if revising methods visuals, restore or recreate the source schematic rather than editing the raster. |
| 03a | `assets/s03_impg_workflow.png` | Exact blob match to `slides/v2/_typst/v1_page_03-03.png` from `28c2337` and review asset on `759b673` (`bb1c203e...`). | Extracted/rasterized by agent-813 from v1 deck page 3; recipe not committed. | Same missing-v1-PDF risk as slide 02. | Prefer a regenerated workflow visual if slide 03 is revised; current PNG is a raster fallback. |
| 03b | `assets/s03_er_callout.png` | Exact blob match to `slides/v2/slide_03_er_callout.png` from `28c2337` and review asset on `759b673` (`dcdb4cd1...`). | `slides/v2/_typst/slide_03_er_callout.R` added in `28c2337`; conceptual external hook `github.com:ekg/erdos_renyi` at `d9ec48f`. | Current callout is a small PNG, not the external erdos_renyi visualization requested by downstream task. | Downstream `review-zoom-03b-erdos-renyi` should inspect `https://github.com/ekg/erdos_renyi` and decide whether to replace this local R callout. |
| 04 | `assets/s04_fig1_panel_a.png`; sibling full `assets/s04_fig1.png` | `s04_fig1.png` exact blob match to `paper_prep/figures/fig1/figure_fig1.png` from `7fde8f2` (`c3c3fcfd...`). Panel crop exact blob match to `slides/v2-zoom/_typst/assets/slide_04_fig1_panel_a.png` from `cb973ab` (`c2418de0...`). | Full figure: `paper_prep/figures/fig1/figure_fig1.R`. Panel crop recipe not committed. | Crop is copied from agent-878 zoom deck; crop boundaries are not reproducible from git. External HPRCv2 has karyogram image history that may supersede this slide. | For `review-zoom-04-hprcv2-karyogram`, compare current Fig 1 panel a against HPRCv2 karyogram outputs from `github.com/pangenome/HPRCv2` commits `ab5c9be` and `4d4b5c3`; do not assume this crop is canonical. |
| 05 | `assets/s05_interchrom.png` | Exact blob match to `slides/v2-review/_typst/assets/s05_interchrom.png` on `759b673` (`b436facd...`). Source figure is root `p_num_chromosomes_wide.png` (`37a8a938...`) from `ee4dfdb`, but current asset is a resized/cropped derivative. | Original generator not obvious in this repo; v2 manifest points to v1 slide 5 / `p_num_chromosomes_wide.pdf`. | Review asset branch is not merged into current `slides/v2-review/_typst/assets`; resize recipe unknown. | Treat as manually prepared review asset; if revising, go back to `p_num_chromosomes_wide.pdf/png` or Fig 1 panel b source rather than scaling this copy. |
| 06a | `assets/s06_length_dist_1_top.png`, `_2_top.png`, `_1_bottom.png`, `_2_bottom.png`; sibling `assets/s06_length_dist.png` | Four visible assets are exact blob matches to `slides/v2-zoom/_typst/assets/slide_06_length_dist_*` from `cb973ab`. `s06_length_dist.png` matches review asset branch `759b673` (`4aad69f...`). | Underlying source documented as `/moosefs/guarracino/HPRCv2/PHR_III/plots/all-vs-all.1Mb.p95.id95.len_length_dist_by_chr_arm.pdf`; crop/split recipe not committed. | Current slide is four crop PNGs from agent-878. It is easy to revise one crop and forget the companion crop set. | Downstream slide 06 task should regenerate the split view from the original PDF/data with a committed script or document the crop geometry. |
| 06b | `assets/s06_clade_callouts.png` | Exact blob match to `slides/v2/slide_06_clade_callouts.png` from `28c2337` and review asset on `759b673` (`d9f8fa2d...`). | `slides/v2/_typst/slide_06_clade_callouts.R` added in `28c2337`. | The callout lists five visible outlier clades while v2 coherence check flags the "six clades" wording issue. | If slide 06 is redesigned, update both the length plot and the callout text together; reuse the R script, not the PNG. |
| 07a | `assets/s07_fig1_panel_c.png`; sibling `assets/s04_fig1.png` | Panel crop exact blob match to `slides/v2-zoom/_typst/assets/slide_07_fig1_panel_c.png` from `cb973ab` (`079f0cf0...`). Full Fig 1 source is `paper_prep/figures/fig1/figure_fig1.png` from `7fde8f2`. | Full figure: `paper_prep/figures/fig1/figure_fig1.R`. Crop recipe not committed. | Crop hides context and has unknown crop geometry. | For heatmap redesign, regenerate from Fig 1 inputs or the R script rather than patching the crop. |
| 07b | `assets/s07b_nj_tree.png` | Exact blob match to `paper_prep/figures/nj_tree_arms/nj_tree_annotated.png` from `dc5032d` / `602a9d3` and review asset on `759b673` (`9cfbdbc5...`). | `paper_prep/figures/nj_tree_arms/nj_tree.R`; README present in same folder. | Source is good, but current deck uses only the PNG. | Reuse upstream NJ tree script if rerendering labels/tree layout. |
| 08a | `assets/s08a_mds_chrom.png` | Exact blob match to `slides/v2/_typst/img/slide_08_mds_chrom.png` from `28c2337` and review asset on `759b673` (`1b580689...`). Off-tree source: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-chromosome.png`. | Off-tree HPRCv2 script documented in slide 08 as `scripts/similarity/plot-similarity-subtelo.R`, `cmdscale(...)`; not in this repo. | PCA/MDS terminology drift. Current file name says MDS, some deck labels say PCA/PC. | Keep "MDS/PCoA" language unless a real PCA artifact is generated. Quantification task should cite the off-tree script/path. |
| 08b | `assets/s08b_mds_superpop.png` | Exact blob match to `slides/v2/_typst/img/slide_08_mds_superpop.png` from `28c2337` and review asset on `759b673` (`445c2e85...`). Off-tree source: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-superpopulation.png`. | Same off-tree HPRCv2 MDS script as 08a. | Same PCA/MDS risk; superpopulation dispersion may need quantitative overlay not present here. | Downstream 08b should regenerate or annotate from the HPRCv2 data, not edit this PNG. |
| 09a | `assets/s09_pca_communities.png` | Exact blob match to `slides/v2/_typst/v1_page_10-10.png` from `28c2337` and review asset on `759b673` (`17c2b7f4...`). | Extracted/rasterized by agent-813 from missing v1 deck page 10; recipe not committed. | It is a v1 raster and uses PCA wording despite MDS/PCoA caveat. | Downstream 09 should consider replacing with a current MDS/community render from HPRCv2 similarity outputs. |
| 09b | `assets/s09_clade_legend.png` | Exact blob match to `slides/v2/slide_09_clade_legend.png` from `28c2337` and review asset on `759b673` (`4371b19b...`). | `slides/v2/_typst/slide_09_clade_legend.R` added in `28c2337`. | Legend may diverge from any regenerated 09a community plot. | If 09a is rerendered, rerun/update the legend script in the same change. |
| 09c | `assets/s09b_communities.png` | Exact blob match to `slides/v2-review/_typst/assets/s09b_communities.png` on `759b673` (`e746f53d...`). Likely from HPRCv2 community-colored MDS alternate named in `figure_manifest.md`, but no committed generator found here. | Unknown; likely off-tree HPRCv2 similarity pipeline. | Generating script and exact source path are not obvious from current repo. | Verify against `/moosefs/guarracino/HPRCv2/PHR_III/similarity/*communities*` before using it as canonical. |
| 10a | `assets/s10_fig3_panel_a.png`; sibling full `assets/s10_fig3.png` | Full `s10_fig3.png` exact blob match to `paper_prep/figures/fig3/figure_fig3.png` from `8a52549` (`d0199104...`). Panel crop exact blob match to `slides/v2-zoom/_typst/assets/slide_10_fig3_panel_a.png` from `cb973ab` (`adedf4bc...`). | Full figure: `paper_prep/figures/fig3/figure_fig3.R`. Crop recipe not committed. | Crop may not be square or may lose matrix context; downstream task explicitly asks for Hi-C visual redesign. | Regenerate from Fig 3 inputs/script or make a new square matrix asset with committed code. |
| 10b | `assets/s10b_ed5.png` | Exact blob match to `paper_prep/figures/ed5/figure_ed5.png` from `09f6e50` and review asset on `759b673` (`b299dc2b...`). | `paper_prep/figures/ed5/figure_ed5.R`. | Current slide shows full ED5 PNG; if only panel b is intended, it is not a focused crop. | Use ED5 script/source tables for Mantel explanation; avoid hand-cropping without recording geometry. |
| 11 | `assets/s10_fig3_panel_c.png`; sibling full `assets/s10_fig3.png` | Panel crop exact blob match to `slides/v2-zoom/_typst/assets/slide_11_fig3_panel_c.png` from `cb973ab` (`b1d3ab4b...`). Full Fig 3 source is `paper_prep/figures/fig3/figure_fig3.png` from `8a52549`. | Full figure: `paper_prep/figures/fig3/figure_fig3.R`. Crop recipe not committed. | Slide number uses `s10_` filename because it is Fig 3 panel c; downstream may confuse 10/11 lineage. | Rename only if changing deck source; otherwise document that slide 11 is Fig 3c. Regenerate from Fig 3 script for redesign. |
| 12a | `assets/s12_fig4_panel_d.png`; sibling full `assets/s12_fig4.png` | Full `s12_fig4.png` exact blob match to `paper_prep/figures/fig4/figure_fig4.png` from `4a1ee16` (`799f9039...`). Panel crop exact blob match to `slides/v2-zoom/_typst/assets/slide_12_fig4_panel_d.png` from `cb973ab` (`d1b70e2f...`). | Full figure: `paper_prep/figures/fig4/figure_fig4.R`. Crop recipe not committed. | Crop recipe is missing; stage trajectory on 12b is separate and can drift from Fig 4 panel. | Keep Fig 4 script as source of truth. If the mouse slide is rebuilt, render panel d and trajectory together. |
| 12b | `assets/s12_trajectory.png` | Exact blob match to `slides/v2/slide_12_stage_trajectory.png` from `28c2337` and review asset on `759b673` (`98dec76c...`). | `slides/v2/_typst/slide_12_stage_trajectory.R` added in `28c2337`. | Load-bearing inset; values are hard-coded in the R script. | Re-run/update the R script if values or visual framing change. |
| 13a | `assets/s13_pedigree_top.png`; sibling `assets/s13_pedigree.png` | Top crop exact blob match to `slides/v2-zoom/_typst/assets/slide_13a_pedigree_top.png` from `cb973ab` (`79bf16f6...`). Full `s13_pedigree.png` likely rasterized from `end-to-end-report/pedigree-plots/washu/PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf` (`d69117e` lineage), but exact raster recipe not committed. | Unknown for rasterization/crop. | Crop boundaries and PDF-to-PNG settings are not reproducible. | Use the WashU untangle PDF as source of truth; commit a crop/raster script if revising. |
| 13b | `assets/s13_pedigree_bottom.png`; sibling `assets/s13_pedigree.png` | Bottom crop exact blob match to `slides/v2-zoom/_typst/assets/slide_13b_pedigree_bottom.png` from `cb973ab` (`03b24ba3...`). | Unknown for rasterization/crop. | Same as 13a; top/bottom crops must stay paired. | Treat 13a/13b as a manually split view; regenerate both halves together. |
| 14a | `assets/s14_dux4.png`; sibling `assets/s14_gene_biology.png` | `s14_gene_biology.png` exact blob match to `slides/v2/slide_14_gene_biology.png` from `863e756` / `6cf052c` and review asset on `759b673` (`e3c08b6c...`). `s14_dux4.png` exact blob match to `slides/v2-zoom/_typst/assets/slide_14a_dux4.png` from `cb973ab` (`495f504a...`). | Full figure: `slides/v2/slide_14_gene_biology.R`. Panel crop recipe not committed. | Cropped panel may decouple from full R script if revised. | Downstream slide 14 work should start from `slide_14_gene_biology.R` and off-tree tables, not this crop. |
| 14b | `assets/s14_or4f.png`; sibling `assets/s14_gene_biology.png` | Exact blob match to `slides/v2-zoom/_typst/assets/slide_14b_or4f.png` from `cb973ab` (`d06d3348...`). | Full figure: `slides/v2/slide_14_gene_biology.R`; panel crop recipe not committed. | OR4F panel is specifically called out by downstream tasks; crop is not a reproducible analysis product. | Rebuild OR4F from the R script/source CSV or replacement enrichment analysis. |
| 14c | `assets/s14_tar1.png`; sibling `assets/s14_gene_biology.png` | Exact blob match to `slides/v2-zoom/_typst/assets/slide_14c_tar1.png` from `cb973ab` (`c842f98c...`). | Full figure: `slides/v2/slide_14_gene_biology.R`; panel crop recipe not committed. | Same crop risk; TAR1 source is off-tree HPRCv2 enrichment table. | Rebuild from script/source table if content changes. |
| 15 | Text-only current page; unused `assets/s15_ed8.png` exists | `assets/s15_ed8.png` exact blob match to `paper_prep/figures/ed8/figure_ed8.png` from `3f30166` and review asset on `759b673` (`cc171b00...`), but current `zoom_review_deck.typ` does not reference it. | `paper_prep/figures/ed8/figure_ed8.R` for unused asset; Typst text slide for current page. | Asset is stale/unused in current review-zoom deck; can mislead downstream workers into thinking slide 15 has a figure. | Leave unused asset alone unless deck source is changed; document it as unused. |

## Current vs Earlier Zoom Deck

`slides/v2-zoom/_typst/zoom_deck.typ` from agent-878 is useful history but
should not be copied forward directly. It embeds absolute paths such as:

- `/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2/_typst/v1_page_02-02.png`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_04_fig1_panel_a.png`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-878/paper_prep/figures/nj_tree_arms/nj_tree_annotated.png`
- `/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2/slide_14_gene_biology.png`

The current deck fixed that path hygiene problem by copying assets into
`slides/v2-review-zoom/_typst/assets` and using relative paths. However, the
copy step did not make the cropped assets reproducible. Exact blob matches show
these current assets are direct copies of agent-878 crops:

- `s04_fig1_panel_a.png`
- `s06_length_dist_1_top.png`
- `s06_length_dist_2_top.png`
- `s06_length_dist_1_bottom.png`
- `s06_length_dist_2_bottom.png`
- `s07_fig1_panel_c.png`
- `s10_fig3_panel_a.png`
- `s10_fig3_panel_c.png`
- `s12_fig4_panel_d.png`
- `s13_pedigree_top.png`
- `s13_pedigree_bottom.png`
- `s14_dux4.png`
- `s14_or4f.png`
- `s14_tar1.png`

If any downstream task edits one of these, it should either record the crop
geometry in a committed script or replace the crop with a regenerated figure.

## External Provenance Hooks

### Slide 03b: `github.com:ekg/erdos_renyi`

The current slide 03b asset is generated locally by
`slides/v2/_typst/slide_03_er_callout.R`, not by the external repository. The
external hook is still relevant for the downstream ER replacement task:

- Remote: `https://github.com/ekg/erdos_renyi.git`
- Current `HEAD` observed: `d9ec48f1945d14f38d0131f56ed288cfc7883e73`
- Files observed: `erdos_renyi_viz.R`, `viz.html`

Recommendation: treat the current callout as a placeholder. If the goal is to
show Erik's ER visualization, inspect and adapt the external repo rather than
polishing `s03_er_callout.png`.

### Slide 04: `github.com/pangenome/HPRCv2` image history

The current slide 04 is a crop from local manuscript Fig 1, not the external
HPRCv2 karyogram. The HPRCv2 repository contains a relevant, newer image lineage:

- Remote: `https://github.com/pangenome/HPRCv2.git`
- Current `main` / `HEAD` observed: `d14883c314e683063abe8b461992f12825ccd5ed`
- Relevant branch: `inter-chr-plot` at `ab5c9bed552ac33dd63919a3b62ab624ceba2f3f`
- Relevant commits:
  - `ab5c9be` `Add inter-chromosomal mapping plots`
  - `4d4b5c3` `Add new inset karyogram plots and adjust inset positioning in R script`
  - `d14883c` `Add supplementary figure for inter-chromosomal matches in HPRC vs CHM13`
- Relevant files on `main`:
  - `scripts/plot-impg-coverage.inter-chr-map.R`
  - `p_interchrom_karyogram_count_inset.100kb.png`
  - `p_interchrom_karyogram_count_rainbow_inset.100kb.png`

The HPRCv2 script reads
`/home/guarracino/Dropbox/git/HPRCv2/data/hprc25272-wf.CHM13.100kb-xm5-id098-l50000.tsv.gz`
and writes the karyogram inset PNG/PDF files. That path is outside this repo.

Recommendation: the slide 04 downstream task should evaluate those HPRCv2
karyogram files as replacements or companions to the current Fig 1 panel a crop.
Do not assume the current crop is the best slide 04 asset.

## Stale / Provenance Risks

| Risk | Evidence | Impact | Recommended mitigation |
|---|---|---|---|
| `.wg-worktrees/agent-878` paths in earlier zoom deck | `rg` finds many such paths in `slides/v2-zoom/_typst/zoom_deck.typ`; current deck render log says none under `slides/v2-review-zoom`. | Copying old Typst snippets can silently reintroduce dead paths. | Use current review-zoom relative paths or regenerate local assets; never copy agent-878 absolute paths. |
| Review assets are on side branch, not current `main` path | `git ls-tree -r 759b673 slides/v2-review/_typst/assets` has the canonical review assets; current `find slides/v2-review/_typst` shows no `assets/` directory. | The phrase "canonical review asset" in `zoom_review_deck.typ` is historically true but not reproducible from current `main` alone. | Use `slides/v2-review-zoom/_typst/assets` as current local copies; if a canonical review deck is needed, merge/restore the review asset directory explicitly. |
| Crop PNGs without scripts | Agent-878 crop blobs match current assets exactly, but no crop script is committed. | Revisions are hard to reproduce and easy to misalign across paired top/bottom panels. | Add a small committed crop/raster script per revised figure or regenerate from source R/PDF. |
| Off-tree HPRCv2 inputs | Slides 06/08/14 and paper figure scripts cite `/moosefs/guarracino/HPRCv2/...`; local `/moosefs/guarracino/HPRCv2` is not itself a git checkout. | Current PNGs may lag Andrea's pipeline and have weak local history. | Record exact off-tree source path, GitHub repo commit where available, and source-data timestamp/checksum when regenerating. |
| PCA/MDS naming drift | v2 coherence check flags that slide 08 is MDS/PCoA while slide 09/v1 assets say PCA. | Audience-facing methodology inconsistency. | Prefer MDS/PCoA labels unless a real PCA artifact is generated. |
| Unused `s15_ed8.png` in current asset directory | Asset exists and matches `paper_prep/figures/ed8/figure_ed8.png`, but `zoom_review_deck.typ` slide 15 is text-only. | Downstream workers may revise an unused asset. | Check Typst references before editing assets; delete only in a deck-cleanup task if requested. |

## Practical Recommendations for Downstream Tasks

1. Treat `slides/v2-review-zoom/_typst/zoom_review_deck.typ` as the current
   page map, but treat PNG assets as cached artifacts unless their source script
   is listed above.
2. For slides 04, 07a, 10a, 11, and 12a, prefer the corresponding
   `paper_prep/figures/*/figure_*.R` scripts over crop PNG edits.
3. For slides 06a, 13a/b, and 14a/b/c, create or recover crop/raster scripts
   before changing visual geometry.
4. For slides 08 and 09, resolve the MDS/PCoA vs PCA wording first; otherwise
   regenerated visuals and labels will drift.
5. For slide 03b and slide 04 replacement work, use the external GitHub hooks
   above as starting points and record the exact external commit in the next
   revision artifact.
