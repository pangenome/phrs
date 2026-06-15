# PPTX slide image provenance

Date: 2026-05-21
Source pptx: `/home/guarracino/Downloads/paper figures for Concerted evolution and unorthodox recombination of human subtelomeres.pptx` (also present untracked as `paper_prep/paper figures for ....pptx`)
Slides audited: 4, 6, 10, 12, 19, 21, 22

## Method

1. Unzipped the pptx and read each `ppt/slides/_rels/slideN.xml.rels` to map slide N to its embedded media file `ppt/media/imageK.png`.
2. Visually inspected each image and the slide text (`<a:t>` runs).
3. Grepped the repo for the verbatim titles, stats, and captions burnt into each image.
4. Confirmed the match by reading each candidate script's `README.md` / `SLIDE_PATCH.md`.

MD5s of the pptx PNGs do not byte-match the candidate PNGs on disk. This is expected: PowerPoint recompresses on insert. Content match is unambiguous (titles, axis labels, sample counts, p-values, statistic values all match exactly).

These are NOT the paper figures in `paper_prep/figures/`. They are deck-specific revision assets for the `slides/v2-review-zoom/` deck (BoG 2026 review zoom deck, revision rounds v3..v9).

## Mapping

### Slide 4 — PHR length distribution by chromosome end

- pptx image: `ppt/media/image14.png`
- Script: `slides/v2-review-zoom/_revision_assets/v9/06a_q_axis_kbp/make_06a_q_axis_kbp.R`
- Output: `phr_length_arm_heatstrip_10kbp.png` and `.pdf` (next to script)
- Inputs:
  - `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv` (override with `PHR_LENGTH_TSV=...`)
- Run:

  ```bash
  Rscript slides/v2-review-zoom/_revision_assets/v9/06a_q_axis_kbp/make_06a_q_axis_kbp.R
  ```

### Slide 6 — PGGB ODGI layout main component (component 8)

- pptx image: `ppt/media/image25.png`
- Script: `slides/v2-review-zoom/_revision_assets/v6/pggb_graph_black/render_pggb_layout_component8_black.R` (wrapper: `render_pggb_graph_black.sh`)
- Output: `pggb_graph_2d_black.png`
- Inputs:
  - ODGI layout TSV: `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.tsv`
  - Component 8 only; 727,156 layout nodes; X/Y swapped for 16:9 fit
- Run:

  ```bash
  ./slides/v2-review-zoom/_revision_assets/v6/pggb_graph_black/render_pggb_graph_black.sh
  ```

### Slide 10 — 07a.1 Tree-ordered arm similarity heatmap

- pptx image: `ppt/media/image7.png`
- Script: `slides/v2-review-zoom/_revision_assets/v5/07a_tree_then_community_heatmap/make_07a_tree_then_community_heatmap.R`
- Output: `07a_tree_ordered_heatmap.png` and `.pdf`
- Inputs (positional args; defaults in script):
  1. arm-level distance matrix `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`
  2. arm-level Leiden k=15 assignment TSV
  3. `paper_prep/figures/fig1/architecture_per_arm.tsv`
  4. CHM13 PHR BED
  5. projected all-vs-all interval TSV
  6. CHM13 chrom.sizes
  7. output directory
- Run:

  ```bash
  Rscript slides/v2-review-zoom/_revision_assets/v5/07a_tree_then_community_heatmap/make_07a_tree_then_community_heatmap.R
  ```

### Slide 12 — 07a.2 Community-ordered arm similarity heatmap

- pptx image: `ppt/media/image4.png`
- Script: same as slide 10 (one run produces both panels)
- Output: `07a_community_ordered_heatmap.png` and `.pdf`
- Inputs: same as slide 10

### Slide 19 — Mouse zygotene scatter + Mantel stage trajectory

- pptx image: `ppt/media/image20.png`
- Script: `slides/v2-review-zoom/_revision_assets/hic_visual_redesign/make_hic_visual_redesign.R` (specifically the block at line 251 that writes the `slide_12_mouse_zygotene_*` PNG)
- Output: `slide_12_mouse_zygotene_trajectory_pairing.png`
- Inputs: mouse Zuo et al. 2021 per-PHR-pair Jaccard + Hi-C zygotene table; per-stage Mantel summary; exact paths inside the R script under `/moosefs/guarracino/HPRCv2/PHR_III/`
- Stats burnt into image: Spearman ρ = 0.715, p = 4.4e-55, n = 344 pairs; per-stage Mantel ρ (lepto 0.681, zygo 0.718, pachy 0.683, diplo 0.577)
- Run:

  ```bash
  Rscript slides/v2-review-zoom/_revision_assets/hic_visual_redesign/make_hic_visual_redesign.R
  ```

### Slide 21 — Human sequence similarity vs 3D contact at 50 kb

- pptx image: `ppt/media/image26.png`
- Script: `slides/v2-review-zoom/_revision_assets/v3/human_3d_dotplot/make_human_3d_dotplot.R`
- Output: `human_arm_pair_dotplot_candidate.png` and `.pdf`
- Inputs:
  - `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_phr_pair_correlation.tsv` (HG002 Pore-C; n=803, ρ=0.485, p=1.6e-48)
  - `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_phr_pair_correlation.tsv` (CHM13 Hi-C; n=688, ρ=0.674, p=2.96e-92)
- Unit: one dot per inter-chromosomal arm pair. Statistic is pointwise Spearman, NOT Mantel.

### Slide 22 — HG002 Pore-C contacts ordered by sequence community (v4 orientation fix)

- pptx image: `ppt/media/image22.png`
- Script: `slides/v2-review-zoom/_revision_assets/v4/10a_xaxis_orientation/make_10a_xaxis_orientation.R`
- Output: `candidate_10a_xaxis_orientation.png` and `.pdf` (1800x1800 square; 10x10 in PDF)
- Inputs:
  - contact matrix: `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_contact_matrix.tsv`
  - community boxes: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`
  - B/W statistic source: `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_global_test.tsv`
- Stats burnt into image: B/W = 0.056, p = 3.9e-85
- Fix vs v3: v3 used `rasterImage(as.raster(t(color_matrix)[, n:1]))` which mirrors the X axis. v4 draws `as.raster(color_matrix)` directly and recomputes box coordinates in the corrected displayed system.

## Summary table

| Slide | pptx image | Script | Output PNG/PDF |
|---|---|---|---|
| 4 | image14.png | `v9/06a_q_axis_kbp/make_06a_q_axis_kbp.R` | `phr_length_arm_heatstrip_10kbp.*` |
| 6 | image25.png | `v6/pggb_graph_black/render_pggb_layout_component8_black.R` | `pggb_graph_2d_black.png` |
| 10 | image7.png | `v5/07a_tree_then_community_heatmap/make_07a_tree_then_community_heatmap.R` | `07a_tree_ordered_heatmap.*` |
| 12 | image4.png | same as slide 10 | `07a_community_ordered_heatmap.*` |
| 19 | image20.png | `hic_visual_redesign/make_hic_visual_redesign.R` | `slide_12_mouse_zygotene_trajectory_pairing.png` |
| 21 | image26.png | `v3/human_3d_dotplot/make_human_3d_dotplot.R` | `human_arm_pair_dotplot_candidate.*` |
| 22 | image22.png | `v4/10a_xaxis_orientation/make_10a_xaxis_orientation.R` | `candidate_10a_xaxis_orientation.*` |

All paths are relative to the repo root `/home/guarracino/git/phrs/`. All `/moosefs/...` inputs require moosefs mount.
