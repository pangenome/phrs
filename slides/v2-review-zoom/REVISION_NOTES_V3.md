# BoG 2026 Review Zoom Deck v3 Revision Notes

Task: `render-bog-annotated-zoom-v3-correction-fanin`

Scope: `slides/v2-review-zoom`

This revision builds `BoG_2026_review_zoom_v3.pdf` from
`_typst/zoom_review_deck.typ`, integrating the corrected v3 fanout artifacts.
The v1 and v2 PDFs remain in place:

- `BoG_2026_review_zoom.pdf`
- `BoG_2026_review_zoom_v2.pdf`
- `BoG_2026_review_zoom_v3.pdf`

## Build

Rendered from `slides/v2-review-zoom/_typst` with Typst 0.13.1:

```bash
typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v3.pdf
typst compile --root .. --ppi 144 zoom_review_deck.typ page-{0p}.png
```

The render log is `slides/v2-review-zoom/_typst/render.log`.

## User Criticism Map

| Complaint / correction request | v3 action | Status |
| --- | --- | --- |
| Slide 06a should be a violin/distribution view and must visibly mark the 500 kb analysis cap. | Replaced the v2 ranked-arm graphic with `v3/06_violin_censor/candidate_06a_named_clade_violin_censor.png`. The slide title and source footer call out the 500 kb cap, and the plot marks the cap as a search/analysis limit rather than measured values beyond 500 kb. | Fixed on page 09 / slide `06a`. |
| Slide 07a heatmap/tree was fuzzy, and tree leaves appeared misaligned with rows. | Replaced the v2 07a PNG with `v3/07a_crisp_aligned/candidate_07a_upgma_crisp_aligned.png`. The fanout `order_validation.tsv` records that tree tip order equals heatmap row and column order. | Fixed on page 12 / slide `07a`. |
| Slide 07c unrooted tree may mislead if acrocentric arms are absent or unreadable. | Replaced the unannotated unrooted tree with `v3/07c_acrocentric_presence/07c_unrooted_acrocentric_status.png`, which keeps 07c as a backup and adds an explicit acrocentric p-arm status table. The audit reports all five acrocentric p arms present and in C7. | Fixed on page 14 / slide `07c`; not promoted to the main similarity slide. |
| Slide 08b should discard the center-based spread metric and use a within-population pairwise variation metric. | Updated the method wording and slide visual to `v3/08b_within_pop_pairwise/within_pop_pairwise_2d_distribution.png`, which summarizes same-superpopulation point-to-point distances in displayed MDS / PCoA coordinates. | Fixed on pages 15 and 17 / slides `08m` and `08b`. |
| Slide 09 should use a square 1:1 MDS / PCoA plot with all communities C1-C15 labeled. | Replaced the v2 community plot with `v3/09_all_communities_1to1/mds_pcoa_all_communities_1to1.png`. The candidate is square, equal-scale, and labels every community C1-C15. | Fixed on page 18 / slide `09`. |
| Slide 10a heatmap axes, transpose, and community boxes needed correction. | Replaced the v2 square matrix with `v3/10a_axis_box_fix/candidate_10a_axis_box_fix.png`. The generator asserts identical row and column order and stores community box coordinates in `sequence_community_boxes.tsv`. | Fixed on page 20 / slide `10a`. |
| Slide 11 should not use unexplained distance-ratio shorthand. | Replaced the v2 candidate with `v3/11_wb_labels/slide11_explicit_distance_labels_candidate.png`, and changed the method slide wording to explicit within-community versus between-community distance language. | Fixed on pages 19 and 22 / slides `10m` and `11`. |
| Slide 12 should keep the strong mouse zygotene slide and consider a human sequence-similarity versus 3D-contact analog if valid. | Kept the mouse zygotene slide unchanged. Added `12b` with `v3/human_3d_dotplot/human_arm_pair_dotplot_candidate.png` because the fanout found a valid conservative human arm-pair analog. | Fixed on pages 23 and 24 / slides `12` and `12b`. |
| Slide 14a / gene block should replace the weak DUX4 context slide with focused systematic genome-browser-style PHR gene annotation panels if valid. | Replaced the old DUX4, OR4F, and TAR1 gene block with five v3 genome-browser panels: DUX4/D4Z4 C1, OR4F C3, OR4F decay C8, TAR1 C2, and acrocentric C7. | Fixed on pages 28-32 / slides `14a`-`14e`. |

## Deferred / Retained Items

No required v3 criticism is deferred. The following choices are intentional
retentions or caveats:

| Item | Decision | Reason |
| --- | --- | --- |
| Slide `07b` rooted readable tree | Retained. | The v3 criticism targeted 07a alignment and 07c acrocentric interpretation; 07b remains a readable backup topology. |
| Slide `10b` Mantel exclusion check | Retained from v2 review redesign. | No v3 correction task replaced 10b, and it still explains exclusion robustness after the corrected 10a matrix. |
| Human PHR sequence-pair control plot | Not integrated into the main deck. | The human fanout recommends the arm-pair analog for the regular Fig. 3 / ED5 data. The sequence-pair plot is a `no_strong` community-free control, so adding it as the main analog would overstate the comparison. |
| Older copied/cropped backup assets outside the v3 criticism list | Retained where still used. | This fan-in task owns the correction deck, not a full re-render of all legacy backup crops. Provenance caveats remain in `REVISION_NOTES.md` and the v3 provenance table below identifies all newly integrated assets. |

## New / Changed Asset Provenance

All paths below are relative to `slides/v2-review-zoom`.

| Slide(s) | Integrated asset | Generator | Source data / audit trail |
| --- | --- | --- | --- |
| `06a` | `_revision_assets/v3/06_violin_censor/candidate_06a_named_clade_violin_censor.png` | `_revision_assets/v3/06_violin_censor/make_06_violin_censor.R` | Reads `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv`; summary in `named_clade_violin_summary.tsv`. |
| `07a` | `_revision_assets/v3/07a_crisp_aligned/candidate_07a_upgma_crisp_aligned.png` | `_revision_assets/v3/07a_crisp_aligned/make_07a_crisp_aligned.R` | Uses the canonical 41 x 41 arm-level Jaccard distance matrix and Leiden k=15 assignments documented in the README; row/tip checks in `order_validation.tsv`; order in `arm_order_upgma.tsv`. |
| `07c` | `_revision_assets/v3/07c_acrocentric_presence/07c_unrooted_acrocentric_status.png` | `_revision_assets/v3/07c_acrocentric_presence/make_07c_acrocentric_presence.R` | Audit tables: `arm_presence_complete.tsv`, `acrocentric_p_status.tsv`, `matrix_audit.tsv`, `acrocentric_p_distance_matrix.tsv`; rebuilt tree in `07c_unrooted_nj_audited.newick`. |
| `08b` | `_revision_assets/v3/08b_within_pop_pairwise/within_pop_pairwise_2d_distribution.png` | `_revision_assets/v3/08b_within_pop_pairwise/make_within_pop_pairwise.R` | Coordinates from `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.full_mds.rds`; labels from HPRC production metadata and sequence assignment TSVs; outputs in `within_pop_pairwise_summary.tsv` and `within_pop_pairwise_distance_sample.tsv`. |
| `09` | `_revision_assets/v3/09_all_communities_1to1/mds_pcoa_all_communities_1to1.png` | `_revision_assets/v3/09_all_communities_1to1/make_all_communities_1to1.R` | Coordinates from `hprcv2.1Mb.subtelo.full_mds.rds`; communities from `hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`; checks in `validation_summary.tsv`; label geometry in `label_positions.tsv`. |
| `10a` | `_revision_assets/v3/10a_axis_box_fix/candidate_10a_axis_box_fix.png` | `_revision_assets/v3/10a_axis_box_fix/make_10a_axis_box_fix.R` | Contact matrix `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/hg002_porec_contact_matrix.tsv`; sequence-community table `hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`; audits in `matrix_order_audit.tsv`, `ordered_arm_haplotypes.tsv`, and `sequence_community_boxes.tsv`. |
| `11` | `_revision_assets/v3/11_wb_labels/slide11_explicit_distance_labels_candidate.png` | `_revision_assets/v3/11_wb_labels/make_slide11_explicit_distance_labels.R` | Uses the same Dip-C and sperm per-cell TSV inputs as the v2 redesign; plotted counts summarized in `slide11_explicit_distance_summary.tsv`. |
| `12b` | `_revision_assets/v3/human_3d_dotplot/human_arm_pair_dotplot_candidate.png` | `_revision_assets/v3/human_3d_dotplot/make_human_3d_dotplot.R` | Human 50 kb arm-pair correlation TSVs under `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/`; statistics and source paths in `human_3d_dotplot_summary.tsv`. |
| `14a`-`14e` | `_revision_assets/v3/gene_browser_panels/panel_01_dux4_d4z4_c1_chr4_chr10.png`; `panel_02_or4f_c3_chr3q.png`; `panel_03_or4f_decay_c8_chr15q.png`; `panel_04_tar1_c2_chr18p.png`; `panel_05_acrocentric_c7_p_arm_group.png` | `_revision_assets/v3/gene_browser_panels/render_gene_browser_panels.R` | Panel inventory and track grammar in `panel_manifest.tsv` and `render_track_schema.tsv`; exact input paths in `input_manifest.tsv`; target loci from `_revision_assets/v3/gene_browser_inventory/target_loci.tsv`. |
| All rendered pages | `BoG_2026_review_zoom_v3.pdf`; `_typst/page-01.png` through `_typst/page-33.png` | `_typst/zoom_review_deck.typ` rendered by Typst 0.13.1 | Deck source plus relative asset paths listed above. |

## Page Map

| Page | Label | Visual focus |
| ---: | --- | --- |
| 01 | 01 | Title text-only focus page |
| 02 | 02m | Method transition: implicit graph over chromosome ends |
| 03 | 02 | Backup interval-tree schematic |
| 04 | 03a | IMPG workflow |
| 05 | 03b | Erdos-Renyi connectivity replacement plot |
| 06 | 04 | HPRCv2 interchrom karyogram |
| 07 | 04b | Backup manuscript Fig 1a genome-wide identity heatmap |
| 08 | 05 | Backup genome-wide count view |
| 09 | 06a | Violin/distribution PHR length view with 500 kb cap |
| 10 | 06b | Clade story matrix |
| 11 | 07m | Method transition: sequence sharing to communities |
| 12 | 07a | Crisp row-matched UPGMA heatmap/tree |
| 13 | 07b | Rooted readable NJ tree backup |
| 14 | 07c | Unrooted NJ with acrocentric p-arm audit |
| 15 | 08m | Method transition: MDS / PCoA and pairwise population variation |
| 16 | 08a | Backup chromosome-colored MDS / PCoA |
| 17 | 08b | Within-population pairwise variation metric |
| 18 | 09 | 1:1 all-community MDS / PCoA labels |
| 19 | 10m | Method transition: explicit 3D contact validation language |
| 20 | 10a | Axis/box-corrected Pore-C community matrix |
| 21 | 10b | Mantel exclusion clarity plot |
| 22 | 11 | Explicit within-community vs between-community single-cell distance plot |
| 23 | 12 | Mouse zygotene plus stage trajectory |
| 24 | 12b | Human sequence-similarity vs 3D-contact arm-pair analog |
| 25 | 13a | Pedigree proof, top readability crop |
| 26 | 13b | Pedigree details, bottom readability crop |
| 27 | 14m | Method transition: copy-aware gene cargo |
| 28 | 14a | DUX4/D4Z4 C1 genome-browser panel |
| 29 | 14b | OR4F-rich C3 genome-browser panel |
| 30 | 14c | OR4F pseudogene endpoint C8 genome-browser panel |
| 31 | 14d | TAR1-rich C2 genome-browser panel |
| 32 | 14e | C7 acrocentric p-arm genome-browser panel |
| 33 | 15 | Closing text-only focus page |

## Validation

- Typst compile succeeded for `BoG_2026_review_zoom_v3.pdf`.
- PNG export succeeded for `_typst/page-01.png` through `_typst/page-33.png`.
- `strings BoG_2026_review_zoom_v3.pdf | rg -c '^  /Type /Page$'` returned `33`.
- `strings BoG_2026_review_zoom_v3.pdf | rg -c '^  /MediaBox \[0 0 959\.76 540\]$'` returned `33`.
- The PDF page size is therefore 959.76 x 540 pt, matching 13.33 in x 7.5 in and 16:9.
- `file _typst/page-01.png _typst/page-33.png` reports 1920 x 1080 RGBA PNGs; `ls _typst/page-*.png | wc -l` reports `33`.
- Standard-library PNG pixel scans found representative corrected pages nonblank: 09, 12, 14, 17, 18, 20, 22, 24, 28, 32, and 33.
- A stale prior-agent worktree path scan over `slides/v2-review-zoom` returns no matches.
- `git diff --check` passes.
