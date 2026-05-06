# BoG 2026 Review Zoom Deck v2 Revision Notes

Task: `render-bog-annotated-zoom-v2-review-fanin`
Date: 2026-05-06 UTC
Scope: `slides/v2-review-zoom`

This revision builds the next review-zoom deck as
`BoG_2026_review_zoom_v2.pdf`. The original current zoom render,
`BoG_2026_review_zoom.pdf`, is preserved as the v1 reference.

The upstream provenance audit was used as the first guide:
`slides/v2-review-zoom/_revision_assets/git_provenance/README.md`. That audit
identifies the current deck as agent-951 commit `10bee88`, squash-merged as
`4862ec7`, and warns that many visible assets are copied/cropped PNGs rather
than reproducible renders. This v2 deck keeps new fanout candidates referenced
from `_revision_assets/...` and records the source path, generating script,
external URL/commit, or fanout commit for every major changed asset below.

## Build Notes

The revised Typst source is:

- `slides/v2-review-zoom/_typst/zoom_review_deck.typ`

Because the v2 Typst source intentionally references sibling provenance folders
such as `../_revision_assets/03b_erdos_renyi/...`, compile with the review-zoom
directory as the Typst root:

```bash
cd slides/v2-review-zoom/_typst
typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v2.pdf
typst compile --root .. --ppi 144 zoom_review_deck.typ page-{0p}.png
```

## Addressed Review Points

| Review area | v2 deck action | Status |
|---|---|---|
| 03b Erdos-Renyi callout | Replaced the weak local two-bar callout with `erdos_renyi_connectivity_candidate.png`, using the external `github.com:ekg/erdos_renyi` threshold logic and the HPRC point (`n = 18,827`, `p ~ 0.12`). | Addressed on slide `03b`. |
| 04 HPRCv2 karyogram | Added the HPRCv2 interchromosomal karyogram as the main slide `04`. Demoted the old Fig 1a heatmap crop to backup slide `04b` and the older count view to backup slide `05`. | Addressed while preserving the old reference views. |
| 06 length distribution density | Replaced the four-crop dense 06a split with the ranked arm summary candidate. | Addressed on slide `06a`. |
| 06 clade callout story | Replaced the old clade callout with the clade story matrix that ties length evidence to named biological/community groups. | Addressed on slide `06b`. |
| 07a heatmap order and p/q labels | Replaced Fig 1c crop with the UPGMA tree-left candidate. The slide source/caption language says UPGMA ordering, p/q label coloring, and Leiden k=15 ticks. | Addressed on slide `07a`. |
| 07b tree completeness / unrooted concern | Added the rooted readable NJ tree as backup `07b` and the unrooted NJ option as backup `07c`. | Addressed, including root-sensitivity option. |
| 08b superpopulation dispersion | Added the RMS-radius dispersion plot. The on-deck framing says population spread is not an ancestry axis. | Addressed on slide `08b`. |
| 09 MDS/community cleanup | Replaced the three-slide PCA/community duplication with one direct-labeled MDS/PCoA community plot. | Addressed on slide `09`. |
| Hi-C methods glossary | Added a 3D validation method transition before the Hi-C block. Detailed O/E, mcool, Mantel, Pore-C, Dip-C, and mouse details are kept as speaker-note material, not separate dense slides. | Addressed on slide `10m` plus speaker-note-only list below. |
| Hi-C visual redesign | Replaced 10a, 10b, 11, and the 12a/12b pair with the four redesigned Hi-C candidate visuals. | Addressed on slides `10a`, `10b`, `11`, and `12`. |
| DUX4/OR4F/TAR1/gene enrichment | Added a gene-cargo method transition, kept DUX4 and TAR1 as backup/context visuals, and replaced OR4F with the regenerated OR4F/gene-family signal candidate. | Addressed on slides `14m`, `14a`, `14b`, `14c`. |
| Density/methodology flow | Added short transition cards before the method, sequence-community, MDS/population, 3D-contact, and gene-cargo blocks. Demoted legacy dense views to backup labels rather than letting them carry the main talk path. | Addressed across slides `02m`, `07m`, `08m`, `10m`, `14m`. |

## Speaker-Note-Only Or Remaining Caveats

These review points are better handled verbally or in backup notes than as more
visible text on already dense slides:

- Erdos-Renyi caveat: wfmash pair selection is k-mer-selective, not literally
  random. The ER curve is a connectivity sanity check, not a biological null.
- HPRCv2 karyogram caveat: the PNG is the better count-landscape view, but it is
  not a 1:1 replacement for the old Fig 1a identity heatmap claim. The old Fig
  1a crop is therefore preserved as backup `04b`.
- Slide 06 caveat: the complete per-arm histogram audit remains useful, but the
  old four-crop dense view is not kept in the spoken path. The candidate assets
  and TSVs preserve the per-arm values for Q&A.
- Slide 07 caveat: do not call the 07a heatmap ordering "NJ". The 07a candidate
  is UPGMA average-linkage; the NJ topology is shown separately on 07b/07c.
- MDS/population caveat: the unit of analysis is a sequence-level PHR flank, not
  one individual, and points are not independent. The displayed 2D MDS/PCoA does
  not support a claim that AFR is unusually dispersed.
- Hi-C glossary details: O/E normalization, `.mcool` storage, Mantel rho, Pore-C,
  Dip-C, sperm, and mouse-stage details should be spoken only when that slide is
  active. The deck keeps the main visual claim uncluttered.
- Gene enrichment caveat: OR4F is a clean gene-family visual, but the
  community-family Fisher rows are not BH-significant. Present OR4F/DUX4/TAR1 as
  copy-aware cargo and markers carried by PHR exchange, not as causal drivers.
- Provenance caveat: retained backup crops from the current review-zoom deck
  still have unknown crop geometry. They are labeled as backup/reference pages
  and are not treated as regenerated analysis products.

## Major Asset Provenance

| Slide(s) | Asset path used by v2 deck | Source / generator / commit |
|---|---|---|
| All | `_typst/zoom_review_deck.typ` | Revised in this fan-in task. Based on current agent-951 deck lineage (`10bee88`, squash `4862ec7`) and the fanout recommendations listed here. |
| `02`, `03a` | `_typst/assets/s02_interval_tree.png`, `_typst/assets/s03_impg_workflow.png` | Retained current review-zoom assets. Provenance audit maps them to v1 page rasters copied through the v2 review asset lineage; extraction recipe is not committed. |
| `03b` | `_revision_assets/03b_erdos_renyi/erdos_renyi_connectivity_candidate.png` | Fanout commit `223eca4`. Generated by `_revision_assets/03b_erdos_renyi/make_03b_erdos_renyi_plot.R`. External source inspected: `https://github.com/ekg/erdos_renyi.git` at `d9ec48f1945d14f38d0131f56ed288cfc7883e73`, especially `erdos_renyi_viz.R` threshold definitions. |
| `04` | `_revision_assets/04_hprcv2_karyogram/p_interchrom_karyogram_count_rainbow_inset.100kb.png` | Fanout commit `f99cf6b`. Copied from `https://github.com/pangenome/HPRCv2` main at `d14883c314e683063abe8b461992f12825ccd5ed`. Relevant generator in HPRCv2: `scripts/plot-impg-coverage.inter-chr-map.R`; source data path is off-repo as documented in the fanout README. |
| `04b` | `_typst/assets/s04_fig1_panel_a.png` | Retained current review-zoom crop. Provenance audit maps full source to `paper_prep/figures/fig1/figure_fig1.png` / `figure_fig1.R`; crop recipe is not committed. |
| `05` | `_typst/assets/s05_interchrom.png` | Retained current review-zoom count-view asset. Provenance audit warns it is a manually prepared/resized derivative, not a newly regenerated v2 asset. |
| `06a`, `06b` | `_revision_assets/06_length_redesign/candidate_06a_ranked_arm_summary.png`, `_revision_assets/06_length_redesign/candidate_06b_clade_story_matrix.png` | Fanout commit `bdd8526`. Generated by `_revision_assets/06_length_redesign/make_06_length_redesign.R` from `arm_length_summary.tsv` and `clade_length_summary.tsv`; original HPRCv2 length-source path is documented in the fanout README. |
| `07a` | `_revision_assets/07a_heatmap_tree_pq/candidate_heatmap_upgma_tree_left_pq.png` | Fanout commit `49a1209`. Generated by `_revision_assets/07a_heatmap_tree_pq/make_candidate_heatmap.R`; ordering table `_revision_assets/07a_heatmap_tree_pq/candidate_upgma_tree_order.tsv`. Uses the canonical 41 x 41 arm-level Jaccard distance matrix and Leiden k=15 assignments documented in the README. |
| `07b`, `07c` | `_revision_assets/07b_tree_options/07b_rooted_acro_readable_large.png`, `_revision_assets/07b_tree_options/07b_unrooted_nj_option.png` | Fanout commit `ba3c52d`. Generated by `_revision_assets/07b_tree_options/make_07b_tree_options.R`. Audit tables: `matrix_audit.tsv`, `arm_presence.tsv`, `clade_recovery.tsv`; seven zero-signal arms are absent by construction. |
| `08a` | `_typst/assets/s08a_mds_chrom.png` | Retained current review-zoom backup asset. Provenance audit maps it to HPRCv2 off-tree `plot-similarity-subtelo.R` / MDS outputs, not to a true PCA. |
| `08b` | `_revision_assets/08b_superpop_dispersion/superpop_dispersion_rms_radius.png` | Fanout commit `979e50a`. Generated by `_revision_assets/08b_superpop_dispersion/make_superpop_dispersion.R`; backing metrics in `superpop_dispersion_metrics.tsv`. |
| `09` | `_revision_assets/09_mds_community/candidate_labeled_mds_community.png` | Fanout commit `4a0a55c`. Generated by `_revision_assets/09_mds_community/make_labeled_mds_community.R`; off-tree HPRCv2 sources documented as `hprcv2.1Mb.subtelo.full_mds.rds` and `hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`. |
| `10m` | Text transition | Draws from fanout commit `2c903c` (`hic_methods/README.md`) and visual-redesign fanout commit `f421645`. No new raster asset. |
| `10a` | `_revision_assets/hic_visual_redesign/slide_10a_square_matrix_candidate.png` | Fanout commit `f421645`. Generated by `_revision_assets/hic_visual_redesign/make_hic_visual_redesign.R` from HG002 Pore-C contact matrix and arm-community TSV sources documented in the README. |
| `10b` | `_revision_assets/hic_visual_redesign/slide_10b_mantel_exclusions_clarity.png` | Fanout commit `f421645`. Generated by `_revision_assets/hic_visual_redesign/make_hic_visual_redesign.R`; focused replacement for ED5 Mantel exclusion view. |
| `11` | `_revision_assets/hic_visual_redesign/slide_11_single_cell_purpose_candidate.png` | Fanout commit `f421645`. Generated by `_revision_assets/hic_visual_redesign/make_hic_visual_redesign.R` from Dip-C and sperm per-cell community enrichment TSVs documented in the README. |
| `12` | `_revision_assets/hic_visual_redesign/slide_12_mouse_zygotene_trajectory_pairing.png` | Fanout commit `f421645`. Generated by `_revision_assets/hic_visual_redesign/make_hic_visual_redesign.R`; combines mouse zygotene scatter and the four-stage trajectory previously split across 12a/12b. |
| `13a`, `13b` | `_typst/assets/s13_pedigree_top.png`, `_typst/assets/s13_pedigree_bottom.png` | Retained current review-zoom crops from `s13_pedigree.png`. Provenance audit maps the likely full-source lineage to the WashU untangle pedigree PDF, but raster/crop recipes are not committed. |
| `14m` | Text transition | Draws from fanout commits `cda3b92` (`14_gene_background`) and `aaac439` (`14_gene_enrichment_or4f`). No new raster asset. |
| `14a`, `14c` | `_typst/assets/s14_dux4.png`, `_typst/assets/s14_tar1.png` | Retained current review-zoom crops. Provenance audit maps full source to `slides/v2/slide_14_gene_biology.R` / `.png`; crop recipes from the earlier zoom split are not committed. |
| `14b` | `_revision_assets/14_gene_enrichment_or4f/or4f_gene_family_signal.png` | Fanout commit `aaac439`. Generated by `_revision_assets/14_gene_enrichment_or4f/make_or4f_gene_family_signal.R` from `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv` plus HPRCv2 community enrichment TSVs documented in the README. |
| `02m`, `07m`, `08m`, `10m`, `14m`, `15` | Typst text only | Flow and wording are synthesized from `narrative_density/README.md` fanout commit `5507949`, the slide-specific fanout READMEs above, and current v2 review deck wording. |

## v2 Page Map

| Page | Label | Role |
|---:|---|---|
| 01 | 01 | Title |
| 02 | 02m | Method transition: implicit graph |
| 03 | 02 | Backup interval-tree schematic |
| 04 | 03a | IMPG workflow |
| 05 | 03b | Erdos-Renyi connectivity threshold plot |
| 06 | 04 | HPRCv2 karyogram |
| 07 | 04b | Backup manuscript Fig 1a heatmap crop |
| 08 | 05 | Backup count-view interchromosomal sharing |
| 09 | 06a | Ranked PHR length summary |
| 10 | 06b | Clade story matrix |
| 11 | 07m | Sequence-community method transition |
| 12 | 07a | UPGMA heatmap with p/q labels |
| 13 | 07b | Rooted readable NJ tree backup |
| 14 | 07c | Unrooted NJ option backup |
| 15 | 08m | MDS/population method transition |
| 16 | 08a | Backup chromosome-colored MDS/PCoA |
| 17 | 08b | Superpopulation dispersion metric |
| 18 | 09 | Direct-labeled MDS/PCoA community plot |
| 19 | 10m | 3D contact validation method transition |
| 20 | 10a | Square Pore-C contact matrix |
| 21 | 10b | Mantel exclusion clarity plot |
| 22 | 11 | Single-cell 3D purpose plot |
| 23 | 12 | Mouse zygotene plus trajectory |
| 24 | 13a | Pedigree proof, top crop |
| 25 | 13b | Pedigree details, bottom crop |
| 26 | 14m | Gene-cargo method transition |
| 27 | 14a | DUX4/D4Z4 context |
| 28 | 14b | OR4F/gene-family signal |
| 29 | 14c | TAR1 repeat-context backup |
| 30 | 15 | Closing |

## Validation Summary

Detailed render and validation output is recorded in:

- `slides/v2-review-zoom/_typst/render.log`
- `slides/v2-review-zoom/_typst/VALIDATION.md`

The intended validation checks are:

- Typst compile succeeds for `BoG_2026_review_zoom_v2.pdf`.
- PNG export produces `page-01.png` through `page-30.png`.
- PDF page size remains 13.33 in x 7.5 in / 16:9.
- All exported pages are nonblank.
- Representative pages from the changed blocks are visually checked:
  `page-05`, `page-06`, `page-09`, `page-12`, `page-17`, `page-20`,
  `page-23`, `page-28`, and `page-30`.
- No stale `.wg-worktrees/agent-878` absolute paths are present under
  `slides/v2-review-zoom`.
