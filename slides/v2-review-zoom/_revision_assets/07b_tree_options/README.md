# Slide 07b Tree Audit And Options

## Short Answer

Slide 07b is not missing any of the 41 expected arm-level subtelomere tips. The current tree is visually small because it is an arm-level, 41-tip tree rendered as a 1600 x 1300 PNG and then fit into a single zoom-review slide. It is not a sequence-level tree of 15,668 PHR sequences or 18,827 telomeric flanks.

The correct data interpretation is:

- 41 PHR-positive chromosome arms are present in the 41 x 41 arm-level distance matrix, the Leiden k=15 assignments, the current NJ tree Newick, and the new candidate tree options.
- 7 of the full 48 human chromosome arms are absent by construction because they have no retained inter-chromosomal PHR signal at the analysis threshold: `chr2_p`, `chr3_p`, `chr5_p`, `chr8_q`, `chr11_q`, `chr14_q`, and `chr18_q`.
- The existing slide planning text that says only six arms are missing is stale. It omits `chr18_q`. The slide 09 legend text that lists `18q (n=1), Xp, Yp` for C15 is also stale against the current matrix/assignment data; current C15 is only `chrX_p`, `chrY_p`.

## Files Produced

All files in this directory were generated without editing the deck source.

| File | Use |
|---|---|
| `07b_rooted_acro_readable_large.png` / `.pdf` | Larger rooted tree option. Same NJ topology, rooted at the C7 acrocentric p-arm MRCA for display continuity with the current slide. |
| `07b_unrooted_nj_option.png` / `.pdf` | Unrooted NJ option. This is defensible because NJ produces an unrooted distance tree and there is no independent biological outgroup among human subtelomeres. |
| `07b_named_clade_legend.png` / `.pdf` | Separate clade legend using the slide 09 community IDs and abstract vocabulary. |
| `07b_rooted_acro.newick` | Rooted-display Newick for the large rooted option. |
| `07b_unrooted_nj.newick` | Unrooted NJ Newick directly from `ape::nj()`. |
| `arm_presence.tsv` | Exhaustive 48-arm presence audit. |
| `matrix_audit.tsv` | Matrix shape/range audit. |
| `clade_legend.tsv` | Community-to-abstract vocabulary map used for color/highlight assignment. |
| `clade_recovery.tsv` | Monophyly check for the six named slide 09 / abstract clades. |
| `make_07b_tree_options.R` | Reproducible renderer for all assets above. |

## Current Asset Audit

Current slide 07b uses `slides/v2-review-zoom/_typst/assets/s07b_nj_tree.png`. The provenance audit reports that this is an exact blob match to `paper_prep/figures/nj_tree_arms/nj_tree_annotated.png`, generated from `paper_prep/figures/nj_tree_arms/nj_tree.R`.

The upstream NJ script documents the same source matrix and output intent at `paper_prep/figures/nj_tree_arms/nj_tree.R:1-16`. It reads `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`, symmetrizes it, zeros the diagonal for NJ, and runs `ape::nj()` at `paper_prep/figures/nj_tree_arms/nj_tree.R:33-49`. It defines the six highlighted abstract clades at `paper_prep/figures/nj_tree_arms/nj_tree.R:51-70`, roots the display at the acrocentric p-arm MRCA at `paper_prep/figures/nj_tree_arms/nj_tree.R:73-100`, and uses a perturbation bootstrap rather than a character-level bootstrap at `paper_prep/figures/nj_tree_arms/nj_tree.R:102-148`.

The current PNG is therefore source-correct, but it is too compressed for slide viewing. It tries to show all 41 tip labels, all internal bootstrap labels, clade coloring, and a legend in one raster asset. The result can look like too few subtelomeres are shown, but the issue is readability, not missing tips.

## Source Data

Primary inputs:

- Distance matrix: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`
- Leiden assignments: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
- Upstream current-tree script: `paper_prep/figures/nj_tree_arms/nj_tree.R`
- Slide 09 legend script: `slides/v2/_typst/slide_09_clade_legend.R`

The report source says the arm-level community step groups 41 arms and excludes 7 zero-signal arms at `end-to-end-report/report/01_pipeline.md:122-126`. It gives the full current Leiden community table at `end-to-end-report/report/01_pipeline.md:133-151`. It documents the `chr18_q` artifact and removal at `end-to-end-report/report/01_pipeline.md:75-85`, and the seven zero-signal arms at `end-to-end-report/report/01_pipeline.md:87-102`.

Matrix audit from `matrix_audit.tsv`:

| Metric | Value |
|---|---:|
| Rows | 41 |
| Columns | 41 |
| Present arms | 41 |
| Missing from 48-arm set | `chr2_p`, `chr3_p`, `chr5_p`, `chr8_q`, `chr11_q`, `chr14_q`, `chr18_q` |
| Raw diagonal range | 0.009725 to 0.597977 |
| Raw max asymmetry | 0 |
| Off-diagonal range after zeroing diagonal | 0.015324 to 1.000000 |

The nonzero raw diagonal is expected for this source matrix because diagonal cells carry within-arm average distances. For NJ, self-distance must be zero, so both the upstream script and this renderer set `diag(D) <- 0` before `ape::nj()`.

## Arm Presence

All 41 expected PHR-positive arms are present. Grouped by current Leiden k=15 community:

| Community | Present arms |
|---|---|
| C1 | `chr4_q`, `chr10_q` |
| C2 | `chr10_p`, `chr18_p` |
| C3 | `chr3_q`, `chr7_p`, `chr9_q`, `chr11_p`, `chr16_q`, `chr19_p` |
| C4 | `chr7_q`, `chr12_q` |
| C5 | `chr6_p`, `chr9_p`, `chr12_p`, `chr20_q` |
| C6 | `chr1_q`, `chr13_q`, `chr17_q`, `chr19_q`, `chr21_q`, `chr22_q` |
| C7 | `chr13_p`, `chr14_p`, `chr15_p`, `chr21_p`, `chr22_p` |
| C8 | `chr15_q` |
| C9 | `chr16_p` |
| C10 | `chr17_p` |
| C11 | `chr1_p`, `chr5_q`, `chr6_q`, `chr8_p` |
| C12 | `chr2_q`, `chr20_p` |
| C13 | `chr4_p` |
| C14 | `chrX_q`, `chrY_q` |
| C15 | `chrX_p`, `chrY_p` |

Absent from the 48-arm human arm set:

| Arm | Status | Reason |
|---|---|---|
| `chr2_p` | absent by construction | Zero inter-chromosomal PHR signal at the >=95% identity PHR threshold. |
| `chr3_p` | absent by construction | Zero inter-chromosomal PHR signal at the >=95% identity PHR threshold. |
| `chr5_p` | absent by construction | Zero inter-chromosomal PHR signal at the >=95% identity PHR threshold. |
| `chr8_q` | absent by construction | Zero inter-chromosomal PHR signal at the >=95% identity PHR threshold. |
| `chr11_q` | absent by construction | Zero inter-chromosomal PHR signal at the >=95% identity PHR threshold. |
| `chr14_q` | absent by construction | Zero inter-chromosomal PHR signal at the >=95% identity PHR threshold. |
| `chr18_q` | absent by construction | Zero retained inter-chromosomal PHR signal after the NA18982#1 `chr18_q` / `chrX` PAR1 scaffold chimera was removed. |

## Tree Algorithm

The generated options use the same algorithmic inputs as the current slide 07b tree:

1. Read the 41 x 41 arm-level Jaccard distance matrix.
2. Symmetrize the matrix defensively.
3. Set the diagonal to zero for valid distance-tree construction.
4. Build the NJ topology with `ape::nj(as.dist(D))`.
5. Write the unrooted Newick directly from `ape::nj()`.
6. For the rooted display only, root the same topology at the C7 acrocentric p-arm MRCA because C7 is monophyletic.

The unrooted option is methodologically cleaner if the slide is meant to avoid implying ancestry or an outgroup. The rooted option is still acceptable as a display layout if its caption says "rooted at C7 acrocentric p-arm MRCA for orientation only."

## Clade Labeling

The candidate assets highlight the same six slide 09 / abstract clades:

| Community | Abstract vocabulary | Members | NJ monophyletic |
|---|---|---|---|
| C1 | 4q-10q DUX4-containing homology | `chr4_q`, `chr10_q` | yes |
| C2 | 10p-18p homology | `chr10_p`, `chr18_p` | yes |
| C6 | tightly linked 22q/21q/19q/1q/13q/17q clade | `chr1_q`, `chr13_q`, `chr17_q`, `chr19_q`, `chr21_q`, `chr22_q` | yes |
| C7 | acrocentric short arms | `chr13_p`, `chr14_p`, `chr15_p`, `chr21_p`, `chr22_p` | yes |
| C14 | Xq/Yq via PAR2 | `chrX_q`, `chrY_q` | yes |
| C15 | Xp/Yp via PAR1 | `chrX_p`, `chrY_p` | yes |

This follows the intended highlighted rows in the slide 09 legend script at `slides/v2/_typst/slide_09_clade_legend.R:7-42`, but with current membership from the actual assignments file. Do not carry forward `18q (n=1)` in C15 unless a separate, newer assignment table is introduced; it is not present in the matrix audited here.

I did not add a seventh highlight for the abstract phrase "large moderate-similarity clade" because slide 09 does not define a specific corresponding community row. If that phrase is retained in the talk, treat it as a backbone-level narrative observation, not as a labeled clade unless a specific arm set is defined.

## Recommendation

Use one tree, not the current small side-by-side heatmap plus tree, if the goal is to answer the reviewer concern at conference distance.

Preferred slide 07b replacement:

- Use `07b_rooted_acro_readable_large.png` plus `07b_named_clade_legend.png`.
- Caption: "NJ tree on the same 41 x 41 arm-level Jaccard distance matrix; all 41 PHR-positive arms are present. Seven zero-signal arms are absent by construction."
- Add one small note if space allows: "Display-rooted at acrocentric p-arm MRCA; topology is NJ from the arm-level distance matrix."

Alternative if the reviewer is sensitive to rooting:

- Use `07b_unrooted_nj_option.png` plus `07b_named_clade_legend.png`.
- Caption: "Unrooted NJ tree on the 41 x 41 arm-level Jaccard distance matrix; no outgroup or ancestry implied."

Avoid repeating the older "six missing arms" callout. The correct missing set has seven arms and includes `chr18_q`.

## Validation

- Rebuilt assets by running `Rscript slides/v2-review-zoom/_revision_assets/07b_tree_options/make_07b_tree_options.R`.
- Verified the generated matrix audit reports 41 rows, 41 columns, and seven absent arms from the 48-arm set.
- Verified all six highlighted slide 09 / abstract clades are monophyletic on the NJ topology in `clade_recovery.tsv`.
- Verified generated PNG dimensions with `file`: rooted asset is 3840 x 2160; unrooted asset is 3200 x 3200; legend asset is 2700 x 1350.
- No files under `slides/v2-review-zoom/_typst/` or any deck source were edited.
