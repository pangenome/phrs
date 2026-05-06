# Review zoom 09: MDS/community recommendation

Task: `review-zoom-09-mds-community-leading`

## Recommendation

Use the community-colored MDS/PCoA as the primary slide and retire the PCA-labeled duplication.

Recommended spoken-deck disposition:

| Current review page | Current label/source | Recommendation | Final-facing label |
| --- | --- | --- | --- |
| 09a | `PCA communities`, `assets/s09_pca_communities.png` | Do not use as the primary community slide. It is a v1 raster with PCA wording and duplicates the actual MDS community story. Keep only as a backup/provenance comparison if needed. | `Backup: legacy v1 community projection` |
| 09b | `Clade legend`, `assets/s09_clade_legend.png` | Do not keep as a stand-alone spoken slide. Integrate the key rows into direct in-plot labels plus a compact side legend. Keep the full table as backup if the audience needs the C1-C15 map. | `Reference: Leiden C1-C15 clade map` |
| 09c | `Community-colored MDS`, `assets/s09b_communities.png` | Promote this content into the primary 09 slot. This is Erik's actual community-colored MDS slide; use it or the candidate direct-labeled version below. | `MDS / PCoA - named clades are Leiden communities` |

If final integration must preserve three review pages for traceability, order them as:

1. Primary: direct-labeled MDS/PCoA community plot.
2. Backup: full C1-C15 clade legend table.
3. Backup: legacy v1 PCA-labeled raster only if provenance comparison is needed.

The talk should not spend three consecutive slides on 09a/09b/09c. The community story should land once, directly on the plot, then move on.

## Candidate asset

Generated candidate:

- `candidate_labeled_mds_community.png`
- `candidate_labeled_mds_community.pdf`
- `make_labeled_mds_community.R`

Design decisions:

- The title and axes use `MDS / PCoA`, not `PCA`.
- Coordinates come from `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.full_mds.rds`.
- Arm-level communities come from `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`.
- The axis labels use the current cached MDS variance values: dimension 1 = 15.55%, dimension 2 = 10.80%.
- The six abstract-named clades are labeled directly on the plot:
  - C15 PAR1, Xp/Yp
  - C14 PAR2, Xq/Yq
  - C7 acrocentric p-arms, 13p/14p/15p/21p/22p
  - C2 10p-18p
  - C6 tight q-arm clade, 1q/13q/17q/19q/21q/22q
  - C1 DUX4/D4Z4, 4q/10q
- The six named communities use the slide-07 NJ-tree named-clade palette:
  - PAR1 red
  - PAR2 blue
  - acrocentric p-arms green
  - 10p-18p orange
  - tight q-arm clade purple
  - DUX4/D4Z4 brown
- Secondary Leiden communities are gray so they remain visible without competing with the abstract/community story.

This candidate intentionally replaces the stand-alone legend-only 09b with in-plot labels plus a compact side legend.

## Audit notes

Current review-zoom mapping shows the duplication plainly: 09a is `PCA communities`, 09b is `Clade legend`, and 09c is `Community-colored MDS` (`slides/v2-review-zoom/_typst/VALIDATION.md` lines 29-31; `slides/v2-review-zoom/_typst/zoom_review_deck.typ` lines 224-247).

The underlying method is MDS/PCoA, not strict PCA. Slide 08 already documents the source: `cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)` on the 15,668 x 15,668 Jaccard distance matrix, with variance labels derived from the MDS eigenvalues (`slides/v2/slide_08_pca_chromosome_superpop.md` lines 22 and 66-70). Final integration should keep that language consistent on slide 09.

The current slide-09 markdown still says `All-vs-all PCA`, `Same PCA layout`, `PC1 / PC2`, and recommends saying `PCA` verbally (`slides/v2/slide_09_pca_communities_clades.md` lines 1-8 and 119-134). Those labels should be rewritten as MDS/PCoA language in the final deck integration, unless someone generates a true PCA artifact from a feature matrix.

The cross-slide community vocabulary to preserve is the slide-07 NJ/Leiden/abstract map: C15 = PAR1, C14 = PAR2, C7 = acrocentric p-arms, C2 = 10p-18p, C6 = tight q-arm clade, and C1 = DUX4/D4Z4 (`slides/v2/slide_07_allvsall_heatmap_nj_clades.md` lines 86-92). Slide 13 then depends on those same Leiden communities for the 494/538 pedigree result and specific C1/C2/C3/C15 examples (`slides/v2/slide_13_pedigree_direct_evidence.md` lines 8-10). Slide 14 specifically uses C1 for the DUX4/D4Z4 biology (`slides/v2/slide_14_gene_biology.md` lines 7 and 49).

The old 09b legend uses light row fills from the slide-06/slide-09 callout palette. That is useful as a table treatment, but for the plotted points the candidate switches to the slide-07 named-clade palette so the color semantics match the NJ tree.

## Final integration checklist

- Replace `PCA communities` with `MDS / PCoA - named clades are Leiden communities`.
- Replace `PC1` / `PC2` with `MDS dimension 1` / `MDS dimension 2` or `PCoA dimension 1` / `PCoA dimension 2`.
- Use the community-colored MDS/PCoA as the primary slide 09 asset.
- Do not keep 09b as a separate spoken legend-only slide unless it is moved to backup.
- Keep the six named community labels and colors consistent with slide 07, and keep C1/C2/C3/C15 vocabulary available for the slide-13 pedigree handoff.
- Do not edit deck source here; final integration owns `zoom_review_deck.typ` and any deck-level renumbering.

## Validation

- Ran `Rscript slides/v2-review-zoom/_revision_assets/09_mds_community/make_labeled_mds_community.R`.
- Confirmed the script generated both PNG and PDF assets in this directory.
- Did not edit deck source.
