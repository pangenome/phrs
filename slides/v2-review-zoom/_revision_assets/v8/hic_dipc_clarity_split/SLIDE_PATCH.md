# Slide patch: v8 Hi-C / Dip-C clarity split

Do not edit the final deck in this task. This file is the handoff for
`review-zoom-v8-fanin-render`.

All asset paths below are relative to
`slides/v2-review-zoom/_typst/zoom_review_deck.typ`.

## 1. Add helper after `figure-slide`

Add this helper after the existing `#let figure-slide(...)` definition. It uses
the deck's existing `header`, `footer`, `col-*` variables.

```typst
#let captioned-figure-slide(num, label, title, path, caption, source: "") = {
  grid(
    rows: (0.34in, 0.34in, 1fr, 0.58in, 0.13in),
    row-gutter: 0.035in,
    align: center,
    header(num, label),
    align(left)[#text(size: 15pt, fill: col-title, weight: "bold")[#title]],
    box(width: 100%, height: 100%)[
      #image(path, width: 100%, height: 100%, fit: "contain")
    ],
    block(
      width: 100%,
      fill: col-note-bg,
      stroke: (left: 3pt + col-note-bar),
      inset: (x: 0.12in, y: 0.055in),
      radius: 2pt,
    )[#text(size: 8.6pt, fill: col-text)[#caption]],
    footer(source),
  )
}
```

## 2. Replace slide `10m.2`

Replace the current `#figure-slide("10m.2", ...)` block at
`zoom_review_deck.typ:683-688` with:

```typst
#captioned-figure-slide(
  "10m.2",
  "Hi-C MDS gives a 3D contact-space view",
  "CHM13 Hi-C MDS is whole-arm contact space",
  "../_revision_assets/v8/hic_dipc_clarity_split/inputs/chm13_hic_mds_3d_coords.png",
  [
    #text(weight: "bold")[Region definition:] whole CHM13 p-arm, centromere, and q-arm intervals from `CHM13_chrom_parts.bed` (72 regions). This is #text(weight: "bold")[not] PHR intervals, #text(weight: "bold")[not] terminal 500 kb, and #text(weight: "bold")[not] 1 Mb PHR flanks/windows. MDS embeds contact frequencies from a broader CHM13 `.mcool` arm/centromere matrix; it is a contact-space summary, not a physical single-cell genome reconstruction.
  ],
  source: "Dataset: CHM13 | Technology: Hi-C | source PNG: /moosefs/guarracino/HPRCv2/PHR_III/HiC/CHM13/plots/MDS_3d_coords.png | source BED: CHM13_chrom_parts.bed | script: analyzer.py",
)
```

Optional physical-coordinate context, only if the deck can afford one extra
slide immediately after `10m.2`:

```typst
#pagebreak()

#captioned-figure-slide(
  "10m.2b",
  "Optional physical-coordinate context",
  "One existing GM12878 Dip-C cell gives a whole-genome 3DG view",
  "../_revision_assets/v8/hic_dipc_clarity_split/plots/gm12878_cell01_whole_genome_3dg_projection.png",
  [
    Optional candidate only: this is a 2D x/y projection of one existing 3DG file (`gm12878_01`) from single-cell Dip-C. It is physical-coordinate context, not CHM13 Hi-C and not a PHR/community validation plot.
  ],
  source: "Dataset: GM12878 | Technology: Dip-C / 3DG | n=1 cell | source 3DG: /moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/3dg/gm12878_01.impute3.round4.clean.3dg.gz",
)
```

Recommendation: use the optional slide only if Erik specifically wants a real
physical-coordinate example. Otherwise keep `10m.2` as the CHM13 contact-space
MDS with the honest region/caveat caption.

## 3. Replace current slide `11a` with four readable slides

Delete/stop using the current `#dipc-validation-panel-slide("11a", ...)` block
at `zoom_review_deck.typ:719-723`. The helper `dipc-validation-panel-slide` may
remain defined; it will simply be unused.

Insert these four slides in its place, after the current slide `11` page and
before slide `11b`.

```typst
#captioned-figure-slide(
  "11a.1",
  "GM12878 Dip-C proximity",
  "GM12878 Dip-C: sequence similarity predicts 3D proximity",
  "../_revision_assets/v8/hic_dipc_clarity_split/plots/gm12878_mantel_proximity.png",
  [
    W/B means within-community divided by between-community 3D distance; W/B < 1 means within-community arms are closer. Proximity convention: `3D proximity = -mean 3D distance`, so higher y-values are closer and a positive rho has the expected direction.
  ],
  source: "Dataset: GM12878 | Technology: Dip-C | n=16 cells | source TSVs: gm12878_mantel_3d.tsv, gm12878_arm_3d_distance_matrix.tsv, hprcv2.1Mb.subtelo.arm_dist_matrix.tsv",
)

#pagebreak()

#captioned-figure-slide(
  "11a.2",
  "GM12878 Dip-C radial",
  "GM12878 Dip-C radial community structure",
  "../_revision_assets/v8/hic_dipc_clarity_split/inputs/gm12878_radial_community.png",
  [
    This radial panel is GM12878 only: per-community normalized radial position plus within-vs-between radial-position similarity. Same-community arms have more similar radial positions than between-community arms.
  ],
  source: "Dataset: GM12878 | Technology: Dip-C | n=16 cells | source PDF/TSV: /moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_radial_community.pdf; gm12878_radial_community.tsv; gm12878_summary.tsv",
)

#pagebreak()

#captioned-figure-slide(
  "11a.3",
  "Sperm scHi-C proximity",
  "Sperm single-cell Hi-C: sequence similarity predicts 3D proximity",
  "../_revision_assets/v8/hic_dipc_clarity_split/plots/sperm_all20_mantel_proximity.png",
  [
    W/B means within-community divided by between-community 3D distance; W/B < 1 means within-community arms are closer. Proximity convention: `3D proximity = -mean 3D distance`, so higher y-values are closer and a positive rho has the expected direction.
  ],
  source: "Dataset: human sperm | Technology: single-cell Hi-C / 3DG | n=20 cells | source TSVs: sperm_all20_mantel_3d.tsv, sperm_all20_arm_3d_distance_matrix.tsv, hprcv2.1Mb.subtelo.arm_dist_matrix.tsv",
)

#pagebreak()

#captioned-figure-slide(
  "11a.4",
  "Sperm scHi-C radial",
  "Sperm single-cell Hi-C radial community structure",
  "../_revision_assets/v8/hic_dipc_clarity_split/inputs/sperm_all20_radial_community.png",
  [
    This radial panel is sperm only: per-community normalized radial position plus within-vs-between radial-position similarity across 20 haploid sperm cells. It measures the same radial statistic as the GM12878 radial panel, in a distinct cell type.
  ],
  source: "Dataset: human sperm | Technology: single-cell Hi-C / 3DG | n=20 cells | source PDF/TSV: /moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_radial_community.pdf; sperm_all20_radial_community.tsv; sperm_all20_summary.tsv",
)
```

This replacement preserves all four concepts from the old tiny multipanel page,
but each panel now occupies its own slide. The old right-side radial plots are
separated and labeled explicitly as GM12878 versus sperm.

## 4. Replace slide `11b`

Replace the current slide `11b` block at `zoom_review_deck.typ:727-732` with:

```typst
#captioned-figure-slide(
  "11b",
  "Negative control: non-sharing S_all arms are farther apart",
  "Negative control: non-sharing S_all arms are farther apart",
  "../_revision_assets/v8/hic_dipc_clarity_split/plots/wb_negative_control_reduced_text.png",
  [
    Reduced text scale. W/B is within-community divided by between-community 3D distance; values below 1 mean closer within-community, while S_all values above 1 show the non-sharing negative control is farther apart.
  ],
  source: "source TSV roots: /moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50 and /moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected; files: *_per_cell.tsv and *_per_community_per_cell.tsv",
)
```

## 5. Replace slide `11c`

Replace the current slide `11c` block at `zoom_review_deck.typ:736-741` with:

```typst
#captioned-figure-slide(
  "11c",
  "Community-free per-cell rho: sequence similarity predicts proximity",
  "Community-free per-cell rho: sequence similarity predicts proximity",
  "../_revision_assets/v8/hic_dipc_clarity_split/plots/community_free_rho_distribution_reduced_text.png",
  [
    Reduced text scale. Rho is computed as sequence similarity versus 3D proximity (`-distance`), so positive rho means more similar sequences are closer in 3D.
  ],
  source: "source TSV roots: /moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50 and /moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected; files: *_community_free_per_cell.tsv and *_community_free_arm.tsv",
)
```

## 6. Preview assets

Standalone previews are available if the fan-in render wants to inspect the
layout before patching the deck:

- `slides/hic_dipc_clarity_split_preview.pdf`
- `slides/slide10m2_chm13_hic_contact_mds.png`
- `slides/slide10m2b_gm12878_whole_genome_3dg_candidate.png`
- `slides/slide11a1_gm12878_mantel_proximity.png`
- `slides/slide11a2_gm12878_radial_community.png`
- `slides/slide11a3_sperm_mantel_proximity.png`
- `slides/slide11a4_sperm_radial_community.png`
- `slides/slide11b_negative_control_reduced_text.png`
- `slides/slide11c_community_free_rho_reduced_text.png`

## 7. Source manifest

Use `Source_manifest.tsv` for panel-level provenance. It includes dataset,
technology, cell/sample count, source TSV/PDF paths, region definitions, and
the distance/proximity convention for every staged panel.
