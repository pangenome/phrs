// Standalone preview pages for review-zoom-v8-hic-dipc-clarity-split.
// These are handoff assets only. Do not compile this file into the final deck.

#set page(
  width: 13.33in,
  height: 7.5in,
  margin: (x: 0.18in, y: 0.12in),
)
#set text(font: "DejaVu Sans", size: 10pt, lang: "en")
#set par(justify: false, leading: 0.62em)

#let col-title = rgb("#1a3a6b")
#let col-hdr-bg = rgb("#dce8f7")
#let col-text = rgb("#222222")
#let col-cap = rgb("#6b6b6b")
#let col-note-bg = rgb("#fff8e8")
#let col-note-bar = rgb("#d4820a")

#let header(num, label) = block(
  fill: col-hdr-bg,
  width: 100%,
  inset: (x: 0.08in, y: 0.035in),
  radius: 2pt,
)[
  #text(size: 8pt, fill: col-title.lighten(15%), weight: "semibold")[
    Slide #num #h(0.12in) #label
  ]
]

#let footer(source) = align(center)[#text(size: 5.6pt, fill: col-cap)[#source]]

#let note-box(body) = block(
  width: 100%,
  fill: col-note-bg,
  stroke: (left: 3pt + col-note-bar),
  inset: (x: 0.12in, y: 0.06in),
  radius: 2pt,
)[#text(size: 9.4pt, fill: col-text)[#body]]

#let caption-box(body) = block(
  width: 100%,
  inset: (x: 0.06in, y: 0.02in),
)[#text(size: 8.4pt, fill: col-text)[#body]]

#let preview-slide(num, label, title, path, caption, source) = {
  grid(
    rows: (0.34in, 0.34in, 1fr, 0.58in, 0.13in),
    row-gutter: 0.035in,
    align: center,
    header(num, label),
    align(left)[#text(size: 15pt, fill: col-title, weight: "bold")[#title]],
    box(width: 100%, height: 100%)[
      #image(path, width: 100%, height: 100%, fit: "contain")
    ],
    caption-box(caption),
    footer(source),
  )
}

#grid(
  rows: (0.34in, 0.34in, 1fr, 0.78in, 0.13in),
  row-gutter: 0.035in,
  align: center,
  header("10m.2", "Hi-C contact-space MDS"),
  align(left)[#text(size: 15pt, fill: col-title, weight: "bold")[CHM13 Hi-C MDS is whole-arm contact space]],
  box(width: 100%, height: 100%)[
    #image("inputs/chm13_hic_mds_3d_coords.png", width: 100%, height: 100%, fit: "contain")
  ],
  note-box[
    #text(weight: "bold")[Region definition:] whole CHM13 p-arm, centromere, and q-arm intervals from `CHM13_chrom_parts.bed` (72 regions). This is #text(weight: "bold")[not] PHR intervals, #text(weight: "bold")[not] terminal 500 kb, and #text(weight: "bold")[not] 1 Mb PHR flanks/windows. MDS embeds contact frequencies from a broader CHM13 `.mcool` arm/centromere matrix; it is a contact-space summary, not a physical single-cell genome reconstruction.
  ],
  footer("Source: /moosefs/guarracino/HPRCv2/PHR_III/HiC/CHM13/plots/MDS_3d_coords.png; analyzer.py + CHM13_chrom_parts.bed"),
)

#pagebreak()

#preview-slide(
  "10m.2b",
  "Optional physical-coordinate context",
  "One existing GM12878 Dip-C cell gives a whole-genome 3DG view",
  "plots/gm12878_cell01_whole_genome_3dg_projection.png",
  [Optional candidate only: this is a 2D x/y projection of one existing 3DG file (`gm12878_01`) from single-cell Dip-C. It is physical-coordinate context, not CHM13 Hi-C and not a PHR/community validation plot.],
  "Dataset: GM12878 | Technology: Dip-C / 3DG | n=1 cell | Source: /moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/3dg/gm12878_01.impute3.round4.clean.3dg.gz",
)

#pagebreak()

#preview-slide(
  "11a.1",
  "GM12878 Dip-C proximity",
  "GM12878 Dip-C: sequence similarity predicts 3D proximity",
  "plots/gm12878_mantel_proximity.png",
  [W/B means within-community divided by between-community 3D distance; W/B < 1 means within-community arms are closer. Proximity convention: `3D proximity = -mean 3D distance`, so higher y-values are closer and a positive rho has the expected direction.],
  "Dataset: GM12878 | Technology: Dip-C | n=16 cells | source TSVs: /moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_mantel_3d.tsv; gm12878_arm_3d_distance_matrix.tsv; /moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv",
)

#pagebreak()

#preview-slide(
  "11a.2",
  "GM12878 Dip-C radial",
  "GM12878 Dip-C radial community structure",
  "inputs/gm12878_radial_community.png",
  [This radial panel is GM12878 only: per-community normalized radial position plus within-vs-between radial similarity. Same-community arms have more similar radial positions than between-community arms.],
  "Dataset: GM12878 | Technology: Dip-C | n=16 cells | source PDF/TSV: /moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_radial_community.pdf; gm12878_radial_community.tsv; gm12878_summary.tsv",
)

#pagebreak()

#preview-slide(
  "11a.3",
  "Sperm scHi-C proximity",
  "Sperm single-cell Hi-C: sequence similarity predicts 3D proximity",
  "plots/sperm_all20_mantel_proximity.png",
  [W/B means within-community divided by between-community 3D distance; W/B < 1 means within-community arms are closer. Proximity convention: `3D proximity = -mean 3D distance`, so higher y-values are closer and a positive rho has the expected direction.],
  "Dataset: human sperm | Technology: single-cell Hi-C / 3DG | n=20 cells | source TSVs: /moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_mantel_3d.tsv; sperm_all20_arm_3d_distance_matrix.tsv; /moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv",
)

#pagebreak()

#preview-slide(
  "11a.4",
  "Sperm scHi-C radial",
  "Sperm single-cell Hi-C radial community structure",
  "inputs/sperm_all20_radial_community.png",
  [This radial panel is sperm only: per-community normalized radial position plus within-vs-between radial similarity across 20 haploid sperm cells. It measures the same radial statistic as the GM12878 radial panel, in a distinct cell type.],
  "Dataset: human sperm | Technology: single-cell Hi-C / 3DG | n=20 cells | source PDF/TSV: /moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_radial_community.pdf; sperm_all20_radial_community.tsv; sperm_all20_summary.tsv",
)

#pagebreak()

#preview-slide(
  "11b",
  "Adjusted text scale",
  "Negative control: non-sharing S_all arms are farther apart",
  "plots/wb_negative_control_reduced_text.png",
  [Reduced text scale version of slide 11b. W/B is within-community divided by between-community 3D distance; values below 1 mean closer within-community, while S_all values above 1 show the non-sharing negative control is farther apart.],
  "source TSV roots: /moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50 and /moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected; files: *_per_cell.tsv and *_per_community_per_cell.tsv",
)

#pagebreak()

#preview-slide(
  "11c",
  "Adjusted text scale",
  "Community-free per-cell rho: sequence similarity predicts proximity",
  "plots/community_free_rho_distribution_reduced_text.png",
  [Reduced text scale version of slide 11c. Rho is computed as sequence similarity versus 3D proximity (`-distance`), so positive rho means more similar sequences are closer in 3D.],
  "source TSV roots: /moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50 and /moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected; files: *_community_free_per_cell.tsv and *_community_free_arm.tsv",
)
