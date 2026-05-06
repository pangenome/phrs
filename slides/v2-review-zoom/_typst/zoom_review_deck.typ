// BoG 2026 review zoom deck, v2 review revision.
// Layout: 16:9 widescreen, one visual focus per page.
// Build:
//   typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v2.pdf
//   typst compile --root .. --ppi 144 zoom_review_deck.typ page-{0p}.png

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
#let col-cap = rgb("#777777")
#let col-note-bg = rgb("#fff8e8")
#let col-note-bar = rgb("#d4820a")
#let col-card-bg = rgb("#f5f8fc")

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

#let footer(source) = {
  if source != "" {
    align(center)[#text(size: 5.6pt, fill: col-cap)[#source]]
  } else {
    align(center)[#text(size: 5.6pt, fill: col-cap)[review zoom v2 focus page]]
  }
}

#let figure-slide(num, label, path, source: "") = {
  grid(
    rows: (0.34in, 1fr, 0.13in),
    row-gutter: 0.02in,
    align: center,
    header(num, label),
    box(width: 100%, height: 100%)[
      #image(path, width: 100%, height: 100%, fit: "contain")
    ],
    footer(source),
  )
}

#let quad-slide(num, label, a, b, c, d, source: "") = {
  grid(
    rows: (0.34in, 1fr, 0.13in),
    row-gutter: 0.02in,
    align: center,
    header(num, label),
    box(width: 100%, height: 100%)[
      #grid(
        columns: (1fr, 1fr),
        rows: (3.27in, 3.27in),
        column-gutter: 0.08in,
        row-gutter: 0.06in,
        align: center,
        box(width: 100%, height: 3.27in)[#image(a, width: 100%, height: 100%, fit: "contain")],
        box(width: 100%, height: 3.27in)[#image(b, width: 100%, height: 100%, fit: "contain")],
        box(width: 100%, height: 3.27in)[#image(c, width: 100%, height: 100%, fit: "contain")],
        box(width: 100%, height: 3.27in)[#image(d, width: 100%, height: 100%, fit: "contain")],
      )
    ],
    footer(source),
  )
}

#let text-slide(num, label, title, body, thesis: none, source: "") = {
  grid(
    rows: (0.34in, 1fr, 0.13in),
    row-gutter: 0.02in,
    align: center,
    header(num, label),
    align(center + horizon)[
      #block(width: 10.8in)[
        #align(center)[
          #text(size: 28pt, weight: "bold", fill: col-title)[#title]
        ]
        #v(0.22in)
        #text(size: 15pt, fill: col-text)[#body]
        #if thesis != none {
          v(0.34in)
          block(
            fill: col-note-bg,
            inset: (x: 0.22in, y: 0.16in),
            radius: 4pt,
            stroke: (left: 4pt + col-note-bar),
            width: 100%,
          )[
            #text(size: 15pt, weight: "bold")[Thesis: ]
            #text(size: 15pt, style: "italic")[#thesis]
          ]
        }
      ]
    ],
    footer(source),
  )
}

#let method-slide(num, label, title, body, takeaway, source: "") = {
  grid(
    rows: (0.34in, 1fr, 0.13in),
    row-gutter: 0.02in,
    align: center,
    header(num, label),
    align(center + horizon)[
      #block(width: 11.1in)[
        #align(center)[
          #text(size: 23pt, weight: "bold", fill: col-title)[#title]
        ]
        #v(0.18in)
        #block(
          width: 100%,
          fill: col-card-bg,
          stroke: 0.8pt + col-hdr-bg.darken(15%),
          inset: (x: 0.24in, y: 0.18in),
          radius: 3pt,
        )[
          #text(size: 15pt, fill: col-text)[#body]
        ]
        #v(0.18in)
        #block(
          width: 100%,
          fill: col-note-bg,
          stroke: (left: 4pt + col-note-bar),
          inset: (x: 0.22in, y: 0.13in),
          radius: 3pt,
        )[
          #text(size: 14pt, weight: "bold", fill: col-text)[Takeaway: ]
          #text(size: 14pt, fill: col-text)[#takeaway]
        ]
      ]
    ],
    footer(source),
  )
}

#text-slide(
  "01",
  "Title - review focus page",
  "Concerted evolution and unorthodox recombination of human subtelomeres",
  [
    Andrea Guarracino and Erik Garrison. Biology of Genomes 2026.

    Companion to HPRC v2: population-scale inter-chromosomal subtelomeric relationships across 466 near-complete haplotypes.
  ],
)

#pagebreak()

#method-slide(
  "02m",
  "Method transition",
  "Method: an implicit graph over chromosome ends",
  [
    Pairwise telomere-anchored alignments are the graph.

    Interval trees make transitive closure queryable without building a GFA or partitioning by chromosome.

    The Erdos-Renyi check asks whether the sampled pair space is dense enough for genome-wide closure.
  ],
  [Sequence-sharing results downstream all come from this implicit pangenome graph over PHR flanks.],
  source: "narrative_density README; current slides 02/03a/03b collapsed into a methods transition",
)

#pagebreak()

#figure-slide(
  "02",
  "Backup: interval trees make alignments queryable",
  "assets/s02_interval_tree.png",
  source: "current review-zoom asset s02_interval_tree.png; provenance audit slide 02 row",
)

#pagebreak()

#figure-slide(
  "03a",
  "IMPG workflow",
  "assets/s03_impg_workflow.png",
  source: "current review-zoom asset s03_impg_workflow.png; provenance audit slide 03a row",
)

#pagebreak()

#figure-slide(
  "03b",
  "Why 12% sampling is enough for genome-wide closure",
  "../_revision_assets/03b_erdos_renyi/erdos_renyi_connectivity_candidate.png",
  source: "03b_erdos_renyi/make_03b_erdos_renyi_plot.R; external erdos_renyi d9ec48f",
)

#pagebreak()

#figure-slide(
  "04",
  "HPRCv2 interchrom karyogram: other chromosomes per 100 kb window",
  "../_revision_assets/04_hprcv2_karyogram/p_interchrom_karyogram_count_rainbow_inset.100kb.png",
  source: "HPRCv2 main d14883c; scripts/plot-impg-coverage.inter-chr-map.R; fanout copied PNG",
)

#pagebreak()

#figure-slide(
  "04b",
  "Backup: manuscript Fig 1a genome-wide identity heatmap",
  "assets/s04_fig1_panel_a.png",
  source: "current review-zoom crop from paper_prep/figures/fig1/figure_fig1.png; crop recipe unknown",
)

#pagebreak()

#figure-slide(
  "05",
  "Backup: genome-wide count view of interchromosomal sharing",
  "assets/s05_interchrom.png",
  source: "current review-zoom asset s05_interchrom.png; demoted after HPRCv2 karyogram",
)

#pagebreak()

#figure-slide(
  "06a",
  "PHR length hierarchy by chromosome arm",
  "../_revision_assets/06_length_redesign/candidate_06a_ranked_arm_summary.png",
  source: "06_length_redesign/make_06_length_redesign.R; arm_length_summary.tsv",
)

#pagebreak()

#figure-slide(
  "06b",
  "PHR length outliers define named biological communities",
  "../_revision_assets/06_length_redesign/candidate_06b_clade_story_matrix.png",
  source: "06_length_redesign/make_06_length_redesign.R; clade_length_summary.tsv",
)

#pagebreak()

#method-slide(
  "07m",
  "Method transition",
  "How sequence sharing becomes communities",
  [
    PHR alignments -> arm-level Jaccard distance matrix.

    Leiden communities, heatmaps, trees, and MDS/PCoA are complementary views of the same sequence-sharing structure.

    Labeling p arms and q arms separately prevents chromosome-end vocabulary from drifting.
  ],
  [The heatmap, tree, and MDS block should be read as one sequence-similarity system, not independent assays.],
  source: "narrative_density README; 07a_heatmap_tree_pq and 09_mds_community READMEs",
)

#pagebreak()

#figure-slide(
  "07a",
  "Arm-level similarity recovers named clades",
  "../_revision_assets/07a_heatmap_tree_pq/candidate_heatmap_upgma_tree_left_pq.png",
  source: "07a_heatmap_tree_pq/make_candidate_heatmap.R; UPGMA order, p/q labels, Leiden k=15 ticks",
)

#pagebreak()

#figure-slide(
  "07b",
  "Backup: tree confirms the same named clades",
  "../_revision_assets/07b_tree_options/07b_rooted_acro_readable_large.png",
  source: "07b_tree_options/make_07b_tree_options.R; NJ on 41x41 arm-level Jaccard matrix",
)

#pagebreak()

#figure-slide(
  "07c",
  "Backup: unrooted NJ option removes rooting implication",
  "../_revision_assets/07b_tree_options/07b_unrooted_nj_option.png",
  source: "07b_tree_options/make_07b_tree_options.R; unrooted NJ option for root-sensitive review",
)

#pagebreak()

#method-slide(
  "08m",
  "Method transition",
  "How to read MDS / PCoA and population spread",
  [
    MDS / PCoA places sequence-level PHR flanks by Jaccard distance.

    Dispersion is spread around a group centroid in the displayed panel, not an ancestry axis.

    The unit is a subtelomeric flank, so arm-community structure dominates the geometry.
  ],
  [Use MDS / PCoA language consistently unless a true feature-matrix PCA is generated.],
  source: "08b_superpop_dispersion README; 09_mds_community README",
)

#pagebreak()

#figure-slide(
  "08a",
  "Backup: MDS / PCoA colored by chromosome",
  "assets/s08a_mds_chrom.png",
  source: "current review-zoom asset; HPRCv2 similarity plot-similarity-subtelo.R off-tree source",
)

#pagebreak()

#figure-slide(
  "08b",
  "Population dispersion on the PHR MDS / PCoA",
  "../_revision_assets/08b_superpop_dispersion/superpop_dispersion_rms_radius.png",
  source: "08b_superpop_dispersion/make_superpop_dispersion.R; superpop_dispersion_metrics.tsv",
)

#pagebreak()

#figure-slide(
  "09",
  "MDS / PCoA: named clades are Leiden communities",
  "../_revision_assets/09_mds_community/candidate_labeled_mds_community.png",
  source: "09_mds_community/make_labeled_mds_community.R; full_mds.rds and arm-leiden-k15 TSV",
)

#pagebreak()

#method-slide(
  "10m",
  "Method transition",
  "How 3D contact validates sequence communities",
  [
    Freeze the sequence communities from the PHR graph.

    Ask independent contact maps whether same-community arms contact each other more than different-community arms.

    B/W ratios, Mantel rho, and exclusion tests are robustness checks on the same question.
  ],
  [The 3D block is validation of the sequence communities, not another way to call them.],
  source: "hic_methods README and hic_visual_redesign README",
)

#pagebreak()

#figure-slide(
  "10a",
  "Sequence communities co-localize in 3D",
  "../_revision_assets/hic_visual_redesign/slide_10a_square_matrix_candidate.png",
  source: "hic_visual_redesign/make_hic_visual_redesign.R; HG002 Pore-C 50 kb matrix",
)

#pagebreak()

#figure-slide(
  "10b",
  "Confound check: the 3D signal survives exclusions",
  "../_revision_assets/hic_visual_redesign/slide_10b_mantel_exclusions_clarity.png",
  source: "hic_visual_redesign/make_hic_visual_redesign.R; focused ED5 Mantel exclusion candidate",
)

#pagebreak()

#figure-slide(
  "11",
  "Single-cell 3D tests whether the bulk signal is per-cell",
  "../_revision_assets/hic_visual_redesign/slide_11_single_cell_purpose_candidate.png",
  source: "hic_visual_redesign/make_hic_visual_redesign.R; Dip-C and sperm per-cell TSVs",
)

#pagebreak()

#figure-slide(
  "12",
  "Mouse zygotene: the bouquet-stage 3D signal",
  "../_revision_assets/hic_visual_redesign/slide_12_mouse_zygotene_trajectory_pairing.png",
  source: "hic_visual_redesign/make_hic_visual_redesign.R; mouse zygotene plus stage trajectory",
)

#pagebreak()

#figure-slide(
  "13a",
  "Pedigree proof: exchange events fall in predicted communities",
  "assets/s13_pedigree_top.png",
  source: "current review-zoom top crop from s13_pedigree.png; crop recipe unknown",
)

#pagebreak()

#figure-slide(
  "13b",
  "Backup: detailed pedigree exchange events",
  "assets/s13_pedigree_bottom.png",
  source: "current review-zoom bottom crop from s13_pedigree.png; crop recipe unknown",
)

#pagebreak()

#method-slide(
  "14m",
  "Method transition",
  "How gene cargo is counted",
  [
    Intersect PHRs with subtelomeric annotation.

    Count copies and gene families across arms and communities.

    Treat DUX4/D4Z4, OR4F, and TAR1 as cargo or markers carried by the exchange network, not as causes of the network.
  ],
  [Gene cargo is a copy-aware catalog and enrichment context layered onto the PHR communities.],
  source: "14_gene_background README; 14_gene_enrichment_or4f README",
)

#pagebreak()

#figure-slide(
  "14a",
  "Gene cargo: DUX4 marks the D4Z4 PHR community",
  "assets/s14_dux4.png",
  source: "current review-zoom crop from slides/v2/slide_14_gene_biology.R; crop recipe unknown",
)

#pagebreak()

#figure-slide(
  "14b",
  "OR4F family copies mark subtelomeric exchange",
  "../_revision_assets/14_gene_enrichment_or4f/or4f_gene_family_signal.png",
  source: "14_gene_enrichment_or4f/make_or4f_gene_family_signal.R; HPRCv2 OR4F/community TSVs",
)

#pagebreak()

#figure-slide(
  "14c",
  "Backup: TAR1 as subtelomeric repeat context",
  "assets/s14_tar1.png",
  source: "current review-zoom crop from slides/v2/slide_14_gene_biology.R; crop recipe unknown",
)

#pagebreak()

#text-slide(
  "15",
  "Closing - review focus page",
  "Closing: sequence sharing, 3D proximity, exchange",
  [
    Method: implicit pangenome graph over telomere-anchored flanks.

    Empirical result: PAR-scale pseudohomology across most chromosome ends.

    Mechanism and proof: 3D proximity predicts exchange, and T2T pedigrees show it is ongoing.
  ],
  thesis: [Subtelomeres concertedly evolve through ongoing inter-chromosomal exchange, observable in pedigrees, predicted by 3D contact, and recovered by an implicit pangenome graph across HPRC v2.],
  source: "closing text updated for v2 review flow",
)
