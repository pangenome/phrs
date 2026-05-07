// BoG 2026 review zoom deck, v3 review revision.
// Layout: 16:9 widescreen, one visual focus per page.
// Build:
//   typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v3.pdf
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
    align(center)[#text(size: 5.6pt, fill: col-cap)[review zoom v3 focus page]]
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
  "PHR length distributions by named clade, censored at 500 kb",
  "../_revision_assets/v3/06_violin_censor/candidate_06a_named_clade_violin_censor.png",
  source: "v3/06_violin_censor/make_06_violin_censor.R; named_clade_violin_summary.tsv; 500 kb analysis cap",
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
  "../_revision_assets/v3/07a_crisp_aligned/candidate_07a_upgma_crisp_aligned.png",
  source: "v3/07a_crisp_aligned/make_07a_crisp_aligned.R; order_validation.tsv confirms tree tips equal heatmap rows",
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
  "Backup: unrooted NJ with acrocentric p-arm audit",
  "../_revision_assets/v3/07c_acrocentric_presence/07c_unrooted_acrocentric_status.png",
  source: "v3/07c_acrocentric_presence/make_07c_acrocentric_presence.R; all five acrocentric p arms present in C7",
)

#pagebreak()

#method-slide(
  "08m",
  "Method transition",
  "How to read MDS / PCoA and population spread",
  [
    MDS / PCoA places sequence-level PHR flanks by Jaccard distance.

    Population variation is summarized with same-superpopulation pairwise distances in the displayed MDS / PCoA space, not a distance-to-center summary.

    The unit is a subtelomeric flank, so arm-community structure dominates the geometry.
  ],
  [Use MDS / PCoA language consistently and keep population spread as a secondary within-population variation check.],
  source: "v3/08b_within_pop_pairwise README; v3/09_all_communities_1to1 README",
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
  "Within-population pairwise variation on the PHR MDS / PCoA",
  "../_revision_assets/v3/08b_within_pop_pairwise/within_pop_pairwise_2d_distribution.png",
  source: "v3/08b_within_pop_pairwise/make_within_pop_pairwise.R; within_pop_pairwise_summary.tsv",
)

#pagebreak()

#figure-slide(
  "09",
  "MDS / PCoA: all Leiden communities C1-C15 labeled",
  "../_revision_assets/v3/09_all_communities_1to1/mds_pcoa_all_communities_1to1.png",
  source: "v3/09_all_communities_1to1/make_all_communities_1to1.R; validation_summary.tsv confirms all C1-C15 labels",
)

#pagebreak()

#method-slide(
  "10m",
  "Method transition",
  "How 3D contact validates sequence communities",
  [
    Freeze the sequence communities from the PHR graph.

    Ask independent contact maps whether same-community arms contact each other more than different-community arms.

    Within-community versus between-community distance ratios, Mantel rho, and exclusion tests are robustness checks on the same question.
  ],
  [The 3D block is validation of the sequence communities, not another way to call them.],
  source: "hic_methods README; v3/10a_axis_box_fix README; v3/11_wb_labels README",
)

#pagebreak()

#figure-slide(
  "10a",
  "Sequence communities co-localize in 3D",
  "../_revision_assets/v3/10a_axis_box_fix/candidate_10a_axis_box_fix.png",
  source: "v3/10a_axis_box_fix/make_10a_axis_box_fix.R; matrix_order_audit.tsv; community boxes aligned to square matrix",
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
  "Single-cell 3D: within-community arms are closer than between-community arms",
  "../_revision_assets/v3/11_wb_labels/slide11_explicit_distance_labels_candidate.png",
  source: "v3/11_wb_labels/make_slide11_explicit_distance_labels.R; explicit within-community vs between-community distance labels",
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
  "12b",
  "Human arms: sequence-similar subtelomeres contact more",
  "../_revision_assets/v3/human_3d_dotplot/human_arm_pair_dotplot_candidate.png",
  source: "v3/human_3d_dotplot/make_human_3d_dotplot.R; human_3d_dotplot_summary.tsv; arm-pair pointwise Spearman",
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
  source: "v3/gene_browser_panels README; v3/gene_browser_inventory README; fixed track grammar",
)

#pagebreak()

#figure-slide(
  "14a",
  "Gene cargo: DUX4/D4Z4 is a paired C1 chr4q/chr10q PHR view",
  "../_revision_assets/v3/gene_browser_panels/panel_01_dux4_d4z4_c1_chr4_chr10.png",
  source: "v3/gene_browser_panels/render_gene_browser_panels.R; panel_manifest.tsv; input_manifest.tsv",
)

#pagebreak()

#figure-slide(
  "14b",
  "Gene cargo: OR4F-rich C3 chr3q PHR",
  "../_revision_assets/v3/gene_browser_panels/panel_02_or4f_c3_chr3q.png",
  source: "v3/gene_browser_panels/render_gene_browser_panels.R; target_loci.tsv; fixed track schema",
)

#pagebreak()

#figure-slide(
  "14c",
  "Gene cargo: OR4F pseudogene endpoint in C8 chr15q",
  "../_revision_assets/v3/gene_browser_panels/panel_03_or4f_decay_c8_chr15q.png",
  source: "v3/gene_browser_panels/render_gene_browser_panels.R; OR4F pseudogene summary in input_manifest.tsv",
)

#pagebreak()

#figure-slide(
  "14d",
  "Gene cargo: TAR1-rich C2 chr18p PHR",
  "../_revision_assets/v3/gene_browser_panels/panel_04_tar1_c2_chr18p.png",
  source: "v3/gene_browser_panels/render_gene_browser_panels.R; TAR1 repeat lane kept separate from gene models",
)

#pagebreak()

#figure-slide(
  "14e",
  "Gene cargo: C7 acrocentric p-arm panel uses one track grammar",
  "../_revision_assets/v3/gene_browser_panels/panel_05_acrocentric_c7_p_arm_group.png",
  source: "v3/gene_browser_panels/render_gene_browser_panels.R; C7 PanSN relative offsets from all-vs-all length TSV",
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
  source: "closing text retained for v3 review flow",
)
