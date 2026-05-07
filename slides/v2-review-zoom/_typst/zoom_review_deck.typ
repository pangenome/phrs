// BoG 2026 review zoom deck, v6 review revision.
// Layout: 16:9 widescreen, one visual focus per page.
// Build:
//   typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v6.pdf
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
    align(center)[#text(size: 5.6pt, fill: col-cap)[review zoom v6 focus page]]
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

#let dipc-metric-card(title, body) = block(
  width: 100%,
  fill: col-card-bg,
  stroke: 0.75pt + col-hdr-bg.darken(15%),
  inset: (x: 0.13in, y: 0.075in),
  radius: 3pt,
)[
  #text(size: 11pt, weight: "bold", fill: col-title)[#title]
  #linebreak()
  #text(size: 8.8pt, fill: col-text)[#body]
]

#let dipc-validation-panel-slide(num, label, source: "") = {
  grid(
    rows: (0.34in, 0.63in, 1fr, 0.13in),
    row-gutter: 0.035in,
    align: center,
    header(num, label),
    grid(
      columns: (1fr, 1fr),
      column-gutter: 0.10in,
      dipc-metric-card(
        "GM12878 Dip-C",
        [W/B = 0.931 (6.9% closer); Fisher p = 2.4e-05; Mantel rho = 0.296, p = 0.002.],
      ),
      dipc-metric-card(
        "Sperm scHi-C",
        [W/B = 0.401 (60% closer); Fisher p = 3.9e-51; Mantel rho = 0.202, p = 0.023.],
      ),
    ),
    box(width: 100%, height: 100%)[
      #grid(
        columns: (1fr, 1fr),
        rows: (2.76in, 2.76in),
        column-gutter: 0.08in,
        row-gutter: 0.055in,
        align: center,
        box(width: 100%, height: 2.76in)[#image("../_revision_assets/v6/dipc_validation/pdf_pngs/gm12878_mantel_scatter.png", width: 100%, height: 100%, fit: "contain")],
        box(width: 100%, height: 2.76in)[#image("../_revision_assets/v6/dipc_validation/pdf_pngs/gm12878_radial_community.png", width: 100%, height: 100%, fit: "contain")],
        box(width: 100%, height: 2.76in)[#image("../_revision_assets/v6/dipc_validation/pdf_pngs/sperm_all20_mantel_scatter.png", width: 100%, height: 100%, fit: "contain")],
        box(width: 100%, height: 2.76in)[#image("../_revision_assets/v6/dipc_validation/pdf_pngs/sperm_all20_radial_community.png", width: 100%, height: 100%, fit: "contain")],
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

#let workflow-step(num, title, body) = block(
  width: 100%,
  fill: rgb("#ffffff"),
  stroke: 0.7pt + col-hdr-bg.darken(15%),
  inset: (x: 0.14in, y: 0.095in),
  radius: 3pt,
)[
  #grid(
    columns: (0.34in, 1fr),
    column-gutter: 0.10in,
    align: top,
    block(
      fill: col-title,
      inset: (x: 0.06in, y: 0.035in),
      radius: 2pt,
    )[
      #align(center)[#text(size: 9pt, weight: "bold", fill: rgb("#ffffff"))[#num]]
    ],
    [
      #text(size: 10.5pt, weight: "bold", fill: col-title)[#title]
      #linebreak()
      #text(size: 8.9pt, fill: col-text)[#body]
    ],
  )
]

#let jaccard-workflow-slide(num, label, source: "") = {
  grid(
    rows: (0.34in, 1fr, 0.13in),
    row-gutter: 0.02in,
    align: center,
    header(num, label),
    align(center + horizon)[
      #block(width: 11.7in)[
        #align(center)[
          #text(size: 22pt, weight: "bold", fill: col-title)[How we turn PHR paths into the similarity matrix]
        ]
        #v(0.14in)
        #grid(
          columns: (1.52fr, 0.95fr),
          column-gutter: 0.18in,
          align: top,
          [
            #grid(
              rows: (auto, auto, auto, auto, auto),
              row-gutter: 0.065in,
              workflow-step(
                "1",
                "IMPG PHR interval calls",
                [Index/query HPRCv2 all-vs-all subtelomeric alignments via IMPG; keep inter-chromosomal sharing inside the terminal window.],
              ),
              workflow-step(
                "2",
                "Arm/haplotype bundles",
                [Extract one PHR sequence per sample, haplotype, and chromosome arm; all paths assigned to chr9q form the chr9q bundle.],
              ),
              workflow-step(
                "3",
                "One PGGB graph",
                [Build one graph over the full collection of PHR paths, not one graph per arm pair.],
              ),
              workflow-step(
                "4",
                "ODGI path Jaccard",
                [`odgi similarity --all -P` reports graph-node Jaccard for every path pair: shared nodes intersection / all traversed nodes union.],
              ),
              workflow-step(
                "5",
                "Arm-level matrix",
                [For every arm A x arm B, average Jaccard over all haplotype-path pairs in bundle A and bundle B.],
              ),
            )
          ],
          [
            #block(
              width: 100%,
              fill: col-card-bg,
              stroke: 0.8pt + col-hdr-bg.darken(15%),
              inset: (x: 0.18in, y: 0.15in),
              radius: 3pt,
            )[
              #text(size: 13pt, weight: "bold", fill: col-title)[Bundle average]
              #v(0.09in)
              #text(size: 10.4pt, fill: col-text)[Cell(A,B) = mean path-pair Jaccard]
              #v(0.13in)
              #text(size: 13pt, weight: "bold", fill: col-title)[Same-arm self-bundle averaging]
              #v(0.08in)
              #text(size: 10.4pt, fill: col-text)[Include same-arm A x A comparisons.]
              #v(0.08in)
              #text(size: 10.4pt, fill: col-text)[A x A can be < 1: it averages distinct haplotypes/paths from the same arm, not only each path compared to itself.]
            ]
            #v(0.16in)
            #block(
              width: 100%,
              fill: col-note-bg,
              stroke: (left: 4pt + col-note-bar),
              inset: (x: 0.17in, y: 0.12in),
              radius: 3pt,
            )[
              #text(size: 11.5pt, weight: "bold", fill: col-text)[Then use:]
              #linebreak()
              #text(size: 10.4pt, fill: col-text)[distance = 1 - Jaccard]
              #linebreak()
              #text(size: 10.4pt, fill: col-text)[arm matrix -> UPGMA tree + Leiden community views]
            ]
          ],
        )
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
  source: "narrative_density README; v5/07a_tree_then_community_heatmap and 09_mds_community READMEs",
)

#pagebreak()

#jaccard-workflow-slide(
  "07j",
  "PHR Jaccard workflow",
  source: "HPRCv2 PHR similarity workflow: IMPG interval calls, PGGB graph over all paths, ODGI path Jaccard, arm-bundle averaging",
)

#pagebreak()

#figure-slide(
  "07j.1",
  "PGGB graph main component: ODGI 2D layout",
  "../_revision_assets/v6/pggb_graph_black/pggb_graph_2d_black.png",
  source: "v6/pggb_graph_black; derived from v5/pggb_graph_odgi; ODGI layout TSV component 8 main component; 727,156 nodes; rotated 16:9; charcoal-on-white recolor",
)

#pagebreak()

#figure-slide(
  "07a.1",
  "Tree-ordered arm similarity heatmap",
  "../_revision_assets/v5/07a_tree_then_community_heatmap/07a_tree_ordered_heatmap.png",
  source: "v5/07a_tree_then_community_heatmap/make_07a_tree_then_community_heatmap.R; UPGMA leaf order drives side tree, rows, and columns",
)

#pagebreak()

#figure-slide(
  "07a.2",
  "Same matrix grouped by Leiden community",
  "../_revision_assets/v5/07a_tree_then_community_heatmap/07a_community_ordered_heatmap.png",
  source: "v5/07a_tree_then_community_heatmap; same Jaccard similarity palette/scale, no side tree, C1-C15 boxes and bands",
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
  source: "hic_methods README; v4/10a_xaxis_orientation README; v3/11_wb_labels README",
)

#pagebreak()

#figure-slide(
  "10a",
  "Sequence communities co-localize in 3D",
  "../_revision_assets/v4/10a_xaxis_orientation/candidate_10a_xaxis_orientation.png",
  source: "v4/10a_xaxis_orientation/make_10a_xaxis_orientation.R; orientation_audit.tsv; corrected X left-to-right orientation",
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

#dipc-validation-panel-slide(
  "11a",
  "Dip-C/sperm validation: Mantel and radial panels",
  source: "v6/dipc_validation/prepare_dipc_validation_assets.sh; existing rendered PDFs converted with Poppler",
)

#pagebreak()

#figure-slide(
  "11b",
  "Negative control: non-sharing S_all arms are farther apart",
  "../_revision_assets/v6/dipc_validation/plots/wb_negative_control.png",
  source: "v6/dipc_validation/make_dipc_validation_summary_plots.R; GM12878 and sperm per-cell/per-community TSVs; no 3D rerun",
)

#pagebreak()

#figure-slide(
  "11c",
  "Community-free per-cell rho: sequence similarity predicts proximity",
  "../_revision_assets/v6/dipc_validation/plots/community_free_rho_distribution.png",
  source: "v6/dipc_validation/make_dipc_validation_summary_plots.R; *_community_free_per_cell.tsv and *_community_free_arm.tsv",
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
  "Copy-number-aware enrichment: what is counted",
  [
    Standard ORA asks whether each gene name is present once; subtelomeric PHR biology repeats gene copies across arms.

    The v5 question is which families recur across gene instances, HPRCv2 arms, and Leiden communities.

    The slide section uses community-arm support and the OR4F copy gradient, not parked genome-wide copy-weighted ORA p-values.

    Caveat: community assignment is not identical to a called CHM13 PHR interval. C5 chr6_p, C7 chr13_p, C14 chrY_q, and C15 chrY_p do not have called rows in chm13.phrs.bed.
  ],
  [Say copy-aware support or copy-number-aware candidate signal; do not present exploratory weighted-ORA fold enrichments or p=0 rows as final statistics.],
  source: "v5/gene_enrichment_narrative/method_box.md and slide_text.md; v5/gene_enrichment_figures README",
)

#pagebreak()

#figure-slide(
  "14a",
  "Copy-number-aware enrichment: method boundary",
  "../_revision_assets/v5/gene_enrichment_figures/copy_aware_method_comparison.png",
  source: "v5/gene_enrichment_figures/copy_aware_method_comparison.png; method_comparison_support.tsv",
)

#pagebreak()

#figure-slide(
  "14b",
  "Copy-aware candidate signals ranked by support",
  "../_revision_assets/v5/gene_enrichment_figures/ranked_copy_aware_gene_signals.png",
  source: "v5/gene_enrichment_figures/ranked_copy_aware_gene_signals.png; ranked_signal_support.tsv; bars are support counts, not q-values",
)

#pagebreak()

#figure-slide(
  "14c",
  "Community/family map separates support from interval scope",
  "../_revision_assets/v5/gene_enrichment_figures/community_family_signal_map.png",
  source: "v5/gene_enrichment_figures/community_family_signal_map.png; community_family_map_support.tsv and community_interval_status.tsv",
)

#pagebreak()

#figure-slide(
  "14d",
  "Backup: DUX4/D4Z4 paired C1 chr4q/chr10q PHR view",
  "../_revision_assets/v3/gene_browser_panels/panel_01_dux4_d4z4_c1_chr4_chr10.png",
  source: "v3/gene_browser_panels/render_gene_browser_panels.R; panel_manifest.tsv; input_manifest.tsv",
)

#pagebreak()

#figure-slide(
  "14e",
  "Backup: OR4F-rich C3 chr3q PHR",
  "../_revision_assets/v3/gene_browser_panels/panel_02_or4f_c3_chr3q.png",
  source: "v3/gene_browser_panels/render_gene_browser_panels.R; target_loci.tsv; fixed track schema",
)

#pagebreak()

#figure-slide(
  "14f",
  "Backup: OR4F pseudogene endpoint in C8 chr15q",
  "../_revision_assets/v3/gene_browser_panels/panel_03_or4f_decay_c8_chr15q.png",
  source: "v3/gene_browser_panels/render_gene_browser_panels.R; OR4F pseudogene summary in input_manifest.tsv",
)

#pagebreak()

#figure-slide(
  "14g",
  "Backup: TAR1-rich C2 chr18p PHR",
  "../_revision_assets/v3/gene_browser_panels/panel_04_tar1_c2_chr18p.png",
  source: "v3/gene_browser_panels/render_gene_browser_panels.R; TAR1 repeat lane kept separate from gene models",
)

#pagebreak()

#figure-slide(
  "14h",
  "Backup: C7 acrocentric p-arm panel uses one track grammar",
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
