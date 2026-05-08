// BoG 2026 review zoom deck, v9 review revision.
// Layout: 16:9 widescreen, one visual focus per page.
// Build:
//   typst compile --root .. zoom_review_deck.typ ../BoG_2026_review_zoom_v9.pdf
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
    align(center)[#text(size: 5.6pt, fill: col-cap)[review zoom v8 focus page]]
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

#let raster-slide(path) = {
  box(width: 100%, height: 100%)[
    #image(path, width: 100%, height: 100%, fit: "contain")
  ]
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

#let community-method-stat(value, label) = block(
  width: 100%,
  fill: col-card-bg,
  stroke: 0.75pt + col-hdr-bg.darken(15%),
  inset: (x: 0.08in, y: 0.06in),
  radius: 3pt,
)[
  #align(center)[
    #text(size: 21.2pt, weight: "bold", fill: col-title)[#value]
    #linebreak()
    #text(size: 9.8pt, fill: col-text)[#label]
  ]
]

#let community-method-card(title, body) = block(
  width: 100%,
  fill: col-card-bg,
  stroke: 0.75pt + col-hdr-bg.darken(15%),
  inset: (x: 0.12in, y: 0.07in),
  radius: 3pt,
)[
  #text(size: 12.8pt, weight: "bold", fill: col-title)[#title]
  #linebreak()
  #text(size: 10.9pt, fill: col-text)[#body]
]

#let community-assignment-method-slide(num, label, source: "") = {
  grid(
    rows: (0.34in, 1fr, 0.13in),
    row-gutter: 0.02in,
    align: center,
    header(num, label),
    align(center + horizon)[
      #block(width: 11.95in)[
        #align(center)[
          #text(size: 27.5pt, weight: "bold", fill: col-title)[How we assigned PHR communities]
        ]
        #v(0.08in)
        #grid(
          columns: (1fr, 1fr, 1fr, 1fr),
          column-gutter: 0.08in,
          community-method-stat("15,668", [PHR paths with inter-chromosomal signal]),
          community-method-stat("41 x 41", [arm matrix; 7 zero-signal arms excluded]),
          community-method-stat("15", [Leiden arm-level communities, C1-C15]),
          community-method-stat("1.16 / 0.347", [resolution / mean silhouette]),
        )
        #v(0.08in)
        #grid(
          columns: (1.13fr, 1fr),
          column-gutter: 0.16in,
          align: top,
          [
            #block(
              width: 100%,
              height: 3.22in,
              fill: rgb("#ffffff"),
              stroke: 0.75pt + col-hdr-bg.darken(15%),
              inset: (x: 0.05in, y: 0.05in),
              radius: 3pt,
            )[
              #image("../_revision_assets/v9/labels_superpop_stats_polish/community_assignment_method_schematic_v9_readable.svg", width: 100%, height: 100%, fit: "contain")
            ]
            #v(0.07in)
            #block(
              width: 100%,
              fill: col-note-bg,
              stroke: (left: 4pt + col-note-bar),
              inset: (x: 0.13in, y: 0.075in),
              radius: 3pt,
            )[
              #text(size: 11.5pt, weight: "bold", fill: col-text)[Community definitions use graph similarity only.]
              #linebreak()
              #text(size: 10.8pt, fill: col-text)[Biological names and 3D contact evidence are annotations interpreted after clustering.]
            ]
          ],
          [
            #grid(
              rows: (auto, auto, auto),
              row-gutter: 0.06in,
              community-method-card(
                "Inputs: graph-path Jaccard",
                [Build one PGGB graph (`pggb -p 95`), then compute path Jaccard with `odgi similarity --all -P`.],
              ),
              community-method-card(
                "Collapse paths to arms",
                [For each arm pair, average haplotype/path pair distances using `distance = 1 - Jaccard`.],
              ),
              community-method-card(
                "Silhouette-selected Leiden",
                [Run Leiden with `w_ij = exp(-d_ij / median(d))`; scan resolution 0.1-3.0 and choose the maximum mean silhouette.],
              ),
            )
            #v(0.06in)
            #block(
              width: 100%,
              fill: rgb("#ffffff"),
              stroke: 0.75pt + col-hdr-bg.darken(12%),
              inset: (x: 0.12in, y: 0.075in),
              radius: 3pt,
            )[
              #text(size: 12.8pt, weight: "bold", fill: col-title)[Robustness check]
              #linebreak()
              #text(size: 10.8pt, fill: col-text)[UPGMA average-linkage on the same matrix gives 14 communities, silhouette 0.342, and exact agreement on 12/15 Leiden communities.]
            ]
          ],
        )
        #v(0.08in)
        #align(center)[
          #text(size: 10.4pt, fill: col-cap)[Arm-level C1-C15 partition across 41 detected-signal arms; the sequence-level 50-community partition is separate.]
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
  "PHR length distribution by chromosome end",
  "../_revision_assets/v9/06a_q_axis_kbp/phr_length_arm_heatstrip_10kbp.png",
  source: "v9/06a_q_axis_kbp/make_06a_q_axis_kbp.R; arm_length_bins_10kbp.tsv; >500 kbp not measured",
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

#captioned-figure-slide(
  "07j.a",
  "Chromosome-end bundle",
  "A chromosome-end bundle is a population collection",
  "../_revision_assets/v9/broom_jaccard_method/broom_1.png",
  [All HPRCv2 haplotype PHR paths assigned to one chromosome end are treated as one subtelomere bundle. The unit is the population collection for that end, not a single reference interval.],
  source: "Erik Garrison broom cartoon copied from repo-root broom_1.png; method wording: subtelomeric_analysis_report.md sections 5 and 6.1",
)

#pagebreak()

#captioned-figure-slide(
  "07j.b",
  "Bundle Jaccard",
  "Compare two chromosome-end bundles by graph-node overlap",
  "../_revision_assets/v9/broom_jaccard_method/brooms_compare.png",
  [For paths drawn from two chromosome-end bundles, Jaccard is intersection / union of the PGGB variation graph nodes they traverse. Averaged path-pair overlaps form the arm-level matrix for heatmaps, Leiden, and MDS; A x A can be below 1 when a bundle has heterogeneous haplotypes/paths.],
  source: "Erik Garrison broom cartoon copied from repo-root brooms_compare.png; method wording: subtelomeric_analysis_report.md sections 5 and 6.1 and current slide 07j",
)

#pagebreak()

#figure-slide(
  "07j.1",
  "PGGB graph main component: ODGI 2D layout",
  "../_revision_assets/v6/pggb_graph_black/pggb_graph_2d_black.png",
  source: "v6/pggb_graph_black; derived from v5/pggb_graph_odgi; ODGI layout TSV component 8 main component; 727,156 nodes; rotated 16:9; charcoal-on-white recolor",
)

#pagebreak()

#community-assignment-method-slide(
  "07j.2",
  "Community assignment method",
  source: "v9/labels_superpop_stats_polish/slide07j2_typst_patch.typ; subtelomeric_analysis_report.md sections 5 and 6.1; HPRCv2 plot-similarity-subtelo.R; arm distance matrix and Leiden k15 assignment TSVs",
)

#pagebreak()

#figure-slide(
  "07a.1",
  "Tree-ordered arm similarity heatmap",
  "../_revision_assets/v5/07a_tree_then_community_heatmap/07a_tree_ordered_heatmap.png",
  source: "v5/07a_tree_then_community_heatmap/make_07a_tree_then_community_heatmap.R; UPGMA leaf order drives side tree, rows, and columns",
)

#pagebreak()

#captioned-figure-slide(
  "07a.1b",
  "Why Leiden for community assignment?",
  "Leiden adds refinement before aggregation",
  "../_revision_assets/v9/leiden_figure_slide/leiden_algorithm_fig3_official.png",
  [Leiden adds a refinement step so communities remain well connected; we use it on the arm-level graph-path Jaccard similarity network.],
  source: "Adapted from Traag, V. A., Waltman, L. & van Eck, N. J. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports 9, 5233 (2019), CC BY 4.0.",
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
  "How to read MDS and population-neighbor spread",
  [
    MDS places sequence-level PHR flanks by graph-path Jaccard distance.

    Population spread is summarized as each point's nearest other same-superpopulation neighbor in displayed D1-D2 MDS space.

    The unit is a subtelomeric flank, so arm-community structure dominates the geometry.
  ],
  [The v8 distance plot is nearest same-superpopulation neighbor only: self is excluded, and it is not centroid or all-pairwise distance.],
  source: "v8/mds_superpop_community_polish/SLIDE_PATCH.md; nearest_same_superpop_mds_summary.tsv; community_mds_label_positions.tsv",
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
  "MDS colored by continental superpopulation",
  "../_revision_assets/v8/mds_superpop_community_polish/superpopulation_mds_08a_matched.png",
  source: "v8/mds_superpop_community_polish/make_mds_superpop_community_polish.R; MDS D1-D2 from hprcv2.1Mb.subtelo.full_mds.rds; 08a-matched 12x10 MDS render style; color is continental superpopulation",
)

#pagebreak()

#captioned-figure-slide(
  "08b.1",
  "Nearest same-superpopulation neighbor in MDS space",
  "Nearest same-superpopulation neighbor in MDS space",
  "../_revision_assets/v9/labels_superpop_stats_polish/nearest_same_superpop_distance_boxplot_bracketed.png",
  [
    For each subtelomeric MDS point, distance is to the nearest other point from the same continental superpopulation in displayed D1-D2 MDS space. Self is excluded. Boxes show distributions; printed values are means. KW is the Kruskal-Wallis global non-parametric test across groups; pairwise Wilcoxon is rank-sum group comparison; BH is Benjamini-Hochberg FDR correction over pairwise tests. Brackets show the five strongest BH-significant contrasts.
  ],
  source: "v9/labels_superpop_stats_polish/nearest_same_superpop_mds_distances.tsv; bracketed boxplot summary; Wilcoxon BH and Cliff delta TSVs",
)

#pagebreak()

#figure-slide(
  "09",
  "MDS: all Leiden communities C1-C15 labeled",
  "../_revision_assets/v8/mds_superpop_community_polish/community_mds_labeled.png",
  source: "v8/mds_superpop_community_polish/make_mds_superpop_community_polish.R; MDS from hprcv2.1Mb.subtelo.full_mds.rds; C1-C15 Leiden labels from community_mds_label_positions.tsv; 1:1 axes; not PCA",
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
  "10m.1",
  "Making Hi-C work at subtelomeric repeats",
  "../_revision_assets/v7/hic_mapq0_methods/hic_mapq0_methods_flow.svg",
  source: "v7/hic_mapq0_methods/README.md; reports 05:5-11, 06:5-7, 07:34-42, 10:35-37; MAPQ0 aggregate signal caveat",
)

#pagebreak()

#figure-slide(
  "10m.2",
  "CHM13 Hi-C contact-space MDS in 3D",
  "../_revision_assets/v9/slide10m2_better_3d_viz/chm13_phr_contact_mds_3d_view.png",
  source: "v9/slide10m2_better_3d_viz/chm13_phr_contact_mds_3d_view.png; same D1-D2-D3 contact-space MDS as next slide; colors are sequence communities; not physical 3DG",
)

#pagebreak()

#figure-slide(
  "10m.3",
  "The same 3D MDS is clearer as 2D projections",
  "../_revision_assets/v9/slide10m2_better_3d_viz/best_replacement_chm13_phr_contact_mds.png",
  source: "v9/slide10m2_better_3d_viz; source matrix: chm13_hic.dist_matrix.tsv; regions: chm13_subtelomeric_regions.bed; color: sequence communities; contact-space MDS, not physical 3DG",
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

#captioned-figure-slide(
  "11a.1",
  "GM12878 Dip-C proximity",
  "GM12878 Dip-C: sequence similarity predicts 3D proximity",
  "../_revision_assets/v8/hic_dipc_clarity_split/plots/gm12878_mantel_proximity.png",
  [
    W/B means within-community divided by between-community 3D distance; W/B < 1 means within-community arms are closer. Proximity convention: `3D proximity = -mean 3D distance`, so higher y-values are closer and positive rho has the expected direction.
  ],
  source: "v8/hic_dipc_clarity_split/Source_manifest.tsv; Dataset GM12878; Technology Dip-C; n=16 cells; source TSVs gm12878_mantel_3d.tsv and arm distance matrices",
)

#pagebreak()

#captioned-figure-slide(
  "11a.2",
  "GM12878 Dip-C radial",
  "GM12878 Dip-C radial community structure",
  "../_revision_assets/v8/hic_dipc_clarity_split/inputs/gm12878_radial_community.png",
  [
    GM12878 only: per-community normalized radial position plus within-vs-between radial-position similarity. Same-community arms have more similar radial positions than between-community arms.
  ],
  source: "v8/hic_dipc_clarity_split/Source_manifest.tsv; Dataset GM12878; Technology Dip-C; n=16 cells; gm12878_radial_community.pdf/TSV",
)

#pagebreak()

#captioned-figure-slide(
  "11a.3",
  "Sperm scHi-C proximity",
  "Sperm single-cell Hi-C: sequence similarity predicts 3D proximity",
  "../_revision_assets/v8/hic_dipc_clarity_split/plots/sperm_all20_mantel_proximity.png",
  [
    W/B means within-community divided by between-community 3D distance; W/B < 1 means within-community arms are closer. Proximity convention: `3D proximity = -mean 3D distance`, so higher y-values are closer and positive rho has the expected direction.
  ],
  source: "v8/hic_dipc_clarity_split/Source_manifest.tsv; Dataset human sperm; Technology single-cell Hi-C/3DG; n=20 cells; sperm_all20_mantel_3d.tsv and arm distance matrices",
)

#pagebreak()

#captioned-figure-slide(
  "11a.4",
  "Sperm scHi-C radial",
  "Sperm single-cell Hi-C radial community structure",
  "../_revision_assets/v8/hic_dipc_clarity_split/inputs/sperm_all20_radial_community.png",
  [
    Sperm only: per-community normalized radial position plus within-vs-between radial-position similarity across 20 haploid sperm cells. It measures the same radial statistic as the GM12878 radial panel in a distinct cell type.
  ],
  source: "v8/hic_dipc_clarity_split/Source_manifest.tsv; Dataset human sperm; Technology single-cell Hi-C/3DG; n=20 cells; sperm_all20_radial_community.pdf/TSV",
)

#pagebreak()

#captioned-figure-slide(
  "11b",
  "Negative control: non-sharing S_all arms are farther apart",
  "Negative control: non-sharing S_all arms are farther apart",
  "../_revision_assets/v8/hic_dipc_clarity_split/plots/wb_negative_control_reduced_text.png",
  [
    Reduced text scale. W/B is within-community divided by between-community 3D distance; values below 1 mean closer within-community, while S_all values above 1 show the non-sharing negative control is farther apart.
  ],
  source: "v8/hic_dipc_clarity_split/Source_manifest.tsv; per-cell and per-community TSV roots for GM12878 Dip-C and sperm single-cell Hi-C",
)

#pagebreak()

#captioned-figure-slide(
  "11c",
  "Community-free per-cell rho: sequence similarity predicts proximity",
  "Community-free per-cell rho: sequence similarity predicts proximity",
  "../_revision_assets/v8/hic_dipc_clarity_split/plots/community_free_rho_distribution_reduced_text.png",
  [
    Reduced text scale. Rho is computed as sequence similarity versus 3D proximity (`-distance`), so positive rho means more similar sequences are closer in 3D.
  ],
  source: "v8/hic_dipc_clarity_split/Source_manifest.tsv; files *_community_free_per_cell.tsv and *_community_free_arm.tsv for GM12878 Dip-C and sperm single-cell Hi-C",
)

#pagebreak()

#figure-slide(
  "12",
  "Mouse zygotene: the bouquet-stage 3D signal",
  "../_revision_assets/v8/typography_legend_cleanup/slide12_mouse_zygotene_large_text.png",
  source: "v8/typography_legend_cleanup/make_typography_legend_cleanup.R; mouse Zuo 2021 stage tables; zygotene pair correlation and stage Mantel rho",
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
  "../_revision_assets/v8/typography_legend_cleanup/slide13b_pedigree_bottom_no_unused_legend.png",
  source: "v8/typography_legend_cleanup/crop_png_top.py; crop of s13_pedigree_bottom.png with unused bottom legend removed; event labels are direct in-panel annotations",
)

#pagebreak()

#method-slide(
  "14m",
  "Method transition",
  "Gene enrichment: report-backed biology, conservative statistics",
  [
    The canonical HPRCv2 community-family Fisher screen tested 116 rows and has 0 BH-significant rows.

    C3 OR and C7 MTCO are near-miss copy-aware candidate signals with BH q = 0.07118, not definitive enriched classes.

    The biology is still report-backed: recurrent OR4F architecture, C1 D4Z4/DUX4L, C5 DDX11L/WASH/FAM138/IQSEC3, C15 PAR1 coding genes, C7 MTCO pseudogenes, and a pseudogene/ncRNA-rich backbone.

    GO/ORA tables are useful context or method contrast only: standard/coding-only outputs collapse copy number, and exploratory copy-weighted ORA still needs calibrated validation.
  ],
  [Say recurrent gene-family architecture, report-backed community-arm support, or copy-aware candidate signal; do not present GO or Fisher rows as BH-significant proof.],
  source: "subtelomeric_analysis_report.md:477-600; end-to-end-report/report/03_gene_enrichment.md; v7/gene_enrichment_report_backed README",
)

#pagebreak()

#figure-slide(
  "14a",
  "Report-backed gene-family architecture",
  "../_revision_assets/v7/gene_enrichment_report_backed/gene_enrichment_report_backed_summary.svg",
  source: "v7/gene_enrichment_report_backed/gene_enrichment_report_backed_summary.svg; report section 9; canonical Fisher screen has 116 rows and 0 BH-significant rows",
)

#pagebreak()

#figure-slide(
  "14b",
  "Copy-aware candidate signals ranked by community-arm support",
  "../_revision_assets/v8/typography_legend_cleanup/slide14b_candidate_signals_talk_ready.png",
  source: "v8/typography_legend_cleanup/make_typography_legend_cleanup.R; v5 ranked_signal_support.tsv; bars are support counts, not q-values or BH-significant effects",
)

#pagebreak()

#figure-slide(
  "14c",
  "Community/family map separates support from statistical proof",
  "../_revision_assets/v8/typography_legend_cleanup/slide14c_community_family_map_talk_ready.png",
  source: "v8/typography_legend_cleanup/make_typography_legend_cleanup.R; v5 community_family_map_support.tsv simplified for talk legibility; direct labels replace legend",
)

#pagebreak()

#raster-slide("../_revision_assets/v8/chm13_ucsc_examples/selected_example-01.png")

#pagebreak()

#raster-slide("../_revision_assets/v8/chm13_ucsc_examples/selected_example-02.png")

#pagebreak()

#raster-slide("../_revision_assets/v8/chm13_ucsc_examples/selected_example-03.png")

#pagebreak()

#raster-slide("../_revision_assets/v8/chm13_ucsc_examples/selected_example-04.png")

#pagebreak()

#raster-slide("../_revision_assets/v8/chm13_ucsc_examples/selected_example-05.png")

#pagebreak()

#raster-slide("../_revision_assets/v8/chm13_ucsc_examples/selected_example-06.png")
