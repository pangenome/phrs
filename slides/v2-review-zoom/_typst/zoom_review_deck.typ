// BoG 2026 review zoom deck
// Layout: 16:9 widescreen, one visual focus per page.
// Build:
//   typst compile zoom_review_deck.typ ../BoG_2026_review_zoom.pdf
//   typst compile --ppi 144 zoom_review_deck.typ page-{0p}.png

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
    align(center)[#text(size: 5.6pt, fill: col-cap)[review zoom focus page]]
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

#let text-slide(num, label, title, body, thesis: none) = {
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
    footer("text-only focus page"),
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

#figure-slide(
  "02",
  "Implicit interval tree",
  "assets/s02_interval_tree.png",
  source: "canonical review asset: s02_interval_tree.png",
)

#pagebreak()

#figure-slide(
  "03a",
  "IMPG workflow",
  "assets/s03_impg_workflow.png",
  source: "canonical review asset: s03_impg_workflow.png",
)

#pagebreak()

#figure-slide(
  "03b",
  "Erdos-Renyi callout",
  "assets/s03_er_callout.png",
  source: "canonical review asset: s03_er_callout.png",
)

#pagebreak()

#figure-slide(
  "04",
  "Fig 1 panel a - genome-wide heatmap",
  "assets/s04_fig1_panel_a.png",
  source: "focused panel crop copied into this deck from local zoom assets; source figure s04_fig1.png",
)

#pagebreak()

#figure-slide(
  "05",
  "Interchromosomal similarities",
  "assets/s05_interchrom.png",
  source: "canonical review asset: s05_interchrom.png",
)

#pagebreak()

#quad-slide(
  "06a",
  "Length distribution - full faceted split view",
  "assets/s06_length_dist_1_top.png",
  "assets/s06_length_dist_2_top.png",
  "assets/s06_length_dist_1_bottom.png",
  "assets/s06_length_dist_2_bottom.png",
  source: "four focused crops copied into this deck from local zoom assets; source figure s06_length_dist.png",
)

#pagebreak()

#figure-slide(
  "06b",
  "Length distribution clade callouts",
  "assets/s06_clade_callouts.png",
  source: "canonical review asset: s06_clade_callouts.png",
)

#pagebreak()

#figure-slide(
  "07a",
  "Fig 1 panel c - arm heatmap",
  "assets/s07_fig1_panel_c.png",
  source: "focused panel crop copied into this deck from local zoom assets; source figure s04_fig1.png",
)

#pagebreak()

#figure-slide(
  "07b",
  "NJ tree with named clades",
  "assets/s07b_nj_tree.png",
  source: "canonical review asset: s07b_nj_tree.png",
)

#pagebreak()

#figure-slide(
  "08a",
  "MDS colored by chromosome",
  "assets/s08a_mds_chrom.png",
  source: "canonical review asset: s08a_mds_chrom.png",
)

#pagebreak()

#figure-slide(
  "08b",
  "MDS colored by superpopulation",
  "assets/s08b_mds_superpop.png",
  source: "canonical review asset: s08b_mds_superpop.png",
)

#pagebreak()

#figure-slide(
  "09a",
  "PCA communities",
  "assets/s09_pca_communities.png",
  source: "canonical review asset: s09_pca_communities.png",
)

#pagebreak()

#figure-slide(
  "09b",
  "Clade legend",
  "assets/s09_clade_legend.png",
  source: "canonical review asset: s09_clade_legend.png",
)

#pagebreak()

#figure-slide(
  "09c",
  "Community-colored MDS",
  "assets/s09b_communities.png",
  source: "canonical review asset: s09b_communities.png",
)

#pagebreak()

#figure-slide(
  "10a",
  "Fig 3 panel a - Hi-C/Pore-C contact matrix",
  "assets/s10_fig3_panel_a.png",
  source: "focused panel crop copied into this deck from local zoom assets; source figure s10_fig3.png",
)

#pagebreak()

#figure-slide(
  "10b",
  "Mantel exclusions",
  "assets/s10b_ed5.png",
  source: "canonical review asset: s10b_ed5.png",
)

#pagebreak()

#figure-slide(
  "11",
  "Fig 3 panel c - single-cell 3D",
  "assets/s10_fig3_panel_c.png",
  source: "focused panel crop copied into this deck from local zoom assets; source figure s10_fig3.png",
)

#pagebreak()

#figure-slide(
  "12a",
  "Fig 4 panel d - mouse meiotic Hi-C",
  "assets/s12_fig4_panel_d.png",
  source: "focused panel crop copied into this deck from local zoom assets; source figure s12_fig4.png",
)

#pagebreak()

#figure-slide(
  "12b",
  "Meiotic stage trajectory inset",
  "assets/s12_trajectory.png",
  source: "canonical review asset: s12_trajectory.png",
)

#pagebreak()

#figure-slide(
  "13a",
  "Pedigree exchange figure - top half",
  "assets/s13_pedigree_top.png",
  source: "top-half readability crop copied into this deck from local zoom assets; source figure s13_pedigree.png",
)

#pagebreak()

#figure-slide(
  "13b",
  "Pedigree exchange figure - bottom half",
  "assets/s13_pedigree_bottom.png",
  source: "bottom-half readability crop copied into this deck from local zoom assets; source figure s13_pedigree.png",
)

#pagebreak()

#figure-slide(
  "14a",
  "Gene biology - DUX4/D4Z4",
  "assets/s14_dux4.png",
  source: "focused panel crop copied into this deck from local zoom assets; source figure s14_gene_biology.png",
)

#pagebreak()

#figure-slide(
  "14b",
  "Gene biology - OR4F",
  "assets/s14_or4f.png",
  source: "focused panel crop copied into this deck from local zoom assets; source figure s14_gene_biology.png",
)

#pagebreak()

#figure-slide(
  "14c",
  "Gene biology - TAR1",
  "assets/s14_tar1.png",
  source: "focused panel crop copied into this deck from local zoom assets; source figure s14_gene_biology.png",
)

#pagebreak()

#text-slide(
  "15",
  "Closing - review focus page",
  "Concerted evolution of human subtelomeres",
  [
    Method: implicit pangenome graph over telomere-anchored flanks.

    Empirical result: PAR-scale pseudohomology across most chromosome ends.

    Mechanism and proof: 3D proximity predicts exchange, and T2T pedigrees show it is ongoing.
  ],
  thesis: [Subtelomeres concertedly evolve through ongoing inter-chromosomal exchange - observable in pedigrees, predicted by 3D, recovered by an implicit pangenome graph across HPRC v2.],
)
