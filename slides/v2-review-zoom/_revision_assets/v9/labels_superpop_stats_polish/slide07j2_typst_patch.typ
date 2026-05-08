// Drop-in replacement guidance for slide 07j.2 in
// slides/v2-review-zoom/_typst/zoom_review_deck.typ.
//
// Replace the existing `community-method-stat`, `community-method-card`, and
// `community-assignment-method-slide` definitions with the block below. Font
// sizes are about 25% larger than the v8/v6 deck macro, assignment wording is
// direct, and the source/provenance footer remains explicit.

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
