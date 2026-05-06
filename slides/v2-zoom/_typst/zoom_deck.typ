#set page(width: 297mm, height: 210mm, margin: 4mm)
#set text(font: "DejaVu Sans", fill: rgb("#111111"))

#let title_style = text.with(size: 18pt, weight: "semibold")
#let footer_style = text.with(size: 8pt, fill: rgb("#777777"))

#let figure_slide(num, title, path, body) = {
  pagebreak(weak: true)
  grid(
    rows: (7%, 89%, 4%),
    row-gutter: 0mm,
    align: center,
    [
      #box(width: 100%, inset: (x: 2mm, y: 0mm))[
        #title_style[#num #title]
      ]
    ],
    [
      #box(width: 100%, height: 100%)[
        #body
      ]
    ],
    [
      #box(width: 100%, inset: (x: 2mm, y: 0mm))[
        #footer_style[#path]
      ]
    ],
  )
}

#let image_slide(num, title, path) = figure_slide(num, title, path, image(path, width: 99%, height: 100%, fit: "contain"))

#let text_slide(num, title, excerpt) = {
  pagebreak(weak: true)
  grid(
    rows: (1fr, auto, 1fr, 7%),
    align: center,
    [],
    [
      #box(width: 82%)[
        #align(center)[
          #text(size: 31pt, weight: "bold")[#title]

          #v(7mm)
          #text(size: 17pt, fill: rgb("#333333"))[#excerpt]
        ]
      ]
    ],
    [],
    [
      #footer_style[[no figure - text-only slide]]
    ],
  )
}

#text_slide(
  "01",
  "Concerted evolution and unorthodox recombination of human subtelomeres",
  [Andrea Guarracino, Erik Garrison. Companion to HPRC v2 (Nature, in submission). Inter-chromosomal subtelomeric relationships at HPRC v2 scale across 466 near-complete haplotypes. Biology of Genomes, Cold Spring Harbor - May 2026.]
)

#image_slide(
  "02",
  "Implicit interval tree - the data structure under the implicit pangenome graph",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2/_typst/v1_page_02-02.png",
)

#image_slide(
  "03a",
  "The implicit pangenome graph - v1 workflow stack",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2/_typst/v1_page_03-03.png",
)

#image_slide(
  "03b",
  "Erdos-Renyi sampling threshold callout",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2/slide_03_er_callout.png",
)

#image_slide(
  "04",
  "Genome-wide identity heatmap - interchromosomal homology at PAR2 scale",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_04_fig1_panel_a.png",
)

#image_slide(
  "05",
  "Interchromosomal similarities - n-chromosomes per region",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/p_num_chromosomes_wide.png",
)

#image_slide(
  "06a",
  "Length distributions of inter-chromosomal matches per arm - panel 1 top",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_06_length_dist_1_top.png",
)

#image_slide(
  "06b",
  "Length distributions of inter-chromosomal matches per arm - panel 1 bottom",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_06_length_dist_1_bottom.png",
)

#image_slide(
  "06c",
  "Length distributions of inter-chromosomal matches per arm - panel 2 top",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_06_length_dist_2_top.png",
)

#image_slide(
  "06d",
  "Length distributions of inter-chromosomal matches per arm - panel 2 bottom",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_06_length_dist_2_bottom.png",
)

#image_slide(
  "06e",
  "Length-distribution clade callouts",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2/slide_06_clade_callouts.png",
)

#image_slide(
  "07a",
  "All-vs-all at the arm level - clustered heatmap",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_07_fig1_panel_c.png",
)

#image_slide(
  "07b",
  "All-vs-all at the arm level - annotated NJ tree",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/paper_prep/figures/nj_tree_arms/nj_tree_annotated.png",
)

#image_slide(
  "08a",
  "All-vs-all in 2D - colored by chromosome",
  "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-chromosome.png",
)

#image_slide(
  "08b",
  "All-vs-all in 2D - colored by superpopulation",
  "/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-superpopulation.png",
)

#image_slide(
  "09a",
  "All-vs-all PCA - 15 arm-level communities",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2/_typst/v1_page_10-10.png",
)

#image_slide(
  "09b",
  "Abstract-anchored clade legend",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2/slide_09_clade_legend.png",
)

#image_slide(
  "10a",
  "Hi-C / Pore-C confirm sequence communities are 3D",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_10_fig3_panel_a.png",
)

#image_slide(
  "10b",
  "Mantel exclusions - signal strengthens after confounds are removed",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/paper_prep/figures/ed5/figure_ed5.png",
)

#image_slide(
  "11",
  "Single-cell 3D - and it works in haploid sperm",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_11_fig3_panel_c.png",
)

#image_slide(
  "12a",
  "Mouse meiotic Hi-C - zygotene per-PHR-pair scatter",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_12_fig4_panel_d.png",
)

#image_slide(
  "12b",
  "Mouse meiotic Hi-C - stage trajectory inset",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2/slide_12_stage_trajectory.png",
)

#image_slide(
  "13a",
  "Caught in the act - PAN027 inherited from PAN010, upper half",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_13a_pedigree_top.png",
)

#image_slide(
  "13b",
  "Caught in the act - PAN027 inherited from PAN010, lower half",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_13b_pedigree_bottom.png",
)

#image_slide(
  "14a",
  "Gene biology aside - DUX4 / D4Z4",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_14a_dux4.png",
)

#image_slide(
  "14b",
  "Gene biology aside - OR4F",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_14b_or4f.png",
)

#image_slide(
  "14c",
  "Gene biology aside - TAR1",
  "/moosefs/erikg/phrs/.wg-worktrees/agent-878/slides/v2-zoom/_typst/assets/slide_14c_tar1.png",
)

#text_slide(
  "15",
  "Concerted evolution of human subtelomeres - what we saw, predicted, and recovered",
  [Thesis: subtelomeres concertedly evolve through ongoing inter-chromosomal exchange - observable in pedigrees, predicted by 3D, recovered by an implicit pangenome graph across HPRC v2.]
)
