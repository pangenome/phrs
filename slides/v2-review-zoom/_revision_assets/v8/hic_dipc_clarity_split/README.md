# Hi-C / Dip-C clarity split for review zoom v8

Task: `review-zoom-v8-hic-dipc-clarity-split`

This directory is a handoff package for the v8 fan-in render. It does not edit
the final Typst deck directly.

## What is staged here

- `inputs/chm13_hic_mds_3d_coords.png`: copied source CHM13 Hi-C MDS render.
- `plots/gm12878_mantel_proximity.{png,pdf}`: replacement GM12878 Mantel-style
  proximity plot with the negative-distance convention made explicit.
- `plots/sperm_all20_mantel_proximity.{png,pdf}`: replacement sperm proximity
  plot with the same sign convention.
- `inputs/gm12878_radial_community.png` and
  `inputs/sperm_all20_radial_community.png`: copied radial-community panels,
  staged here so slide 11a can be split without depending on v6 paths.
- `plots/wb_negative_control_reduced_text.{png,pdf}` and
  `plots/community_free_rho_distribution_reduced_text.{png,pdf}`: adjusted
  slide 11b/11c figure versions with smaller text.
- `plots/gm12878_cell01_whole_genome_3dg_projection.{png,pdf}`: optional
  whole-genome physical-coordinate context from one existing GM12878 3DG file.
- `hic_dipc_clarity_split_slides.typ` and `slides/hic_dipc_clarity_split_preview.pdf`:
  standalone previews of the recommended replacement pages.

## Slide 10m.2 audit answer

The current CHM13 3D MDS asset is neither PHR intervals nor 1 Mb flanks/windows.
It is a contact-space embedding computed from a broader CHM13 contact matrix over
whole p-arm, centromere, and q-arm intervals.

Concrete trace:

- The current deck uses `../_revision_assets/v7/hic_3d_plots/pngs/chm13_hic_mds_3d_coords.png`
  on slide `10m.2` in `slides/v2-review-zoom/_typst/zoom_review_deck.typ:683-688`.
  That staged PNG is byte-identical to
  `/moosefs/guarracino/HPRCv2/PHR_III/HiC/CHM13/plots/MDS_3d_coords.png`
  (`sha256 0d9a5f92ac3cc662efb1a7caad463c8181b10105de75a89ac245f618b5cd1d93`).
- The source BED is
  `/moosefs/guarracino/HPRCv2/PHR_III/HiC/CHM13/annotations/CHM13_chrom_parts.bed`.
  Its first chromosome is represented as `chr1 p`, `chr1 c`, `chr1 q`
  (`CHM13_chrom_parts.bed:1-3`), and every chromosome is similarly partitioned
  into p/centromere/q rows through `chrY` (`CHM13_chrom_parts.bed:70-72`).
  The file has 72 regions: 24 chromosomes times p/c/q.
- The analyzer parses the 4-column BED into region divisions
  (`/moosefs/guarracino/HPRCv2/submission_Randiak/scripts/analyzer.py:30-38`),
  uses those BED rows as `region_coords` and labels (`analyzer.py:95-104`),
  fetches each pair of regions from the cooler matrix and sums balanced contact
  values (`analyzer.py:129-133`), normalizes and returns the inter-region
  matrix (`analyzer.py:168-171`), converts contact to `1 - contact` distance
  (`analyzer.py:737`), and runs 3D MDS with `dissimilarity="precomputed"`
  (`analyzer.py:192-200`, `analyzer.py:741-750`).
- The project command block for this NOR/contact MDS analysis uses
  `scripts/hic/analyzer.py`, `--method MDS`, `--matrix-norm Sum`,
  `--bed-file $out/annotations/${sample}_chrom_parts.bed`, and
  `--norm-sizes $out/annotations/chrom.sizes`
  (`/moosefs/guarracino/HPRCv2/code.md:6039-6058` and duplicated in
  `/moosefs/guarracino/HPRCv2/code_PHR-and-3D.md:3313-3332`).
- NOR/acropcentric highlighting comes from
  `/moosefs/guarracino/HPRCv2/PHR_III/HiC/CHM13/annotations/CHM13_nor.chroms`,
  which lists `chr13`, `chr14`, `chr15`, `chr21`, and `chr22`
  (`CHM13_nor.chroms:1-5`).

Caption to use on slide 10m.2:

> CHM13 Hi-C MDS over whole p-arm, centromere, and q-arm contact regions from
> `CHM13_chrom_parts.bed`; not PHR intervals, not terminal 500 kb, and not 1 Mb
> PHR flanks/windows. This is a 3D contact-space summary, not a physical
> single-cell genome reconstruction.

## Whole-genome 3D feasibility

There is no existing rendered CHM13 whole-genome physical 3D reconstruction in
the audited slide assets. The CHM13 figure is an MDS embedding of contact
frequencies from a bulk `.mcool`, so it should not be described as a physical
genome structure.

Existing single-cell 3DG coordinate files do make a physical-coordinate context
plot feasible without rerunning hickit or remapping. I prepared one optional
candidate:

- `plots/gm12878_cell01_whole_genome_3dg_projection.png`
- source 3DG:
  `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/3dg/gm12878_01.impute3.round4.clean.3dg.gz`

Use this only as optional context. It is one GM12878 Dip-C cell, not CHM13 Hi-C,
and not a PHR/community validation summary. If the deck needs only the current
Hi-C visual answer, keep slide 10m.2 as the contact-space MDS with the honest
caption above.

## Slide 11a split

The old slide 11a placed four panels on one page. This package splits it into
four readable one-figure pages:

1. `11a.1`: GM12878 Dip-C proximity plot.
2. `11a.2`: GM12878 Dip-C radial-community plot.
3. `11a.3`: sperm single-cell Hi-C / 3DG proximity plot.
4. `11a.4`: sperm radial-community plot.

This keeps each panel large and labels the two radial plots by dataset and
measurement. The radial panels measure per-community normalized radial position
and within-vs-between radial-position similarity. The proximity panels measure
arm-level sequence similarity against 3D proximity.

## W/B definition

`W/B` means within-community divided by between-community distance:

`W/B = mean within-community 3D distance / mean between-community 3D distance`

For distance-based Dip-C/sperm panels, `W/B < 1` means within-community arms are
closer than between-community arms. For the negative-control slide, `S_all`
values above 1 mean non-sharing arms are farther apart, as expected.

## Distance and proximity sign convention

The original Mantel scatter was visually confusing because the plotted y-values
were negative distances while the axis text still said distance. The source
script already converts distance to proximity before computing/plotting the
Mantel result (`community_3d_enrichment.py:765-768`), but the axis labels still
say `Sequence similarity distance` and `Mean 3D distance`
(`community_3d_enrichment.py:855-856`).

The replacement `*_mantel_proximity` plots make the convention explicit:

- x-axis: `Sequence similarity (1 - Jaccard distance)`, so higher x means more
  shared sequence.
- y-axis: `3D proximity (-mean Euclidean distance)`, so higher y means closer.
- Positive rho now reads in the expected direction: more similar sequence is
  closer in 3D.

Report-backed metrics:

- GM12878 Dip-C: `W/B = 0.931`, 6.9% closer, Fisher `p = 2.4e-05`, Mantel
  `rho = 0.296`, `p = 0.002` (`end-to-end-report/report/06_dipc_validation.md:7-19`).
- Sperm single-cell Hi-C / 3DG: `W/B = 0.401`, 60% closer, Fisher
  `p = 3.9e-51`, Mantel `rho = 0.202`, `p = 0.023`
  (`end-to-end-report/report/06_dipc_validation.md:71-87`).

## Region definitions for Dip-C and sperm panels

The GM12878 T2T Dip-C analysis uses 16 cells after excluding cell 12 as a
duplicate of cell 10. Particles were selected using per-arm CHM13 PHR
coordinates for 38 C-community arms and terminal 500 kb regions for 8 arms
without CHM13 PHR support (`end-to-end-report/report/06_dipc_validation.md:7`).

The sperm analysis uses 20 haploid sperm cells and the corrected PHR-specific
single-cell 3D enrichment output (`end-to-end-report/report/06_dipc_validation.md:71-87`).
Source scripts and output roots are listed at
`end-to-end-report/report/06_dipc_validation.md:102-114`, including
`community_3d_enrichment.py`, the GM12878 output root, the sperm output root,
and `phr_and_500kb_regions.bed`.

## Validation notes

- Generated assets with `Rscript make_hic_dipc_clarity_assets.R`.
- Rendered standalone preview PDF and PNG pages with
  `typst compile hic_dipc_clarity_split_slides.typ slides/hic_dipc_clarity_split_preview.pdf`
  and `typst compile --ppi 144 hic_dipc_clarity_split_slides.typ slides/preview-{0p}.png`.
- The final deck was not edited.
