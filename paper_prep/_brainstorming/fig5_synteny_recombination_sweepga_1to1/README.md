# Fig5 Direct sweepGA 1:1 No-Scaffold Schematic Comparison

This directory is a sibling comparison to
`paper_prep/_brainstorming/fig5_synteny_recombination_schematic/`. It reuses
the same schematic renderer, layout constants, flow styling, native 0-based
half-open coordinate conventions, and event formatting as
`plot_synteny_recombination_schematic.py`. The original schematic directory is
not modified or overwritten.

## Regenerate

Run from the repository root:

```bash
python3 paper_prep/_brainstorming/fig5_synteny_recombination_sweepga_1to1/build_selected_segments_from_direct_paf.py
python3 paper_prep/_brainstorming/fig5_synteny_recombination_sweepga_1to1/plot_synteny_recombination_sweepga_1to1.py
```

The plot script also runs the builder automatically when
`selected_segments.sweepga_1to1.tsv` is absent.

## Outputs

- `selected_segments.sweepga_1to1.tsv` - selected-segments-like geometry table
  derived from direct sweepGA filtered PAF rows.
- `fig5_synteny_recombination_sweepga_1to1_full.svg` - full schematic SVG.
- `fig5_synteny_recombination_sweepga_1to1_full.pdf` - PDF converted from the
  SVG with Guix `librsvg` / `rsvg-convert` when available.
- `pdf_conversion_status.txt` - converter provenance for the rendered PDF.

## Inputs And Filters

Geometry is read only from the direct sweepGA filtered PAFs in
`paper_prep/_brainstorming/pedigree_direct_sweepga_concordance/filtered_paf/`
whose suffix is `one_one_noscaffold.paf.gz`:

- `PAN027pat_vs_PAN011_hap1.one_one_noscaffold.paf.gz`
- `PAN027pat_vs_PAN011_hap2.one_one_noscaffold.paf.gz`
- `PAN028mat_vs_PAN027_hap1.one_one_noscaffold.paf.gz`
- `PAN028mat_vs_PAN027_hap2.one_one_noscaffold.paf.gz`
- `PAN027mat_vs_PAN010_hap1.one_one_noscaffold.paf.gz`
- `PAN027mat_vs_PAN010_hap2.one_one_noscaffold.paf.gz`

Only the four event-relevant PAFs contribute rows to this schematic:
`PAN027pat_vs_PAN011_hap1`, `PAN027pat_vs_PAN011_hap2`,
`PAN028mat_vs_PAN027_hap1`, and `PAN028mat_vs_PAN027_hap2`. The PAN027 maternal
vs PAN010 files are listed here because they are part of the complete direct
1:1 no-scaffold filter package, but the current Fig5 event manifest has no
PAN010 event.

The filter is the direct package `filter_id=one_one_noscaffold`, implemented as
`sweepGA --num-mappings 1:1 --scaffold-jump 0` in
`../pedigree_direct_sweepga_concordance/scripts/filter_paf.py` and recorded in
`../pedigree_direct_sweepga_concordance/config/filter_matrix.tsv`. No identity,
alignment-length, query-coverage, graph path, `nb:i:1`, or CHM13/reference-space
projection filter is applied by this directory.

## Difference From Graph/Untangle `selected_segments.tsv`

The original
`../fig5_synteny_recombination_schematic/selected_segments.tsv` is a curated
graph/`odgi untangle` strict-path segment table. It represents selected
conservative intervals, recovers target-side intervals from a prior native
graph/untangle strict PAF, and partitions the child/query window into schematic
segments with roles such as same-chromosome context, primary donor, side
fragment, and low-confidence tail.

`selected_segments.sweepga_1to1.tsv` instead contains the native PAF alignment
rows emitted by the direct haplotype-to-parent-haplotype sweepGA comparison for
the same three event query windows. Query and target coordinates are parsed
directly from the PAF sequence names and local PAF start/end columns. They are
native assembly-window coordinates and are not projected into CHM13/reference
space. Direct PAF rows may overlap and do not form the same graph-derived
partition of the 500 kb query window.

Event roles in the direct table are assigned only for visual grouping in the
existing renderer:

- same target arm as the event query arm -> `same-chromosome context`
- PAR1 event target `chrYp` -> `PAR positive control`
- autosomal target `chr3q` -> `primary donor`
- autosomal target `chr20q` -> `low-confidence tail`
- other non-query target arms -> `side fragment`

Optional community/PHR annotation fields are copied from the original
`selected_segments.tsv` only by `(event_id, target_arm, event_role)` when a
matching label exists. They are not used to alter direct PAF geometry.

Current direct-row counts:

| event | rows | role totals from direct rows |
| --- | ---: | --- |
| `PAR1_XY_positive_control` | 56 | same-chromosome context 488,791 bp; PAR positive control 514,995 bp |
| `PAN027_chr9q_chr3q_PHR_candidate` | 8 | same-chromosome context 795,424 bp; side fragment 129,499 bp |
| `PAN028_chr9q_chr3q_PHR_candidate` | 10 | same-chromosome context 798,497 bp |

The autosomal totals exceed 500 kb where direct 1:1 rows overlap on the query
axis. This is intentionally left as direct PAF-native geometry rather than
post-processed into the graph/untangle selected-segment partition.
