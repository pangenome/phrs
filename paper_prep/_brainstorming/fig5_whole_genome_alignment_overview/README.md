# Fig5 whole-genome alignment overview

This package generates whole-genome overview panels for the Fig5 pedigree
alignment evidence. The main panel is not a candidate-only zoom: it lays out all
query chromosomes in native query-coordinate order and draws one retained-target
track per method / comparison / chop length.

## Inputs

The generator uses only already-produced upstream artifacts.

- Untangle strict/corrected overview:
  `paper_prep/_brainstorming/fig5_untangle_whole_genome_overview/untangle_whole_genome_segments.tsv`
- SweepGA/FastGA f16 query-grid chopped and SweepGA `1:1` ANI-filtered PAFs:
  `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/summaries/query_grid_chop_filter_manifest.tsv`
- SweepGA/FastGA f32 query-grid chopped and SweepGA `1:1` ANI-filtered PAFs:
  `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/summaries/query_grid_chop_filter_manifest.tsv`
- wfmash `-p95` updated-bin PAFs after the same query-grid chopping and
  SweepGA `1:1` ANI filtering:
  `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/query_grid_filter_manifest.tsv`
- Candidate callout coordinates copied from the existing query-grid panels:
  `paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_query_grid_chop_filter_panels/config/panel_windows.tsv`

The output manifest records the exact filtered PAF paths and whether each input
was present when the tables were generated:
`whole_genome_method_manifest.tsv`.

## Binning and scoring

PAF inputs are parsed after upstream query-grid chopping and upstream SweepGA
`--num-mappings 1:1 --scoring ani --overlap 0` filtering. This script does not
run raw alignment, raw PAF chopping, or raw many-to-many filtering.

For each retained alignment row, query coordinates are split into 1 Mb bins.
The script sums query-overlap support per retained target chromosome within each
bin. The plotted winner for a bin is the target chromosome with the largest
retained support. Same-query/target-chromosome bins are assigned to
`same_chromosome`; bins with no retained support are emitted explicitly as
`no_support`. The binned table also reports
`top_interchrom_target_chrom` and `top_interchrom_support_bp`, so expected PAR1
chrY support and chr9q/chr3q support remain visible even when same-chromosome
support is the dominant retained target in a 1 Mb bin.

Mean identity for a winning bin is support-weighted from PAF matches/alignment
length. Untangle rows use the strict/corrected segment identity already present
in the untangle source table.

The support matrix in `whole_genome_support_matrix.tsv` is query-arm by
target-arm, with query/target arms assigned as a native-coordinate p/q proxy by
chromosome midpoint. Exact cytoband centromere coordinates are not available for
all native sample assemblies represented in these PAFs, so these are labelled
`*_arm_proxy` rather than cytogenetic arms. The figure's compact heatmap
aggregates the full matrix to target chromosome totals by method/chop setting so
that the visual remains legible.

## Figure layout

The rendered figure has three stacked sections:

1. Whole-genome binned tracks. Every query chromosome is visible on the x-axis.
   Rows are method/comparison/chop settings. Neutral fill marks no retained
   support, gray marks same-chromosome retained support, and colored bins mark
   the dominant interchromosomal retained target chromosome. Recurrent
   interchromosomal targets are labelled directly above rows, so interpretation
   does not rely on color alone.
2. Compact interchromosomal support heatmap. Cell text gives support totals; the
   full query-arm by target-arm table is committed separately.
3. Secondary callouts for PAR1 and the chr9q/chr3q candidate windows, using the
   same query-grid coordinates as the current Fig5 query-grid panels. These
   callouts are subordinate to the whole-genome tracks.

## Regeneration

From this directory:

```bash
./run_overview.sh
./validate_outputs.sh
```

`run_overview.sh` invokes:

```bash
python3 scripts/build_whole_genome_alignment_overview.py
Rscript scripts/plot_whole_genome_alignment_overview.R .
```

The generator scans only the already-filtered PAFs listed in the upstream
manifests. The committed deliverables are the lightweight binned/matrix tables,
method manifest, validation script, and rendered PDF/PNG/SVG.

## Outputs

- `fig5_whole_genome_alignment_overview.pdf`
- `fig5_whole_genome_alignment_overview.png`
- `fig5_whole_genome_alignment_overview.svg`
- `whole_genome_binned_support.tsv`
- `whole_genome_support_matrix.tsv`
- `whole_genome_method_manifest.tsv`
- `validate_outputs.sh`
