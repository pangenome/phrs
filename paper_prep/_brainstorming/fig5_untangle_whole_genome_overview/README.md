# Fig5 Strict Untangle Whole-Genome Overview

This package renders a whole-genome overview of the strict untangle/sweepGA
primary-path signal used for the Fig5 candidate panels. It is designed as a
genome-scale barcode, not as a full alignment ribbon plot.

## Inputs

- `paper_prep/_brainstorming/fig5_sweepga_1to1_redraw/conservative_segments.tsv`
  provides the whole strict `nb=1` primary-path geometry.
- `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/selected_segments.tsv`
  and `event_manifest.tsv` provide corrected callout labels and candidate
  roles for PAR1, PAN027 chr9q->chr3q, and PAN028 chr9q->chr3q.
- `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv`
  is optional metadata context only. Patch extents are not used to draw
  geometry.

## Coordinate System

All plotted and tabulated coordinates are native sample assembly window
coordinates parsed from the untangle query path names. They are 0-based,
half-open native windows and are not CHM13-projected coordinates.

For readability at whole-genome scale, the top barcode plot allocates one equal
block to each queried terminal arm and scales segments within that arm by local
position in the native terminal source window. Exact parsed native coordinates
are retained in `untangle_whole_genome_segments.tsv`.

## Figure Encoding

- Rows are transmitted child haplotype/comparison tracks.
- Columns are queried arms ordered from chr1p/chr1q through chrYp/chrYq where
  present in the strict table.
- Gray intervals remain on the same chromosome/arm primary path.
- Colored intervals are interchromosomal strict primary-path switches, colored
  by target chromosome/arm context.
- Candidate callout boxes mark PAR1, PAN027 chr9q->chr3q, and PAN028
  chr9q->chr3q in the genome-wide context.
- Side fragments and the low-confidence PAN027 tail are retained as dashed,
  lower-lane caveat intervals with triangular markers. They are not promoted to
  primary donor calls.

The lower panels summarize switch support as a compact query-arm vs target-arm
heatmap and as the largest query-arm interchromosomal support burdens.

## Outputs

- `fig5_untangle_whole_genome_overview.pdf`
- `fig5_untangle_whole_genome_overview.png`
- `fig5_untangle_whole_genome_overview.svg`
- `untangle_whole_genome_segments.tsv`
- `untangle_whole_genome_summary.tsv`
- `README.md`
- `validate_outputs.sh`

## Regeneration

Run from the repository root:

```bash
panel_dir=paper_prep/_brainstorming/fig5_untangle_whole_genome_overview
python3 "$panel_dir/scripts/extract_untangle_whole_genome_segments.py" --panel-dir "$panel_dir"
Rscript "$panel_dir/scripts/plot_untangle_whole_genome_overview.R" "$panel_dir"
bash "$panel_dir/validate_outputs.sh"
```

The generated TSVs and rendered assets are committed so downstream review does
not require rebuilding external untangle intermediates.
