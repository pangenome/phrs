# Fig5 strict untangle query-grid-style panels

This package renders a compact Fig5 panel set in the same visual language as
the query-grid raw FASTA panels, but the geometry is not raw FASTA f16/f32
query-grid output.

The displayed intervals are strict primary-path untangle/sweepGA geometry from:

- `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/selected_segments.tsv`
- `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/event_manifest.tsv`

The optional WashU `patches.tsv` metadata is already represented only as
upstream labels where it joined exactly to strict segments. Patch extents are
not used to draw panel geometry.

## Events

Rows are ordered to match the raw FASTA query-grid reference panel:

1. `PAR1_XY_positive_control`
2. `PAN027_chr9q_chr3q_PHR_candidate`
3. `PAN028_chr9q_chr3q_PHR_candidate`

The PAN028 row is the corrected review-facing strict chr9q path from
`PAN028#1#chr9.haplotype1:134380985-134880984_chr9_qarm`. It replaces the
older PAN028 chr3q side-fragment view and should not be interpreted as that
superseded panel.

## Geometry Rules

The figure uses native sample assembly query coordinates from
`selected_segments.tsv`:

- `native_query_start_0based`
- `native_query_end_0based_exclusive`
- `native_query_interval_0based_half_open`

All rows retain `nb=1` strict primary-path support. Permissive multimap or
nth-best rows are excluded, and raw FASTA query-grid PAF chunks are not used for
geometry. The output manifest records this as `raw_fasta_query_grid_source=not_used`.

Primary donor and same-chromosome context intervals are drawn as the main
compact row. Side fragments and the low-confidence PAN027 tail are preserved as
dashed caveat markers below the row, not as equivalent primary donors.

## Outputs

- `fig5_untangle_query_grid_style_panels.pdf`
- `fig5_untangle_query_grid_style_panels.png`
- `fig5_untangle_query_grid_style_panels.svg`
- `untangle_panel_segments.tsv`
- `untangle_panel_summary.tsv`
- `untangle_panel_manifest.tsv`
- `validate_outputs.sh`

Regenerate from the repository root:

```bash
panel_dir=paper_prep/_brainstorming/fig5_untangle_query_grid_style_panels
python3 "$panel_dir/scripts/extract_untangle_panel_segments.py" --panel-dir "$panel_dir"
Rscript "$panel_dir/scripts/plot_untangle_query_grid_style_panel.R" "$panel_dir"
bash "$panel_dir/validate_outputs.sh"
```

The generated TSVs are intentionally committed with the rendered assets so the
panel can be reviewed without rebuilding external untangle intermediates.
