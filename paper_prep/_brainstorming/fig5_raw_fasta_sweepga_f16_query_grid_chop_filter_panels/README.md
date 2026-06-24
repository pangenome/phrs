# Fig5 query-grid chopped raw FASTA SweepGA/FastGA f16 panels

This package generates the corrective Fig5 query-grid panel set from the
merged f16 query-grid rerun and overlap audit outputs.

The displayed rows are:

- PAR1 X/Y positive control from `PAN027pat_vs_PAN011_joint`
- PAN027 chr9q -> chr3q candidate window from `PAN027pat_vs_PAN011_joint`
- PAN028 chr9q -> chr3q candidate window from `PAN028mat_vs_PAN027_joint`

For each row, the package shows the completed query-grid chop lengths: 10 kb,
5 kb, and 2 kb. The failed/cancelled 1 kb jobs are intentionally excluded.

## Inputs

The source summary files are read from:

- `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/summaries/query_grid_chop_filter_manifest.tsv`
- `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/summaries/query_grid_overlap_audit.tsv`

The filtered PAF paths are taken from the manifest. In this worktree the large
query-grid PAFs are not vendored under a relative directory, so the script uses
the manifest's resolved `/moosefs/.../agent-2700/...` PAF paths when present.
The package records both the manifest path and resolved path in
`query_grid_panel_manifest.tsv`.

## Query-grid chunking and filtering

These panels use PAF rows chopped by query coordinate grid, not by each source
alignment row's start. Query-grid chopping places each retained chunk on a fixed
grid for the query sequence at the requested length, with zero chunk overlap.
The PAF tags include `zm:Z:query-grid`, `zl:i:<length>`, and `zo:i:0`.

After chopping, the final SweepGA filter is the no-merge 1:1 configuration:

```text
SweepGA 1:1 --overlap 0 --scoring ani --scaffold-jump 0
```

That is why every summary and segment row carries:

- `filter_id=one_to_one_ani_o0`
- `num_mappings=1:1`
- `scaffold_jump=0`
- `scoring=ani`
- `filter_overlap=0`
- `chunk_mode=query-grid`

## Difference from the older row-start panels

The older row-start chop panels used chunks anchored to each input alignment
row. When homologous haplotypes had shifted row starts, this could leave
alternate haplotypes offset from one another and create apparent redundant chr3
support inside the same query window.

The query-grid rerun instead aligns all chunks to the same genomic query
coordinate grid before the 1:1 SweepGA filter. The audit table reports zero
boundary violations and zero retained chr3 query redundancy for the completed
2 kb, 5 kb, and 10 kb query-grid jobs. The figure uses genomic query
coordinates directly, not only relative 0-500 kb window coordinates.

## Outputs

- `fig5_raw_fasta_sweepga_f16_query_grid_chop_filter_panels.pdf`
- `fig5_raw_fasta_sweepga_f16_query_grid_chop_filter_panels.png`
- `fig5_raw_fasta_sweepga_f16_query_grid_chop_filter_panels.svg`
- `query_grid_panel_segments.tsv`
- `query_grid_panel_summary.tsv`
- `query_grid_panel_manifest.tsv`

Regenerate from the repository root:

```bash
panel_dir=paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_query_grid_chop_filter_panels
python3 "$panel_dir/scripts/extract_query_grid_panel_segments.py" --panel-dir "$panel_dir"
Rscript "$panel_dir/scripts/plot_query_grid_panel.R" "$panel_dir"
bash "$panel_dir/validate_outputs.sh"
```
