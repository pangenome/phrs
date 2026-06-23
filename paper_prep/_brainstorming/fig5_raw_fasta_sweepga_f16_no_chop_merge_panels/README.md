# Fig5 Raw FASTA SweepGA f16 No-Chop Merge Panels

This package generates candidate Fig5 evidence panels from the completed raw whole-genome
f16 SweepGA/FastGA many:many PAFs without running `pafchop`.

The panel comparison is deliberately narrow:

- raw no-chop `1:1`, no scaffold merge, ANI scoring;
- raw no-chop `1:1`, 50 kb scaffold merging, ANI scoring;
- raw no-chop `1:1`, 50 kb scaffold merging, log-length-ANI scoring.

Panel coordinates are plotted in absolute query-chromosome coordinates, not window-relative
coordinates. The summary table reports both summed and unioned expected-target overlap, and
the chr9q -> chr3q rows provide the chr3-survival comparison requested for the Fig5 audit.

## Inputs

Raw PAFs come from:

`/moosefs/erikg/phrs/.wg-worktrees/agent-2649/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/raw_paf`

Each raw PAF is the direct `sweepga --fastga-frequency 16 --num-mappings many:many --scaffold-jump 0`
whole-genome product recorded in that package's `summaries/slurm_jobs.tsv`.

## Rebuild

From the repository root:

```bash
PANEL_DIR=paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_no_chop_merge_panels
SOURCE_PACKAGE=/moosefs/erikg/phrs/.wg-worktrees/agent-2649/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16

bash "$PANEL_DIR/scripts/run_no_chop_merge_panels.sh" "$PANEL_DIR" "$SOURCE_PACKAGE"
bash "$PANEL_DIR/scripts/validate_outputs.sh" "$PANEL_DIR"
```

The runner writes:

- `filtered_paf/*.paf.gz`: raw no-chop sweepga-filtered PAFs for each comparison/mode;
- `raw_merge_panel_segments.tsv`: panel-overlapping filtered alignments;
- `raw_merge_panel_summary.tsv`: event/mode support summary and chr3 survival status;
- `raw_merge_panel_manifest.tsv`: source, command, checksum, and output provenance;
- `fig5_raw_fasta_sweepga_f16_no_chop_merge_panels.{pdf,svg}`: comparison panels.
