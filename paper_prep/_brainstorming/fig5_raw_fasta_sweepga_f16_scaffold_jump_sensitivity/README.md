# Fig5 raw FASTA f16 SweepGA scaffold-jump sensitivity

This package contains the final PAF-filtering sensitivity sweep for the Fig5
raw-FASTA f16 evidence. The source alignments are the current whole-genome raw
f16 many:many PAFs from:

`/moosefs/erikg/phrs/.wg-worktrees/agent-2649/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/raw_paf/*.sweepga_frequency16_many_many_j0.paf.gz`

The final filtering matrix is defined in `config/matrix.tsv` and explicitly
includes:

- `--scaffold-jump`: `0`, `10k`, `20k`, `50k`
- `--num-mappings`: `1:1`, `4:many`
- `--scoring`: `ani`, `log-length-ani`
- `--min-aln-length`: unset/default, `1000`, `5000`, `10000`
- `--scaffold-mass`: fixed at the SweepGA default `10k`
- `--overlap`: SweepGA default, not varied

The unfiltered/multiway baseline is summarized directly from the raw
`many:many` source PAFs as `raw_many_many_unfiltered`; no additional SweepGA
filtering is applied to that baseline.

## Reproduce

Build the task manifest:

```bash
paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_scaffold_jump_sensitivity/scripts/build_filter_tasks.py
```

Run the heavy final PAF filtering on Slurm with bounded parallelism. The array
runner decompresses each source PAF to `/dev/shm`, runs the updated
`/home/erikg/.cargo/bin/sweepga`, compresses the filtered PAF, and writes a
`.sha256` checksum next to each heavy output:

```bash
paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_scaffold_jump_sensitivity/scripts/submit_filter_array.sh
```

After all array tasks complete, summarize and render the panel:

```bash
paper_prep/_brainstorming/fig5_raw_fasta_sweepga_f16_scaffold_jump_sensitivity/scripts/finalize_outputs.sh
```

## Outputs

- `filter_tasks.tsv`: complete Slurm task manifest and command templates.
- `filtered_paf_manifest.tsv`: source and filtered PAF checksums plus command lines.
- `candidate_window_segments.tsv`: per-alignment rows overlapping the absolute query-coordinate candidate windows.
- `candidate_window_summary.tsv`: expected-target rows, expected-target sum/union bp, all target-chrom union bp, row counts, and chr3 survival status per cell.
- `target_chrom_breakdown.tsv`: per-target-chromosome row and union-bp breakdown.
- `figures/fig5_raw_fasta_sweepga_f16_scaffold_jump_sensitivity.{pdf,svg,png}`: compact heatmap/table panel encoding chr3 union bp and off-target-only cells.

Heavy `filtered_paf/*.paf.gz` files and Slurm logs are intentionally ignored by
git. Their paths and checksums are captured in `filtered_paf_manifest.tsv`.
