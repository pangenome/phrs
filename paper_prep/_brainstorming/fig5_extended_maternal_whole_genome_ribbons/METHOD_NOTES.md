# Method Notes: Maternal Fig. 5 Candidates

These notes document the exact maternal whole-genome ribbon workflow for
Fig. 5C-D.

## Alignment pre-filter and IMPG scan

The maternal Fig. 5 candidates use existing raw SweepGA/FastGA f32
many:many PAFs. Before any IMPG similarity calculation, each raw PAF was reduced
with the required SweepGA pre-filter:

```bash
/home/erikg/.cargo/bin/sweepga \
  --num-mappings 10:10 \
  --scaffold-jump 0 \
  --scoring ani \
  --temp-dir "$SCRATCH" \
  --output-file "$FILTERED_TMP" \
  "$RAW_UNCOMPRESSED"
```

This 10:10 SweepGA reduction is a required pre-filter for the plotted maternal
outputs. The rendered figures do not use raw many:many PAFs directly.

The child query haplotype was tiled into 2 kb windows. Windows were removed
before IMPG if they overlapped exact CHM13 centromere intervals from
`data/chm13-annotations.bed`, had no interchromosomal PAF support, or exceeded
the maximum pre-IMPG interchromosomal PAF depth. The exact filter parameters
were `--window-size 2000`, `--min-depth 1`, `--max-depth 100`, and
`--interchrom-only`.

`impg similarity` was then run on only the surviving query-window BED, using the
10:10-filtered PAF, the child query FASTA, and the maternal donor FASTA:

```bash
/home/erikg/.cargo/bin/impg similarity \
  --alignment-files "$FILTERED_PAF" \
  --target-bed "$QUERY_BED" \
  --sequence-files "$QUERY_FASTA" "$TARGET_FASTA" \
  --gfa-engine poa \
  --no-merge \
  --num-mappings many:many \
  --scaffold-jump 0 \
  --threads "$THREADS"
```

The IMPG stream was reduced to class winners per 2 kb query window with
`filter_impg_similarity_class_winners.py`. For each retained query window, the
downstream ribbon plot compares the best interchromosomal/non-homologous class
winner against the best same-chromosome homologous class winner. Colored ribbons
are drawn only for windows where the interchromosomal winner beats the
homologous winner.

## Runtime and versions

The successful maternal jobs ran from `/moosefs/erikg/phrs` on Slurm partition
`tux`, one node, `--cpus-per-task=96`, and `--mem=700G`. The wrapper used
`THREADS="${SLURM_CPUS_PER_TASK:-96}"`, exported `RAYON_NUM_THREADS` and
`OMP_NUM_THREADS` to that value, passed `--threads "$THREADS"` to IMPG, and
used `bgzip -@ "$THREADS"` plus `pigz -p "$THREADS"` for compression.

Temporary uncompressed PAFs and SweepGA scratch files were created below
`/dev/shm` with per-job directories named
`fig5_pre_impg.<comparison>.<jobid>.*`; the wrapper removed scratch on exit.

Successful runtime metadata:

| Comparison | Slurm job | Host | SweepGA | IMPG | Threads |
| --- | ---: | --- | --- | --- | ---: |
| `PAN027mat_vs_PAN010_joint` | 1708164 | `tux05` | `/home/erikg/.cargo/bin/sweepga`, `sweepga 0.1.1` | `/home/erikg/.cargo/bin/impg`, `impg 0.4.1` | 96 |
| `PAN028mat_vs_PAN027_joint` | 1708167 | `tux05` | `/home/erikg/.cargo/bin/sweepga`, `sweepga 0.1.1` | `/home/erikg/.cargo/bin/impg`, `impg 0.4.1` | 96 |

The maternal monitor validated both class-winner TSVs with `gzip -t` and SHA256
checks. The successful depth summaries confirmed 2 kb windows,
`interchrom_only=1`, no high-depth windows, and 48,332 centromere windows
removed in each comparison.

## Ribbon display grouping

The renderer is:

```bash
bash paper_prep/_brainstorming/fig5_extended_maternal_whole_genome_ribbons/scripts/render_maternal_whole_genome_ribbons.sh
```

It calls
`paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft/scripts/plot_whole_genome_ribbon_draft.py`
for each maternal comparison.

Chromosomes are concatenated in chromosome order and scaled by native length.
The top and bottom tracks are the two maternal donor haplotypes; the middle
track is the child query haplotype.

Raw 2 kb windows are grouped only by exact end-to-end adjacency. A grouped run
requires the same child sequence and donor sequence, query endpoints that touch
exactly, and donor endpoints that touch exactly in a consistent donor direction.
The configured and audited merge gap is `end_to_end_merge_gap_bp = 0`. The
merge audits for both comparisons report
`max_absorbed_query_endpoint_gap_bp = 0` and
`max_absorbed_donor_endpoint_gap_bp = 0` for both the interchromosomal and
homologous layers.

Colored interchromosomal/non-homologous ribbons are drawn for grouped runs of
at least 10 kb with mean interchromosomal identity at least 0.95, after requiring
the interchromosomal class winner to beat the homologous class winner. The
homologous-context render adds same-chromosome mother-child chains in light gray
for grouped runs of at least 10 kb and identity at least 0.95. Display-only
minimum widths may be applied to homologous ribbons and endpoint marks for
legibility; the TSV coordinates, base-pair totals, and window counts retain the
exact grouped intervals.

## Rendered outputs

For each comparison, the directory contains SVG, PDF, and PNG renders for:

- `*.whole_genome_ribbon.*`: colored interchromosomal/non-homologous winners.
- `*.whole_genome_homologous_context_ribbon.*`: colored winners over light gray
  homologous context.

The conversion status files record `rsvg-convert version 2.54.5`. The SVG
viewBox is 3600 by 840 pixels, and the manuscript-facing PNG exports are 7200
by 1680 pixels.
