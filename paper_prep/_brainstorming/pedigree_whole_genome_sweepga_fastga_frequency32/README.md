# Fig5 whole-genome sweepGA/FastGA frequency-32 raw iteration

Date: 2026-06-24

This package is the explicit FastGA k-mer occurrence threshold 32 iteration for
the Fig5 PAN027/PAN028 whole-genome raw sweepGA/FastGA matrix. It is deliberately
separate from the frequency-16 package at
`paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/`,
so downstream query-grid chop/filter panels can compare f16 and f32 without
overwriting either raw PAF set.

Primary command shape:

```bash
/home/erikg/.cargo/bin/sweepga \
  --fastga \
  --fastga-frequency 32 \
  --num-mappings many:many \
  --scaffold-jump 0 \
  --temp-dir /dev/shm/... \
  --output-file ... \
  QUERY.fa TARGET.fa
```

The three comparisons and full whole-genome FASTA inputs mirror the f16 package:

- `PAN027pat_vs_PAN011_joint`
- `PAN027mat_vs_PAN010_joint`
- `PAN028mat_vs_PAN027_joint`

All raw alignment work was submitted through Slurm with `/dev/shm` scratch. The
command logs under `logs/` record the compute node, scratch path, explicit
sweepGA path, `which`, realpath, sha256, version/help text, `--check-fastga`,
and exact command line.

## Workflow

Raw frequency-32 Slurm jobs:

- `1706271`: `PAN027pat_vs_PAN011_joint` RAW_PAF_OK 12:23:02 on `octopus07`
- `1706272`: `PAN027mat_vs_PAN010_joint` RAW_PAF_OK 06:15:09 on `octopus08`
- `1706273`: `PAN028mat_vs_PAN027_joint` RAW_PAF_OK 13:03:38 on `octopus09`

The logs prove the required command semantics, for example:

```bash
/home/erikg/.cargo/bin/sweepga --fastga --fastga-frequency 32 \
  --num-mappings many:many --scaffold-jump 0 \
  --temp-dir /dev/shm/sg.freq32.<job>.<suffix> \
  --output-file /dev/shm/sg.freq32.<job>.<suffix>/<comparison>.frequency32.paf \
  /dev/shm/sg.freq32.<job>.<suffix>/q.fa \
  /dev/shm/sg.freq32.<job>.<suffix>/t.fa
```

## Result

All three raw frequency-32 whole-genome many:many scaffold-jump 0 PAFs exist,
validate with `gzip -t`, and have sha256 sidecars:

- `PAN027pat_vs_PAN011_joint.sweepga_frequency32_many_many_j0.paf.gz`: 5,161,676,206 bytes.
- `PAN027mat_vs_PAN010_joint.sweepga_frequency32_many_many_j0.paf.gz`: 5,665,017,941 bytes.
- `PAN028mat_vs_PAN027_joint.sweepga_frequency32_many_many_j0.paf.gz`: 9,744,517,614 bytes.

Raw chr3 support comparable to the f16 package:

- `PAN027_chr9q_chr3q_PHR_candidate`: 45 raw chr3 rows, 545,640 bp summed query-window overlap, 262,875 bp query-union coverage.
- `PAN028_chr9q_chr3q_PHR_candidate`: 43 raw chr3 rows, 546,514 bp summed query-window overlap, 262,821 bp query-union coverage.

Comparator summary:

- `PAN027_chr9q_chr3q_PHR_candidate`: frequency32 raw chr3 `yes`, prior sweepGA raw `no`, wfmash p95 `yes`.
- `PAN028_chr9q_chr3q_PHR_candidate`: frequency32 raw chr3 `yes`, prior sweepGA raw `no`, wfmash p95 `yes`.

This f32 iteration is slower and emits more raw mappings than the f16 package,
but it preserves the chr3 support needed for downstream f16-vs-f32 query-grid
chop/filter panel comparisons.

Required summaries:

- `summaries/slurm_jobs.tsv`
- `summaries/sweepga_binary.tsv`
- `summaries/fastga_binary.tsv`
- `summaries/input_manifest.tsv`
- `summaries/raw_chr3_support.tsv`
- `summaries/paf_file_summary.tsv`
- `summaries/frequency_sensitivity_summary.tsv`

Raw PAF outputs are named with `sweepga_frequency32_many_many_j0.paf.gz`; each
has a `.sha256` sidecar and is validated with `gzip -t`/`sha256sum -c`.

Validated with:

```bash
paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/scripts/validate_outputs.sh
```

## Query-grid chop/filter rerun

The f32 query-grid chop/filter rerun is separate from the f16 query-grid
package at
`paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/`.
It uses the f32 raw many:many PAFs directly from `raw_paf/`, not older chopped
outputs.

Required command shape:

```bash
pafchop --length <10000|5000|2000> --overlap 0 \
  --chunk-mode query-grid --comparison-id <comparison>.f32

/home/erikg/.cargo/bin/sweepga \
  --num-mappings 1:1 \
  --scaffold-jump 0 \
  --scoring ani \
  --overlap 0 \
  --output-file /dev/shm/.../filtered.paf \
  /dev/shm/.../input.paf
```

Outputs are distinct from f16:

- chopped PAFs: `chopped_paf_qgrid_l<N>_o0/`
- filtered PAFs: `filtered_paf_chop_sensitivity_query_grid/l<N>/`
- manifest: `summaries/query_grid_chop_filter_manifest.tsv`
- Slurm/status summary: `summaries/query_grid_chop_filter_slurm.tsv`
- shifted-boundary proof: `summaries/query_grid_shifted_boundary_audit.tsv`

Run status as of 2026-06-25:

- Slurm array `1706550` completed 8 of 9 required f32 cells on `octopus07`.
- Completed and `pigz -t` validated cells: all three lengths for
  `PAN027pat_vs_PAN011_joint`, all three lengths for
  `PAN027mat_vs_PAN010_joint`, and 10000/5000 bp for
  `PAN028mat_vs_PAN027_joint`.
- Missing cell: `PAN028mat_vs_PAN027_joint` at 2000 bp. Array element
  `1706550_9` was cancelled after 05:54:42 before finalizing the chopped PAF.
  A high-resource single-cell retry, array job `1706559_9`, was cancelled after
  00:31:23. A standard single-cell retry, job `1706560`, was also cancelled
  after 00:00:53.
- `summaries/query_grid_chop_filter_manifest.tsv` records this final cell as
  `MISSING_OR_INVALID`; all other rows are `OK`.

The f16 package remains the comparison point for the completed query-grid
sensitivity matrix, including the note that f16 1 kb was cancelled for runtime.
