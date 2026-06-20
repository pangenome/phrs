# Direct sweepGA parental haplotype concordance

Scratch package for testing whether direct child-haplotype to transmitting-parent-haplotype `sweepGA`/FastGA alignments recover the same inheritance and recombination structure as the graph/`odgi untangle` analysis used for the Fig. 5 schematic.

This is review-only evidence under `paper_prep/_brainstorming/`; no manuscript or `submission/` files are edited.

## Inputs

- Window FASTA: `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/washu.1Mb.telo_500kb_trimmed.fa.gz`
- Query lists reused from prior native untangle/sweepGA scratch:
  `/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/queries_*.txt`
- Parent target lists reused from upstream untangle:
  `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/targets_*.txt`
- Graph-derived comparison targets:
  `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/event_manifest.tsv`
  and `selected_segments.tsv`
- Prior strict direct/untangle target:
  `paper_prep/_brainstorming/fig5_sweepga_1to1_redraw/conservative_segments.tsv`
  plus `/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm/`

All coordinates are native sample assembly/window coordinates parsed from FASTA record names and graph-derived selected segment tables. No reference-projected coordinates are introduced.

## sweepGA binary

The installed binary used for submission is `/home/erikg/.cargo/bin/sweepga`, resolving to `/export/local/home/erikg/.cargo/bin/sweepga`, with `sweepga --version` reporting `sweepga 0.1.1`.

The current CLI supports the required spelling directly:

```bash
sweepga --num-mappings many:many --scaffold-jump 0
```

Raw jobs use:

```bash
sweepga --fastga --num-mappings many:many --scaffold-jump 0 --temp-dir "$TMPDIR" --output-file "$TMPDIR/<comparison>.paf" query.fa target.fa
```

`TMPDIR` is a per-job directory under `${SLURM_TMPDIR:-/tmp}` and is removed by a shell `trap`. `/tmp` is used as node-local scratch when Slurm does not set `SLURM_TMPDIR`; an attempted `/dev/shm` run reached FastGA but failed during `FAtoGDB` temporary-file cleanup before producing alignments.

## Comparisons

Configured in `config/comparisons.tsv`:

- `PAN027pat_vs_PAN011_hap1`
- `PAN027pat_vs_PAN011_hap2`
- `PAN027mat_vs_PAN010_hap1`
- `PAN027mat_vs_PAN010_hap2`
- `PAN028mat_vs_PAN027_hap1`
- `PAN028mat_vs_PAN027_hap2`

These cover the requested transmitting-parent checks: PAN027 paternal product versus PAN011 hap1/hap2, PAN027 maternal product versus PAN010 hap1/hap2, and PAN028 maternal product versus PAN027 hap1/hap2. Inspecting the graph manifest did not reveal an additional directly required transmitting-parent comparison for this first pass.

## Commands

Run from this package directory:

```bash
python3 scripts/prepare_inputs.py
bash scripts/submit_raw_many_many.sh
squeue -j "$(tail -n +2 summaries/slurm_jobs.tsv | cut -f2 | paste -sd, -)"
```

After raw PAFs exist:

```bash
bash scripts/run_filter_matrix.sh
python3 scripts/summarize_paf.py
python3 scripts/plot_concordance.py
```

The filter matrix in `config/filter_matrix.tsv` preserves raw `many:many`/no-scaffold PAFs and derives:

- `1:1` no-scaffold
- `1:many` no-scaffold
- `2:many` no-scaffold
- `4:many` no-scaffold
- simple identity/length/query-coverage threshold output

Thresholds are parameterized in the TSV so simple PAF filters can be changed without rerunning expensive alignment.

## Outputs

- Extracted FASTA and name lists: `inputs/`
- Raw many:many/no-scaffold PAFs: `raw_paf/*.paf.gz`
- Derived filtered PAFs: `filtered_paf/*.paf.gz`
- Slurm and sweepGA logs: `logs/`
- Job table: `summaries/slurm_jobs.tsv`
- Input manifest: `summaries/input_manifest.tsv`
- Per-PAF summary: `summaries/paf_file_summary.tsv`
- Graph-vs-direct concordance: `summaries/direct_vs_graph_concordance.tsv`
- Review-only concordance SVG: `plots/direct_sweepga_concordance_review.svg`
- Report: `REPORT.md`
- Focused review plot: `plots/direct_sweepga_concordance_focused.{svg,pdf}`
- Full-window overview plot: `plots/direct_sweepga_full_genome_overview.{svg,pdf}`

`direct_vs_graph_concordance.tsv` includes `evidence_source_recommendation` so downstream organization can decide whether a clean direct sweepGA segment should become the primary Fig. 5/pedigree evidence source. Recommendations are intentionally per segment/event-role, not a blanket replacement of graph/untangle output.

## Current status

Final active Slurm jobs are recorded in `summaries/slurm_jobs.tsv` and completed successfully:

- `1704265` PAN027 paternal hap2 vs PAN011 hap1
- `1704266` PAN027 paternal hap2 vs PAN011 hap2
- `1704267` PAN027 maternal hap1 vs PAN010 hap1
- `1704268` PAN027 maternal hap1 vs PAN010 hap2
- `1704269` PAN028 maternal hap1 vs PAN027 hap1
- `1704270` PAN028 maternal hap1 vs PAN027 hap2

Use `sacct -j 1704265,1704266,1704267,1704268,1704269,1704270` to inspect the completed Slurm run. To resume or add comparisons, edit `config/comparisons.tsv`, run `python3 scripts/prepare_inputs.py`, then `bash scripts/submit_raw_many_many.sh`.

Initial submission `1704247`-`1704252` failed before alignment because Slurm ran the batch script from its spool copy and the wrapper inferred the package directory incorrectly. The wrapper now exports `PACKAGE_DIR` explicitly from `submit_raw_many_many.sh`; those failed jobs produced no raw PAFs and are not analysis artifacts.

Second submission `1704253`-`1704258` reached FastGA but failed in `FAtoGDB` when using `/dev/shm` scratch (`FATAL ERROR: failed to remove temporary file errno 2`). Empty 20-byte gzip placeholders were removed before resubmission with node-local `/tmp` scratch.

Third submission `1704259`-`1704264` showed that `sweepga --output-file -` is not stdout-streaming in this installed binary: FastGA completed, but stdout was empty and gzip placeholders were invalid analysis artifacts. The wrapper now writes a real temporary PAF under per-job scratch and compresses that file.

Fourth submission `1704265`-`1704270` is the valid completed raw run. It produced all six required compressed raw PAFs, the filter matrix, summaries, concordance table, and review-only SVG/PDF plots.
