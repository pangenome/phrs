# Fig5 Raw Many:Many IMPG Similarity Scan

This directory is the only valid output area for the clean raw many:many IMPG
replacement of the invalid existing-PAF reducer task.

The pipeline uses the raw PAF evidence layers only:

- WFMASH updated-bin raw PAFs from
  `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/query_grid_filter_manifest.tsv`
- SweepGA/FastGA frequency-32 raw PAFs from
  `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/summaries/query_grid_chop_filter_manifest.tsv`
- query/target FASTA names from
  `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/input_manifest.tsv`

The old directory
`paper_prep/_brainstorming/fig5_whole_genome_existing_paf_impg_like_scan/` is
not used here except as a failure record.

## Current Status

The 10 kb full-BED run completed for all three WFMASH raw PAF comparisons.
SweepGA/FastGA f32 full-BED jobs `1706581` and `1706582` timed out at the
24-hour Slurm walltime, and sibling job `1706583` was cancelled after that
blocker was confirmed. Do not run the summary script until complete
SweepGA `.tsv.gz` outputs are regenerated; the current SweepGA TSV files are
partial or empty. See `REPORT.md` for exact job states and output paths.

## Commands

Generate full-genome 10 kb target BEDs, six Slurm scripts, and manifests:

```bash
python3 paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/generate_fullbed_impg_jobs.py
```

Submit the six primary jobs:

```bash
python3 paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/generate_fullbed_impg_jobs.py --submit
```

Summarize completed IMPG outputs:

```bash
python3 paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_scan/summarize_impg_similarity.py
```

## IMPG CLI Notes

The installed `/home/erikg/.cargo/bin/impg` is `impg 0.4.1`. Its
`similarity --help` documents `--merge-distance 0` and `--no-merge`, but the
binary rejects using both flags together. These jobs use `--no-merge`, which is
the exact accepted no-merging mode.

This IMPG build also requires `--sequence-files QUERY.fa TARGET.fa` and
`--gfa-engine poa` for `similarity` tabular output when using precomputed PAF
alignments. The generated Slurm scripts therefore extend the requested command
shape with those two required arguments while preserving raw PAF, full BED,
many:many mapping, no merge, no scaffold, and all allocated threads.

The WFMASH raw PAFs are already acceptable to IMPG. The SweepGA/FastGA f32 raw
PAFs are regular gzip streams; IMPG 0.4.1 requires BGZF for random-access PAF
indexing. The pipeline therefore writes BGZF-normalized copies under
`bgzf_raw_paf/` and records both the original `source_raw_paf` and the
`impg_alignment_paf` consumed by IMPG in `run_manifest.tsv`.
