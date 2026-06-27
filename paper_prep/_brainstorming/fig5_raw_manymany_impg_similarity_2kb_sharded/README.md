# Fig5 Raw Many:Many IMPG Similarity 2 kb Sharded

This directory contains the replacement execution scaffold for Fig5 IMPG
similarity over exact full-genome 2 kb target windows.

Primary evidence is raw, unfiltered many:many PAF only:

- WFMASH updated-bin raw PAFs from
  `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/query_grid_filter_manifest.tsv`
- SweepGA/FastGA f32 raw PAFs from
  `paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency32/summaries/query_grid_chop_filter_manifest.tsv`
- Query and target FASTAs from
  `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/input_manifest.tsv`

The pipeline does not run WFMASH, SweepGA, FastGA, minimap2, seqwish, odgi, or
new graph construction. SweepGA PAFs are passed to IMPG through BGZF-normalized
copies from the previous full-BED attempt after validating them with
`bgzip -t`.

## Files

- `scripts/build_2kb_sharded_impg_pipeline.py` builds exact 2 kb target BEDs
  from target FASTA `.fai` files, splits them into shard BEDs, writes shard and
  Slurm manifests, and optionally submits Slurm arrays.
- `scripts/finalize_2kb_sharded_impg.py` validates completed shard metadata,
  assembles one compressed TSV per method x comparison with a single header,
  and writes the requested summary tables.
- `manifests/target_bed_shards.tsv` records every FAI-derived BED shard.
- `manifests/shard_manifest.tsv` records source raw PAF, IMPG alignment PAF,
  FASTAs, BED shard, command, Slurm IDs, IMPG version/path, and output path for
  every shard.
- `manifests/slurm_array_manifest.tsv` records array-level submission state.
- `manifests/shard_completion_manifest.tsv` is refreshed by the finalizer.

## Submitted Command Shape

Each array task executes:

```bash
/home/erikg/.cargo/bin/impg similarity \
  --alignment-files EXISTING_RAW_OR_BGZF_PAF \
  --target-bed SHARD_2KB.bed \
  --sequence-files QUERY.fa TARGET.fa \
  --gfa-engine poa \
  --no-merge \
  --num-mappings many:many \
  --scaffold-jump 0 \
  --threads "${SLURM_CPUS_PER_TASK}"
```

The manifest records the command text with the literal
`${SLURM_CPUS_PER_TASK}` placeholder; the Slurm script expands it and passes
exactly the job allocation to IMPG.

