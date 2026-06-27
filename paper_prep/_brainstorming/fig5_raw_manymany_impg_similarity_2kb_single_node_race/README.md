# Fig5 Raw Many:Many IMPG Similarity 2 kb Single-Node Race

This scratch pipeline races the existing sharded Fig5 raw many:many IMPG
similarity scan with six no-array Slurm jobs: one job per method and comparison.
It consumes only existing raw unfiltered many:many PAFs and the existing full
2 kb target BEDs recorded by
`paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_sharded/`.

It does not run WFMASH, SweepGA/FastGA, minimap2, seqwish, odgi, or any new
alignment. Each Slurm job calls `/home/erikg/.cargo/bin/impg similarity` once
against the full target BED and passes `--threads "${SLURM_CPUS_PER_TASK}"`.

Submit or regenerate manifests with:

```bash
python3 paper_prep/_brainstorming/fig5_raw_manymany_impg_similarity_2kb_single_node_race/scripts/build_single_node_race.py --submit
```

The build script submits a Slurm dependency finalizer after the six race jobs.
The finalizer writes completion and best-per-window plotting summaries under
`summaries/` after all six jobs finish.
