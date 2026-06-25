# Fig5 Whole-Genome Graph Untangle Pilot

This directory records the pilot graph/untangle workflow for `PAN027` paternal hap2
against the `PAN011` joint parent. The requested primary input is the existing
filtered WFMASH 2 kb query-grid 1:1 ANI PAF, built from the same whole-genome
FASTAs that WFMASH aligned.

## Inputs

`config/pilot_sources.tsv` is the source manifest. The primary run is:

- `run_id`: `wfmash_p95_updated_bin_l2000`
- comparison: `PAN027pat_vs_PAN011_joint`
- query FASTA: `/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027pat_vs_PAN011_joint.query.fa`
- target FASTA: `/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs/PAN027pat_vs_PAN011_joint.target.fa`
- filtered PAF: `/moosefs/erikg/phrs/.wg-worktrees/agent-2719/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/query_grid_filter/filtered_paf_qgrid_l2000_o0/PAN027pat_vs_PAN011_joint.one_to_one_ani_o0.chopped_l2000_o0_query_grid.paf.gz`
- filtered PAF SHA256: `35ec3bd8a491777444ce8e81c8ae817c96cd2657058953d184ebd62677ae2622`
- filtered row count: `1,507,938`

The manifest also includes `sweepga_fastga_f32_l2000`, which uses the analogous
SweepGA/FastGA frequency-32 filtered PAF. That row is a follow-on recipe, not the
first validity gate.

## Workflow

Validate the setup:

```bash
bash paper_prep/_brainstorming/fig5_whole_genome_graph_untangle_pilot/scripts/validate_pilot_setup.sh
```

Submit the WFMASH graph pilot:

```bash
sbatch paper_prep/_brainstorming/fig5_whole_genome_graph_untangle_pilot/scripts/run_graph_untangle_one.sbatch wfmash_p95_updated_bin_l2000
```

If the WFMASH graph builds cleanly and `summaries/wfmash_p95_updated_bin_l2000.focus_summary.tsv`
shows interpretable support, submit the f32 comparison using the same code path:

```bash
sbatch paper_prep/_brainstorming/fig5_whole_genome_graph_untangle_pilot/scripts/run_graph_untangle_one.sbatch sweepga_fastga_f32_l2000
```

The Slurm script stages FASTAs and the decompressed PAF under
`work/tmp` by default, verifies the filtered PAF SHA256 and row count when
recorded, and runs:

```bash
seqwish --seqs seqs.txt --paf-alns filtered.paf --gfa RUN.gfa --threads N --temp-dir SCRATCH --show-progress
odgi build --gfa RUN.gfa --out RUN.og --threads N --progress
odgi paths --idx RUN.og --list-paths
odgi untangle --idx RUN.og --query-paths query_paths.txt --target-paths target_paths.txt --merge-dist 50000 --n-best 4
odgi untangle --idx RUN.og --query-paths query_paths.txt --target-paths target_paths.txt --merge-dist 50000 --n-best 4 --paf-output
```

Heavy outputs are ignored by git under `work/`. Runtime command manifests and
focused summaries are also ignored while Slurm is writing them; inspect them in
`summaries/` after the job completes, then promote stable result tables in a
separate commit if they are needed as permanent figure inputs.

To force node-local scratch for a smaller rerun, submit with
`GRAPH_SCRATCH_ROOT=/dev/shm`. The first pilot attempts showed that `/dev/shm`
was not reliable for concatenating the full PAN027 paternal plus PAN011 joint
FASTAs on the worker node, so the committed default favors capacity over speed.

## Focus Checks

`config/focus_patterns.tsv` defines the comparison bins:

- `chr9q_chr3q_candidate`: PAN027 paternal chr9 path against PAN011 chr3 paths.
- `chr9q_native_context`: PAN027 paternal chr9 path against PAN011 chr9 paths.
- `par_xy_positive_control`: PAN027 paternal chrX path against PAN011 chrY paths.
- `acrocentric_controls`: PAN027 paternal acrocentric chromosomes against PAN011
  acrocentric chromosomes.

Interpretation rule for this pilot:

1. The graph is valid only if `seqwish`, `odgi build`, path listing, and untangle
   all complete without path-name loss.
2. The chr9q/chr3q candidate is useful only if it is not overwhelmed by generic
   all-to-all acro/PAR-like projections and if native chr9 support remains the
   major context.
3. PAR and acro bins are controls, not Fig5 evidence by themselves.
4. Run the f32 SweepGA/FastGA row only after the WFMASH graph passes the validity
   gate, so both methods are compared under the same graph/untangle extraction
   logic.
