# Fig5 pre-IMPG depth-filtered similarity pilot

Narrow one-node pilot for the paternal PAN027 haplotype against PAN011 father:

1. Filter the existing raw SweepGA/FastGA `-f32` PAF with `sweepga` using
   `--num-mappings 1:1`, `4:4`, and `10:10`, with scaffold chaining disabled.
2. Build 2 kb query-space windows from the child haplotype FASTA.
3. Drop exact CHM13 centromere windows and windows with no interchromosomal
   support or excessive interchromosomal PAF depth before calling IMPG.
4. Run `impg similarity` on the surviving BED only, then keep bounded top hits
   per query window for plotting.
5. For the best-period Fig5 check, rerun only `impg similarity` from the
   existing filtered PAFs and pre-IMPG BEDs, retaining the best same-chromosome
   homolog candidate and the best interchromosomal candidate per 2 kb query
   window. The best-all reducer writes the winning class and
   `delta_interchrom_minus_homolog` so the query-grid panel can mark windows
   where interchromosomal similarity beats the homolog.

This uses existing alignments only. It does not run FastGA/WFMASH alignment or
build a graph.

Submit:

```bash
bash paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/scripts/submit_one_node_pilot.sh
```

Rerun only the best-period IMPG similarity from already-generated
SweepGA f32 PAFs and pre-IMPG query BEDs:

```bash
bash paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/scripts/run_best_period_similarity.sh
bash paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_query_grid_panels/scripts/make_panels.sh
```

Default knobs:

- comparison: `PAN027pat_vs_PAN011_joint`
- window: `2000` bp
- max pre-IMPG interchromosomal PAF depth per window: `100`
- top post-IMPG rows per window: `20`
- Slurm: `tux`, `96` CPUs, one node, no array
