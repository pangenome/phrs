# Fig5 pre-IMPG depth-filtered similarity pilot

Narrow one-node pilot for the paternal PAN027 haplotype against PAN011 father:

1. Filter the existing raw SweepGA/FastGA `-f32` PAF with `sweepga` using
   `--num-mappings 1:1`, `4:4`, and `10:10`, with scaffold chaining disabled.
2. Build 2 kb query-space windows from the child haplotype FASTA.
3. Drop exact CHM13 centromere windows and windows with no interchromosomal
   support or excessive interchromosomal PAF depth before calling IMPG.
4. Run `impg similarity` on the surviving BED only, then keep bounded top hits
   per query window for plotting.

This uses existing alignments only. It does not run FastGA/WFMASH alignment or
build a graph.

Submit:

```bash
bash paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/scripts/submit_one_node_pilot.sh
```

Default knobs:

- comparison: `PAN027pat_vs_PAN011_joint`
- window: `2000` bp
- max pre-IMPG interchromosomal PAF depth per window: `100`
- top post-IMPG rows per window: `20`
- Slurm: `tux`, `96` CPUs, one node, no array

