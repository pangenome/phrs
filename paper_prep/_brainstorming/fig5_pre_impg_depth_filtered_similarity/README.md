# Fig5 pre-IMPG depth-filtered similarity pilot

Narrow one-node pilot for direct-alignment recombination detection in the
WashU pedigree. The canonical Fig. 5/Extended Data input is not the raw
many:many PAF directly: the existing raw SweepGA/FastGA `-f32` PAF must first
be reduced with `sweepga --num-mappings 10:10 --scaffold-jump 0 --scoring ani`.
Without this `10:10` pre-filter, the whole-genome 2 kb IMPG similarity scan is
not tractable and is not the analysis used for the figure.

The workflow is:

1. Filter the existing raw SweepGA/FastGA `-f32` PAF with `sweepga` using
   `--num-mappings 10:10`, with scaffold chaining disabled. Older pilots also
   tested `1:1` and `4:4`.
2. Build 2 kb query-space windows from the child haplotype FASTA.
3. Drop exact CHM13 centromere windows and windows with no interchromosomal
   support or excessive interchromosomal PAF depth before calling IMPG.
4. Run `impg similarity` on the surviving BED only, then keep the best
   same-chromosome and best interchromosomal class winners per query window for
   homologous-vs-non-homologous recombination plotting. The older top-N output
   path is retained for auditing.

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

For the current class-winner path used by the whole-genome ribbon plots:

```bash
COMPARISON_ID=PAN027mat_vs_PAN010_joint BASES=10:10 RUN_TOPN=0 RUN_CLASS_WINNERS=1 \
  sbatch paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/scripts/run_one_node_pilot.sh

COMPARISON_ID=PAN028mat_vs_PAN027_joint BASES=10:10 RUN_TOPN=0 RUN_CLASS_WINNERS=1 \
  sbatch paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/scripts/run_one_node_pilot.sh
```
