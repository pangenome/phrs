# Maternal Extended Data Fig. 5 whole-genome ribbons

This directory contains whole-genome homologous-competition ribbon renders for
the validated maternal 10:10 IMPG class-winner outputs from
`fig5-ed-maternal-10to10-impg-monitor`.

The renderer is the generalized Fig. 5 draft script:

```bash
bash paper_prep/_brainstorming/fig5_extended_maternal_whole_genome_ribbons/scripts/render_maternal_whole_genome_ribbons.sh
```

Inputs:

- `PAN027mat_vs_PAN010_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz`
- `PAN028mat_vs_PAN027_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz`
- query and target FASTA indexes from
  `/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs`

The filtering and display thresholds are inherited unchanged from the accepted
Fig. 5 whole-genome ribbon draft: 2 kb query windows, exact end-to-end grouping
only (`end_to_end_merge_gap_bp = 0`), homologous same-chromosome context chains
drawn in light gray when grouped runs are at least 10 kb and identity is at
least 0.95, and colored interchromosomal/non-homologous ribbons drawn only when
the interchromosomal class winner beats the homologous same-chromosome winner.

Each comparison subdirectory contains:

- `*.whole_genome_ribbon.svg`
- `*.whole_genome_homologous_context_ribbon.svg`
- `*.whole_genome_ribbon_runs.tsv`
- `*.whole_genome_ribbon_summary.tsv`
- `*.whole_genome_homologous_context_runs.tsv`
- `*.whole_genome_homologous_context_summary.tsv`
- `*.whole_genome_ribbon_merge_audit.tsv`
- `*.conversion_status.txt`

If `rsvg-convert` is available, the renderer also emits PDF and PNG siblings for
both SVGs. Otherwise the conversion status records that SVG-only output was
produced.
