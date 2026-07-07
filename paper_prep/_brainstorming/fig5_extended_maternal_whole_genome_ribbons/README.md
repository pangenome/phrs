# Maternal Fig. 5 whole-genome ribbons

This directory contains maternal companions for the Fig. 5 whole-genome
recombination ribbon view. The renders use the validated maternal
10:10 IMPG class-winner outputs from `fig5-ed-maternal-10to10-impg-monitor`:

- `PAN027mat_vs_PAN010_joint`: child maternal haplotype compared with mother
  PAN010.
- `PAN028mat_vs_PAN027_joint`: child maternal haplotype compared with mother
  PAN027.

These are the source renders for the maternal Fig. 5C-D panels. The
homologous-context PNGs are also copied into
`submission/fig/ExtendedDataFigures` as retained source assets, while the
manuscript-facing Fig. 5 uses the combined PNG in
`submission/fig/MainFigures/Fig5_whole_genome_recombination.png`.

## Build

The renderer is the generalized Fig. 5 draft script:

```bash
bash paper_prep/_brainstorming/fig5_extended_maternal_whole_genome_ribbons/scripts/render_maternal_whole_genome_ribbons.sh
```

The renderer reads class-winner IMPG outputs from:

- `PAN027mat_vs_PAN010_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz`
- `PAN028mat_vs_PAN027_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz`

The script currently resolves those inputs under
`/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs`
and query/target FASTA indexes under
`/moosefs/erikg/phrs/.wg-worktrees/agent-2636/paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/inputs`.

## Upstream IMPG method

The rendered inputs are not raw many:many alignment output. They come from the
Fig. 5 pre-IMPG depth-filtered similarity workflow:

1. Start from existing raw SweepGA/FastGA `-f32` many:many PAFs.
2. Apply the required SweepGA pre-filter before IMPG:
   `sweepga --num-mappings 10:10 --scaffold-jump 0 --scoring ani`.
3. Split the child query haplotype into 2 kb windows.
4. Remove exact CHM13 centromere-overlapping windows using
   `data/chm13-annotations.bed`.
5. Keep only query windows with interchromosomal support and with
   pre-IMPG interchromosomal PAF depth between 1 and 100 inclusive
   (`--min-depth 1`, `--max-depth 100`, `--interchrom-only`).
6. Run `impg similarity` on the surviving BED with `--gfa-engine poa`,
   `--no-merge`, `--num-mappings many:many`, `--scaffold-jump 0`, and the Slurm
   thread count.
7. Keep the class winners per 2 kb query window: the best same-chromosome
   homologous match and the best interchromosomal/non-homologous match.

For the successful maternal jobs, the runtime metadata recorded:

| Comparison | Slurm job | Host | Threads | SweepGA | IMPG | Scratch |
| --- | ---: | --- | ---: | --- | --- | --- |
| `PAN027mat_vs_PAN010_joint` | 1708164 | `tux05` | 96 | `/home/erikg/.cargo/bin/sweepga`, `sweepga 0.1.1` | `/home/erikg/.cargo/bin/impg`, `impg 0.4.1` | `/dev/shm/fig5_pre_impg.PAN027mat_vs_PAN010_joint.1708164.*` |
| `PAN028mat_vs_PAN027_joint` | 1708167 | `tux05` | 96 | `/home/erikg/.cargo/bin/sweepga`, `sweepga 0.1.1` | `/home/erikg/.cargo/bin/impg`, `impg 0.4.1` | `/dev/shm/fig5_pre_impg.PAN028mat_vs_PAN027_joint.1708167.*` |

Slurm requested one node on the `tux` partition with `--cpus-per-task=96` and
`--mem=700G`. The wrapper exported `RAYON_NUM_THREADS` and `OMP_NUM_THREADS` to
the same value and passed `--threads 96` to IMPG. Temporary uncompressed PAFs
and SweepGA temp files were staged in per-job `/dev/shm` directories, removed by
the wrapper's exit trap.

The maternal monitor validated the final class-winner outputs with `gzip -t`
and SHA256 checks, and confirmed that the class-winner skip reports contained
only header rows.

## Display method

The filtering and display thresholds are inherited unchanged from the accepted
Fig. 5 whole-genome ribbon draft:

- Chromosomes are concatenated in chromosome order with actual
  chromosome-length scaling within each genome track.
- The upper and lower tracks are the two mother donor haplotypes; the middle
  track is the child query haplotype.
- Colored ribbons are drawn from adjacent 2 kb query windows only when the
  interchromosomal/non-homologous class winner beats the best homologous
  same-chromosome class winner for that window.
- Colored ribbons are restricted to grouped runs of at least 10 kb with mean
  interchromosomal identity at least 0.95.
- The homologous-context render adds same-chromosome mother-child homologous
  chains in light gray when grouped runs are at least 10 kb and identity is at
  least 0.95.
- Raw 2 kb windows are grouped by exact end-to-end adjacency only:
  `end_to_end_merge_gap_bp = 0`. A run is merged only when the query endpoints
  touch exactly and the donor endpoints also touch exactly in a consistent
  donor direction for the same query and donor sequences.
- The plotted homologous ribbons and endpoint marks use a display-only minimum
  width for legibility. The run tables and summaries retain exact native query
  and donor coordinates.

The merge audits for both comparisons confirm zero absorbed query or donor
endpoint gaps in the interchromosomal and homologous layers.

## Current metrics

| Comparison | Interchrom 2 kb segments | End-to-end interchrom runs | Drawn high-confidence runs | Drawn high-confidence bp | Dominant class |
| --- | ---: | ---: | ---: | ---: | --- |
| `PAN027mat_vs_PAN010_joint` | 1,036 | 660 | 25 | 366,000 | acrocentric-acrocentric only |
| `PAN028mat_vs_PAN027_joint` | 4,290 | 2,400 | 104 | 2,070,000 | 102 acrocentric-acrocentric runs plus one chr5/chr1 candidate and one acrocentric-other run |

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

Both rendered comparisons converted successfully with `rsvg-convert version
2.54.5`; the SVG viewBox is 3600 by 840 pixels and the exported PNGs are
7200 by 1680 pixels.

## Companion notes

- `caption.md` contains the source caption note for the Fig. 5C-D maternal panels.
- `METHOD_NOTES.md` contains manuscript-ready methods text.
- `PROMOTION_DECISION.md` records the decision to fold the maternal panels into
  the main Fig. 5 manuscript figure.
