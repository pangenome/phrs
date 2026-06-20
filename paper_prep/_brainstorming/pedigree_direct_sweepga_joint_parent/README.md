# Joint-parent direct sweepGA Fig5 comparison

This package reruns the Fig5 direct sweepGA comparison as a joint parent-choice
problem. The previous direct package aligned each transmitting child haplotype
against parent hap1 and parent hap2 in separate jobs and then compared or
combined separately filtered PAFs. That is wrong for a parent-choice display:
`1:1`, `1:many`, `2:many`, and `4:many` filters must choose among both parental
haplotypes at the same time, otherwise a mapping can survive independently in
both hap-specific target sets and make the `1:1` result look more meaningful
than it is.

## Inputs and coordinate convention

Source FASTA:

`/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/washu.1Mb.telo_500kb_trimmed.fa.gz`

This is a telomeric-window FASTA, not a whole-genome alignment input. PAF query
and target offsets are local 0-based half-open coordinates inside each trimmed
0-500 kb telomeric source window. Native assembly coordinates in downstream
tables are reconstructed from FASTA record names for labels only; they are not
CHM13 projections.

The joint comparisons are:

| comparison_id | query FASTA | combined target FASTA |
| --- | --- | --- |
| `PAN027pat_vs_PAN011_joint` | PAN027 paternal hap2 windows | PAN011 hap1 + hap2 windows |
| `PAN027mat_vs_PAN010_joint` | PAN027 maternal hap1 windows | PAN010 hap1 + hap2 windows |
| `PAN028mat_vs_PAN027_joint` | PAN028 maternal hap1 windows | PAN027 hap1 + hap2 windows |

The target FASTA headers are collapsed to one sweepGA target genome group
(`PAN011#joint#...`, `PAN010#joint#...`, `PAN027#joint#...`) so sweepGA performs
a pairwise child-vs-combined-parent run rather than a three-genome all-vs-all
run. The original parent haplotype remains recoverable from the record label
(`haplotype1`, `haplotype2`, `maternal`, `paternal`) and is restored in the
schematic segment tables.

## Commands

From the repository root:

```bash
python3 paper_prep/_brainstorming/pedigree_direct_sweepga_joint_parent/scripts/prepare_inputs.py
bash paper_prep/_brainstorming/pedigree_direct_sweepga_joint_parent/scripts/submit_raw_many_many.sh
bash paper_prep/_brainstorming/pedigree_direct_sweepga_joint_parent/scripts/run_filter_matrix.sh
python3 paper_prep/_brainstorming/pedigree_direct_sweepga_joint_parent/scripts/summarize_paf.py
```

The heavy FastGA/sweepGA alignments were submitted through Slurm, not run on the
head node. Successful jobs:

| comparison_id | Slurm job | node | raw records |
| --- | ---: | --- | ---: |
| `PAN027pat_vs_PAN011_joint` | 1704274 | octopus07 | 72 |
| `PAN027mat_vs_PAN010_joint` | 1704275 | octopus07 | 64 |
| `PAN028mat_vs_PAN027_joint` | 1704276 | octopus07 | 694 |

The first submitted jobs (`1704271`-`1704273`) intentionally remain in `logs/`
as provenance for a failed attempt: sweepGA detected the combined target as two
target genome groups and failed during `ALNtoPAF`. The corrected header-collapse
preparation produced the successful jobs above.

## Outputs

Raw many:many/no-scaffold PAFs are preserved in `raw_paf/`. Joint filtered PAFs
are in `filtered_paf/` for:

- `many_many_noscaffold`
- `one_one_noscaffold`
- `one_many_noscaffold`
- `two_many_noscaffold`
- `four_many_noscaffold`
- `simple_i95_l1k_q80` as a carry-over summary threshold control

`filtered_paf/*many_many_noscaffold.paf.gz` is an exact raw-PAF copy for
filter-matrix completeness; the source-of-truth raw files remain in `raw_paf/`.

Summary tables:

- `summaries/input_manifest.tsv`
- `summaries/slurm_jobs.tsv`
- `summaries/paf_file_summary.tsv`
- `summaries/direct_vs_graph_concordance.tsv`

Fig5 schematic sibling outputs are in:

`paper_prep/_brainstorming/fig5_synteny_recombination_joint_parent/`

Key files:

- `selected_segments.joint_parent_one_one_noscaffold.tsv`
- `selected_segments.joint_parent_four_many_noscaffold.tsv`
- `selected_segments.joint_parent_many_many_noscaffold.tsv`
- `fig5_synteny_recombination_joint_parent_1to1_full.{svg,pdf}`
- `fig5_synteny_recombination_joint_parent_4many_full.{svg,pdf}`

The old schematic directory and the prior direct-sweepGA-1to1 sibling directory
were not overwritten.

## Filter interpretation

`1:1` is a diagnostic control, not a good display model for multimap PHRs. In
this joint run it still crushes donor signal: the selected Fig5 events contain
9 rendered rows, with no chr3q donor rows. It retains PAR1 (`chrYp`, 144,103 bp)
and same-chromosome context, but it collapses the autosomal PHR donor structure.

`4:many` is the best current schematic filter. It produces the same selected
Fig5 rows as raw many:many for the plotted events (17 rows), while keeping the
display less permissive than raw. It preserves PAR1 (`chrYp`, 144,103 bp) and
recovers a PAN028 chr3q donor segment (`chr3q`, 7,765 bp) plus a chr16p side
fragment (`14,898 bp`). It does not recover a PAN027 chr3q donor in the direct
joint PAF; the PAN027 chr9q query raw rows map to chr9q and chr9p only, so that
absence is not caused by downstream joint filtering.

Raw many:many/no-scaffold is useful as the audit layer and is preserved as the
primary alignment evidence. For the selected Fig5 schematic events, raw
many:many and `4:many` give the same rendered segment set. For package-wide
summaries, raw many:many keeps the most inter-arm structure, especially in
`PAN028mat_vs_PAN027_joint` (161 inter-arm records; 4,845,500 inter-arm query bp)
relative to joint `1:1` (25 inter-arm records; 1,143,287 inter-arm query bp).

Conclusion: use joint `4:many` for the visible Fig5 direct-sweepGA comparison,
keep joint raw many:many as the audit/provenance layer, and treat joint `1:1`
only as a diagnostic control for how strongly strict one-to-one filtering
suppresses multimap subtelomeric PHR signal.
