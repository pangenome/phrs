# Query-grid chop/filter overlap audit

## Conclusion

The query-grid chopper fixed the row-start phase artifact for the completed
Fig5 f16 candidate-window reruns. For every completed comparison/chop length in
`summaries/query_grid_chop_filter_manifest.tsv` (`l10000`, `l5000`, `l2000`),
the chopped PAF boundary audit found zero rows whose query boundaries violate
the absolute `k * N` grid after allowing clipped first/last fragments at raw
alignment edges. The retained candidate-window rows after SweepGA
`--num-mappings 1:1 --overlap 0 --scoring ani` have zero query-space redundant
bp, and the old alternating h1/h2 offset-overlap pattern does not persist.

The remaining ambiguity is biological/input support, not chunk phase: query-grid
filtering keeps only sparse chr3 chunks, while most retained window sequence is
still chr9 self/near-self support. The prior row-start panels had many more chr3
rows at 2 kb/5 kb because shifted raw mappings were chopped on different phases;
query-grid chopping collapses those into non-overlapping absolute bins before
SweepGA filtering.

`l1000` was not rerun here. The existing status files report failed jobs for
`PAN027pat_vs_PAN011_joint.l1000` and `PAN027mat_vs_PAN010_joint.l1000`, and no
completed manifest row exists for `PAN028mat_vs_PAN027_joint.l1000`; those are
reported as missing/status-only rows in the audit TSV.

## Artifacts

- `summaries/query_grid_overlap_audit.tsv`: one row per completed
  event/comparison/chop length, plus status-only rows for missing failed `l1000`
  outputs. Columns include grid-boundary counts, query-overlap metrics,
  chr3 retained support, target chromosome/haplotype summaries, retained chr3
  coordinates, and prior row-start comparison metrics.
- `summaries/query_grid_overlap_audit_examples.tsv`: example retained rows from
  the query-grid filtered outputs and chr3 rows from the prior row-start panel
  table for debugging.
- `scripts/query_grid_overlap_audit.py`: reproducible generator for both TSVs.

## Boundary Check

Completed query-grid chopped PAFs all pass the absolute-bin boundary test:

| comparison | chop bp | chopped rows | non-grid starts | non-grid ends | violations |
|---|---:|---:|---:|---:|---:|
| PAN027mat_vs_PAN010_joint | 2000 | 15760745 | 657895 | 658030 | 0 |
| PAN027mat_vs_PAN010_joint | 5000 | 6699237 | 658052 | 658227 | 0 |
| PAN027mat_vs_PAN010_joint | 10000 | 3679079 | 658114 | 658277 | 0 |
| PAN027pat_vs_PAN011_joint | 2000 | 13058186 | 659814 | 660020 | 0 |
| PAN027pat_vs_PAN011_joint | 5000 | 5621262 | 660021 | 660138 | 0 |
| PAN027pat_vs_PAN011_joint | 10000 | 3140024 | 660066 | 660222 | 0 |
| PAN028mat_vs_PAN027_joint | 2000 | 26495029 | 733927 | 734316 | 0 |
| PAN028mat_vs_PAN027_joint | 5000 | 11038200 | 734217 | 734508 | 0 |
| PAN028mat_vs_PAN027_joint | 10000 | 5884305 | 734265 | 734584 | 0 |

The non-grid starts/ends are clipped raw-alignment edges. They are not
violations because each such boundary matches the `zs:i` or `ze:i` raw edge
tag emitted by `pafchop-rs`; all interior chunk boundaries are on the absolute
query grid.

## Candidate-Window Overlap Behavior

| event | chop bp | rows in window | chr3 rows | chr3 sum bp | chr3 union bp | chr3 haplotypes | chr3 chunks | prior row-start chr3 rows/sum/union |
|---|---:|---:|---:|---:|---:|---|---|---|
| PAN027_chr9q_chr3q_PHR_candidate | 2000 | 214 | 2 | 4000 | 4000 | h2:2 | 136168000-136170000->PAN011#joint#h2_chr3:205538813-205540813(+)<br>136172000-136174000->PAN011#joint#h2_chr3:205542815-205544815(+) | 38/76000/65562 |
| PAN027_chr9q_chr3q_PHR_candidate | 5000 | 93 | 1 | 5000 | 5000 | h2:1 | 136165000-136170000->PAN011#joint#h2_chr3:205535813-205540813(+) | 6/30000/27219 |
| PAN027_chr9q_chr3q_PHR_candidate | 10000 | 47 | 1 | 10000 | 10000 | h2:1 | 136160000-136170000->PAN011#joint#h2_chr3:205530811-205540813(+) | 1/10000/10000 |
| PAN028_chr9q_chr3q_PHR_candidate | 2000 | 189 | 3 | 6000 | 6000 | h1:3 | 134710000-134712000->PAN027#joint#h1_chr3:201235052-201237052(+)<br>134750000-134752000->PAN027#joint#h1_chr3:201274685-201276685(+)<br>134830000-134832000->PAN027#joint#h1_chr3:201354720-201356720(+) | 25/48009/43121 |
| PAN028_chr9q_chr3q_PHR_candidate | 5000 | 80 | 2 | 10000 | 10000 | h1:1;h2:1 | 134735000-134740000->PAN027#joint#h2_chr3:202402192-202407192(+)<br>134750000-134755000->PAN027#joint#h1_chr3:201274685-201279685(+) | 4/20000/20000 |
| PAN028_chr9q_chr3q_PHR_candidate | 10000 | 45 | 1 | 10000 | 10000 | h1:1 | 134750000-134760000->PAN027#joint#h1_chr3:201274685-201284685(+) | 2/17232/17232 |

For all completed candidate rows:

- `query_redundant_bp = 0` and `chr3_query_redundant_bp = 0`.
- `chr3_overlap_hap_switch_pairs = 0`.
- `alternating_h1_h2_offset_persists = no`.

This is the expected behavior if the old overlap was caused by row-start phase:
the 2 kb PAN027 row-start panel had 38 chr3 rows with 76,000 summed bp but only
65,562 query-union bp, and 10 overlapping h1/h2 switch pairs. Query-grid
filtering at 2 kb retains only two chr3 rows, both h2, with 4,000 summed bp and
4,000 union bp. PAN028 shows the same direction: row-start 2 kb had 25 chr3
rows, 48,009 summed bp, 43,121 union bp, and four overlapping h1/h2 switch
pairs; query-grid 2 kb retains three non-overlapping h1 chr3 chunks.

## Debug Examples

The examples table includes early query-grid rows and prior row-start chr3 rows
for direct inspection. Representative rows:

| event | chop bp | source | query interval | target |
|---|---:|---|---|---|
| PAN027_chr9q_chr3q_PHR_candidate | 2000 | query_grid_filtered | 135704825-135706000 | PAN011#joint#h2_chr9:133238040-133240040 |
| PAN027_chr9q_chr3q_PHR_candidate | 2000 | row_start_filtered_no_merge_ani | 136037407-136039407 | PAN011#joint#h1_chr3:202381955-202383955 |
| PAN027_chr9q_chr3q_PHR_candidate | 2000 | row_start_filtered_no_merge_ani | 136041407-136043407 | PAN011#joint#h2_chr3:205412541-205414541 |
| PAN028_chr9q_chr3q_PHR_candidate | 2000 | query_grid_filtered | 134710000-134712000 | PAN027#joint#h1_chr3:201235052-201237052 |

Use `summaries/query_grid_overlap_audit_examples.tsv` for more rows and source
line numbers.
