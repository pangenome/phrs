# f32 query-grid chop/filter overlap audit

## Conclusion

The f32 query-grid chop/filter rerun preserves the no-overlap behavior seen in
the f16 query-grid audit. All nine completed f32 comparison/chop-length outputs
in `summaries/query_grid_chop_filter_manifest.tsv` pass the absolute query-grid
boundary audit: every interior `q_start`/`q_end` boundary is on `k * N`, and all
non-grid boundaries are raw-alignment edge clips matching `zs:i` or `ze:i`.
There are zero boundary violations.

The f32 mapping-frequency change does not recover additional chr3 support in
the Fig5 chr9q -> chr3q candidate windows. Relative to the f16 query-grid audit,
f32 is identical for PAN027 at 2 kb, 5 kb, and 10 kb; identical for PAN028 at
2 kb and 10 kb; and lower for PAN028 at 5 kb, where the f16 chr3 h2 chunk is no
longer retained. No candidate row has query-space redundant bp after filtering,
and no target-switch pair overlaps in query space. Thus, the mapping-frequency
effect is sparse and does not add chr3 evidence, while the chunk-phase/filtering
effect remains fixed by query-grid chopping.

## Artifacts

- `scripts/query_grid_overlap_audit.py`: reproducible generator for both TSVs.
- `summaries/query_grid_overlap_audit.tsv`: one row per completed f32
  comparison/chop length, including boundary counts, candidate-window overlap
  metrics, target-switch counts, chr3 support, exact retained chr3 chunks, and
  direct f16-vs-f32 comparison columns.
- `summaries/query_grid_overlap_audit_examples.tsv`: retained query-grid rows
  and row-start comparison examples for debugging exact candidate-window
  intervals.

The f16 baseline is read from:
`../pedigree_whole_genome_sweepga_fastga_frequency16/summaries/query_grid_overlap_audit.tsv`.

## Boundary Audit

All completed f32 rows have zero boundary violations.

| comparison | chop bp | chopped rows | non-grid starts | non-grid ends | violations |
|---|---:|---:|---:|---:|---:|
| PAN027mat_vs_PAN010_joint | 2000 | 29,850,810 | 1,468,432 | 1,468,778 | 0 |
| PAN027mat_vs_PAN010_joint | 5000 | 12,822,289 | 1,468,788 | 1,469,200 | 0 |
| PAN027mat_vs_PAN010_joint | 10000 | 7,146,577 | 1,468,930 | 1,469,304 | 0 |
| PAN027pat_vs_PAN011_joint | 2000 | 25,335,045 | 1,452,910 | 1,453,180 | 0 |
| PAN027pat_vs_PAN011_joint | 5000 | 11,009,828 | 1,453,305 | 1,453,501 | 0 |
| PAN027pat_vs_PAN011_joint | 10000 | 6,230,158 | 1,453,400 | 1,453,667 | 0 |
| PAN028mat_vs_PAN027_joint | 2000 | 47,650,121 | 1,590,136 | 1,590,923 | 0 |
| PAN028mat_vs_PAN027_joint | 5000 | 20,014,329 | 1,590,801 | 1,591,367 | 0 |
| PAN028mat_vs_PAN027_joint | 10000 | 10,798,729 | 1,590,922 | 1,591,507 | 0 |

The non-grid starts/ends are not failures. The audit only accepts them when the
chopped boundary equals the raw-edge tag carried by the PAF row; otherwise the
row would be counted in `boundary_violation_rows`.

## Candidate Windows

All candidate-window rows have `query_redundant_bp = 0`,
`chr3_query_redundant_bp = 0`, `query_overlap_max_depth = 1`, and
`target_overlap_switch_pairs = 0`. Adjacent target switches remain because the
filter keeps different non-overlapping query bins that map to different
target chromosome/haplotype groups; those switches do not overlap in query
space.

| event | chop bp | rows | target switches | chr3 rows | chr3 union bp | f16 chr3 union bp | f32-f16 chr3 union bp | f32 vs f16 |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| PAN027_chr9q_chr3q_PHR_candidate | 2000 | 214 | 12 | 2 | 4,000 | 4,000 | 0 | same |
| PAN027_chr9q_chr3q_PHR_candidate | 5000 | 93 | 3 | 1 | 5,000 | 5,000 | 0 | same |
| PAN027_chr9q_chr3q_PHR_candidate | 10000 | 47 | 3 | 1 | 10,000 | 10,000 | 0 | same |
| PAN028_chr9q_chr3q_PHR_candidate | 2000 | 191 | 51 | 3 | 6,000 | 6,000 | 0 | same |
| PAN028_chr9q_chr3q_PHR_candidate | 5000 | 80 | 22 | 1 | 5,000 | 10,000 | -5,000 | less |
| PAN028_chr9q_chr3q_PHR_candidate | 10000 | 45 | 14 | 1 | 10,000 | 10,000 | 0 | same |

## Exact chr3 Chunks

PAN027 chr9q -> chr3q:

| chop bp | retained chr3 chunks |
|---:|---|
| 2000 | `136168000-136170000->PAN011#joint#h2_chr3:205538813-205540813(+)`; `136172000-136174000->PAN011#joint#h2_chr3:205542815-205544815(+)` |
| 5000 | `136165000-136170000->PAN011#joint#h2_chr3:205535813-205540813(+)` |
| 10000 | `136160000-136170000->PAN011#joint#h2_chr3:205530811-205540813(+)` |

PAN028 chr9q -> chr3q:

| chop bp | retained chr3 chunks |
|---:|---|
| 2000 | `134710000-134712000->PAN027#joint#h1_chr3:201235052-201237052(+)`; `134750000-134752000->PAN027#joint#h1_chr3:201274685-201276685(+)`; `134830000-134832000->PAN027#joint#h1_chr3:201354720-201356720(+)` |
| 5000 | `134750000-134755000->PAN027#joint#h1_chr3:201274685-201279685(+)` |
| 10000 | `134750000-134760000->PAN027#joint#h1_chr3:201274685-201284685(+)` |

The PAN028 5 kb f16-only chr3 chunk that drops out at f32 was:
`134735000-134740000->PAN027#joint#h2_chr3:202402192-202407192(+)`.

## Interpretation

The f16-to-f32 mapping-frequency change increases the global chopped row counts,
but it does not increase chr3 candidate support after the same query-grid
chop/filter procedure. The only chr3 support change in the audited candidate
rows is a decrease: PAN028 at 5 kb goes from two f16 chr3 chunks and 10 kb
query-union bp to one f32 chr3 chunk and 5 kb query-union bp.

The chunk-phase/filtering behavior is stable. Query-grid chunks still enter
SweepGA as absolute query bins, and the post-filter candidate windows remain
non-overlapping. The old row-start fixed-offset artifact, where shifted h1/h2
rows could survive as overlapping query-space support, does not recur in f32.
