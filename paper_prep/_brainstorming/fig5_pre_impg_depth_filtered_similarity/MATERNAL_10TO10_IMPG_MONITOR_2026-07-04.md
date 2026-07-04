# Maternal 10:10 IMPG monitor, 2026-07-04

Task: `fig5-ed-maternal-10to10-impg-monitor`

## Scope

Monitored the maternal class-winner IMPG scans submitted from commit
`993c8a4b0d6dc7097d196c2cae1cca206ab2970b` with:

- script: `paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/scripts/run_one_node_pilot.sh`
- method: existing f32 raw many:many SweepGA/FastGA PAF -> `sweepga --num-mappings 10:10 --scaffold-jump 0 --scoring ani` -> 2 kb query BED with centromere/depth filtering -> `impg similarity` class winners
- knobs: `BASES=10:10`, `RUN_TOPN=0`, `RUN_CLASS_WINNERS=1`, `MAX_POST_IMPG_CANDIDATES=5000`
- submit/output root observed in Slurm logs: `/moosefs/erikg/phrs`

No duplicate jobs were submitted while an existing viable job was running.

## Slurm states

Collected with:

```bash
sacct -j 1708164,1708165,1708166,1708167 \
  --format=JobID,JobName%45,State,ExitCode,Elapsed,NodeList%30,Submit,Start,End -P
```

| JobID | JobName | State | ExitCode | Elapsed | NodeList | Submit | Start | End |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1708164 | fig5ed-pan027mat | COMPLETED | 0:0 | 00:30:15 | tux05 | 2026-07-04T13:09:56 | 2026-07-04T13:09:56 | 2026-07-04T13:40:11 |
| 1708165 | fig5ed-pan028mat | CANCELLED by 1001 | 0:0 | 00:00:00 | None assigned | 2026-07-04T13:09:56 | None | 2026-07-04T13:38:59 |
| 1708166 | fig5ed-pan028mat-48 | FAILED | 1:0 | 00:02:28 | octopus02 | 2026-07-04T13:38:59 | 2026-07-04T13:38:59 | 2026-07-04T13:41:27 |
| 1708167 | fig5ed-pan028mat | COMPLETED | 0:0 | 00:31:05 | tux05 | 2026-07-04T13:42:01 | 2026-07-04T13:42:01 | 2026-07-04T14:13:06 |

`1708165` was cancelled while pending with zero elapsed time and no node assigned.
`1708166` was a short-lived replacement attempt on `octopus02`; its stderr ended
with `Error: No such file or directory (os error 2)` when `sweepga` attempted to
read the raw PAF staged in `/dev/shm`. The successful replacement for the
PAN028 maternal scan was `1708167` on `tux05`.

## Binary versions

The successful runtime metadata recorded:

| JobID | Host | sweepga | impg | Threads |
| --- | --- | --- | --- | --- |
| 1708164 | tux05 | `/home/erikg/.cargo/bin/sweepga`, `sweepga 0.1.1` | `/home/erikg/.cargo/bin/impg`, `impg 0.4.1` | 96 |
| 1708167 | tux05 | `/home/erikg/.cargo/bin/sweepga`, `sweepga 0.1.1` | `/home/erikg/.cargo/bin/impg`, `impg 0.4.1` | 96 |

## Output validation

All paths below are under
`/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity`.

| Comparison | Final output | gzip -t | Size bytes | Lines | SHA256 |
| --- | --- | --- | ---: | ---: | --- |
| `PAN027mat_vs_PAN010_joint` | `outputs/PAN027mat_vs_PAN010_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz` | OK | 13,514,559 | 342,316 | `f0609f47282589cdc8d827f504f50e4a7c012a3b924a43231b6f127002de6667` |
| `PAN028mat_vs_PAN027_joint` | `outputs/PAN028mat_vs_PAN027_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz` | OK | 14,004,182 | 340,454 | `123df76b8203728038a1fd796f970e0db1c3da8402419e6faf488031de0d27e4` |

Both outputs have the expected class-winner schema:

```text
chrom	start	end	group.a	group.b	group.a.length	group.b.length	intersection	jaccard.similarity	cosine.similarity	dice.similarity	estimated.identity	winner_class	query_seq	other_seq	query_chrom	other_chrom	raw_candidate_count	class_candidate_count
```

The intermediate filtered PAFs also passed gzip integrity checks:

- `filtered_paf/PAN027mat_vs_PAN010_joint.sweepga_f32.10to10.noscaffold.ani.paf.gz`
- `filtered_paf/PAN028mat_vs_PAN027_joint.sweepga_f32.10to10.noscaffold.ani.paf.gz`

## Manifest validation

The successful manifests exist, are non-empty, and contain one `10:10` row with
`status=OK`, `topn_impg_similarity_tsv_gz=NA`, and the expected class-winner
output path:

- `summaries/pre_impg_depth_filtered_manifest.PAN027mat_vs_PAN010_joint.1708164.tsv`
- `summaries/pre_impg_depth_filtered_manifest.PAN028mat_vs_PAN027_joint.1708167.tsv`

The depth summaries confirm 2 kb windows, `min_depth=1`, `max_depth=100`,
`interchrom_only=1`, and no high-depth windows:

| Comparison | PAF rows | Counted rows | Total windows | Kept windows | Zero-depth windows | High-depth windows | Centromere windows | Max observed depth | Mean kept depth |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `PAN027mat_vs_PAN010_joint` | 458,813 | 341,547 | 1,509,456 | 171,936 | 1,289,188 | 0 | 48,332 | 33 | 4.298826 |
| `PAN028mat_vs_PAN027_joint` | 473,579 | 341,974 | 1,518,303 | 173,098 | 1,296,873 | 0 | 48,332 | 48 | 4.270858 |

The class-winner skip reports contain only their header rows for both scans,
indicating no query windows were skipped after the class-winner filter.

## Conclusion

Both required maternal `10:10` class-winner IMPG similarity outputs exist and
pass `gzip -t`. The final successful Slurm jobs were `1708164` for
`PAN027mat_vs_PAN010_joint` and `1708167` for `PAN028mat_vs_PAN027_joint`.
