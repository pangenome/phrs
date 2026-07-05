# Maternal 10:10 IMPG reconciliation, 2026-07-05

Task: `fig5-ed-maternal-10to10-impg-fixed`

## Decision

No Slurm resubmission was needed. The required canonical class-winner outputs
already exist under the shared repo root at:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs/PAN027mat_vs_PAN010_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz`
- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/outputs/PAN028mat_vs_PAN027_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz`

## Provenance

Required method confirmed from
`paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/scripts/run_one_node_pilot.sh`:

`f32 raw many:many SweepGA/FastGA PAF -> sweepga --num-mappings 10:10 --scaffold-jump 0 --scoring ani -> 2 kb query BED with centromere/depth filtering -> impg similarity class winners`

Successful jobs validated:

- `1708164` for `PAN027mat_vs_PAN010_joint`
- `1708167` for `PAN028mat_vs_PAN027_joint`

Runtime metadata files:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/metadata/runtime.1708164.txt`
- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/metadata/runtime.1708167.txt`

Recorded binary versions:

- `sweepga=/home/erikg/.cargo/bin/sweepga`
- `sweepga 0.1.1`
- `impg=/home/erikg/.cargo/bin/impg`
- `impg 0.4.1`

Both successful jobs ran on `tux05` with `slurm_cpus_per_task=96`.

## Validation

### Output integrity

`gzip -t` passed for both canonical output files:

| Comparison | Job | Gzip | Size bytes | Lines | SHA256 |
| --- | --- | --- | ---: | ---: | --- |
| `PAN027mat_vs_PAN010_joint` | `1708164` | `OK` | 13,514,559 | 342,316 | `f0609f47282589cdc8d827f504f50e4a7c012a3b924a43231b6f127002de6667` |
| `PAN028mat_vs_PAN027_joint` | `1708167` | `OK` | 14,004,182 | 340,454 | `123df76b8203728038a1fd796f970e0db1c3da8402419e6faf488031de0d27e4` |

### Manifest status

Successful manifests checked:

- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/summaries/pre_impg_depth_filtered_manifest.PAN027mat_vs_PAN010_joint.1708164.tsv`
- `/moosefs/erikg/phrs/paper_prep/_brainstorming/fig5_pre_impg_depth_filtered_similarity/summaries/pre_impg_depth_filtered_manifest.PAN028mat_vs_PAN027_joint.1708167.tsv`

Each manifest contains a single `10:10` row with:

- `status=OK`
- `topn_impg_similarity_tsv_gz=NA`
- the expected `class_winner_tsv_gz` path

### Supporting checks

- Filtered PAF gzip validation passed for both `10:10` `.noscaffold.ani.paf.gz`
  files.
- Depth summaries confirm `window_size=2000`, `min_depth=1`,
  `max_depth=100`, `interchrom_only=1`, and `high_depth_windows=0` for both
  comparisons.
- Class-winner skip reports contain header rows only, indicating no windows
  were dropped by the class-winner filter.

## Conclusion

The canonical maternal `10:10` IMPG class-winner outputs for
`PAN027mat_vs_PAN010_joint` and `PAN028mat_vs_PAN027_joint` are complete,
gzip-valid, and consistent with the successful Slurm runs `1708164` and
`1708167`. No rerun was required.
