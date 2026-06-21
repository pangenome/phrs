# Fig5 whole-genome wfmash -p95 evidence review

Date: 2026-06-21

Task: `fig5-whole-genome-wfmash-p95-evidence-review`

Reviewed package:

`paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95/`

Prior comparison points:

- `paper_prep/_brainstorming/fig5_whole_genome_sweepga_evidence_review/REPORT.md`
- `paper_prep/_brainstorming/fig5_whole_genome_sweepga_evidence_review/segment_support.tsv`
- `paper_prep/_brainstorming/fig5_whole_genome_sweepga_closeout/QA_REPORT.md`

## Direct answer

Yes. Raw whole-genome wfmash `v0.24.2` with literal `-p 95` recovered chr3-target PAF rows overlapping both Fig5 chr9 candidate windows before any filtering or one-to-one reduction:

- `PAN027_chr9q_chr3q_PHR_candidate`: 2 chr3 rows, 73,000 bp summed and union query overlap.
- `PAN028_chr9q_chr3q_PHR_candidate`: 3 chr3 rows, 152,860 bp summed query overlap and 150,862 bp union query overlap.

This changes the raw direct-alignment evidence status relative to the prior whole-genome sweepGA/FastGA raw `many:many -j 0` review. The prior raw sweepGA table had no chr3 rows for either autosomal chr9 candidate window; it had chr9 only for PAN027 and chr9 plus one chr16 side fragment for PAN028.

This is a technical provenance and raw-evidence review only. I did not modify `submission/` and did not create or update a Fig5 schematic.

## Technical provenance

### Whole-genome inputs

PASS. The wfmash package used the same three pedigree comparison IDs as the corrected sweepGA/FastGA package:

- `PAN027pat_vs_PAN011_joint`
- `PAN027mat_vs_PAN010_joint`
- `PAN028mat_vs_PAN027_joint`

The configured source FASTAs are the recovered full WashU v1.1 assemblies under `/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga/`, not window-only, chromosome-only, arm-only, or historical trimmed telomeric inputs. See `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95/config/comparisons.tsv:2`, `:3`, and `:4`.

The prepared input manifest records full haplotype FASTA scope for all three rows: 23 full-chromosome query records and 46 joint-parent target records where both target haplotypes are combined. The scope field explicitly says "full whole-genome haplotype FASTA records from source assembly .fai; no telomeric-window FASTA" on `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95/summaries/input_manifest.tsv:2`, `:3`, and `:4`.

### Binary and CLI semantics

PASS. The primary completed PAFs were generated with a source-built upstream current release, not the stale local Guix `0.12.5` binary. `summaries/wfmash_binary.tsv:2` records:

- role: `primary_binary`
- version: `v0.24.2-0-g774c01ff`
- release tag: `v0.24.2`
- tag commit: `774c01ffb9df010d6b529520033ed7dce0cb95d5`
- status: `built_from_source`

The same table records that the downloaded GitHub Linux binary was rejected for host glibc incompatibility, not used as a primary binary (`summaries/wfmash_binary.tsv:4`).

The captured `v0.24.2` help confirms current CLI semantics: `-s` is sketch size, `-w` is window size, `-p` is minimum mapping identity, and `-B` is temp directory. See `logs/wfmash_v0.24.2.help.txt:12`, `:13`, `:17`, and `:61`. The literal result reviewed here does not rely on `-s` or the permissive `-w 1k` configuration.

### Commands, -p 95, and scratch

PASS. `config/wfmash_parameters.tsv:2` defines the primary `literal_p95` parameter set as exactly `-p 95`. The three completed primary job rows in `summaries/wfmash_jobs.tsv:2`, `:4`, and `:6` all record:

- run label: `current_v0.24.2`
- parameter set: `literal_p95`
- status: `COMPLETED_OUTPUT_PRESENT`
- `wfmash_options`: `-p 95`
- `scratch_base`: `/dev/shm`
- binary: the source-built `tools/wfmash-v0.24.2-built/wfmash`

The exact command logs independently confirm `-p 95` and `-B /dev/shm/.../tmp`:

- `logs/current_v0.24.2.PAN027pat_vs_PAN011_joint.literal_p95.1704318.command.log:17`
- `logs/current_v0.24.2.PAN027mat_vs_PAN010_joint.literal_p95.1704320.command.log:17`
- `logs/current_v0.24.2.PAN028mat_vs_PAN027_joint.literal_p95.1704322.command.log:17`

The same logs record `/dev/shm` scratch and node-local staged target/query FASTAs on lines `:8` through `:11`, plus the output PAF and checksum on lines `:12`, `:18`, and `:19`.

The execution script matches those logs: it defaults `WFMASH_SCRATCH_BASE` to `/dev/shm`, stages `target.fa` and `query.fa` under the job scratch directory, runs wfmash with `-B "$TMP_BASE"`, bgzips the local PAF, and copies only the gzipped raw PAF back to package storage (`scripts/run_wfmash_one.sh:13`, `:41`-`:47`, `:75`-`:82`, and `:84`-`:86`).

### Raw PAF handling

PASS. The primary evidence below is from the raw gzipped PAFs listed in `summaries/wfmash_jobs.tsv`, not from filtered, one-to-one, or schematic-selected subsets. The review inspected the manifest-listed raw PAFs directly at their absolute paths under the upstream execution worktree:

- `PAN027pat_vs_PAN011_joint.literal_p95.wfmash-v0.24.2.paf.gz`, SHA-256 `96312a525e02e7e3ed0f0a73ee3fda22df2977e0c213789084952bcf8ed2c0af`.
- `PAN027mat_vs_PAN010_joint.literal_p95.wfmash-v0.24.2.paf.gz`, SHA-256 `6e63ec7bbc172f1ec334e6de6c8eb018876a9010ee3ad3c2522522bad3da6b8e`.
- `PAN028mat_vs_PAN027_joint.literal_p95.wfmash-v0.24.2.paf.gz`, SHA-256 `c20e9d49a958faf89ab71198ddd3d961e2d23a86803e9865f01d9e5d237db2c4`.

The checksums I recomputed from the raw PAF files matched `summaries/wfmash_jobs.tsv:2`, `:4`, and `:6` and the three command logs. The middle PAN027 maternal comparison is relevant for provenance completeness but not for either chr9/chr3 candidate window in `config/candidate_windows.tsv`.

## Candidate-window evidence

The candidate windows reviewed were:

- `PAN027_chr9q_chr3q_PHR_candidate`: `PAN027#2#chr9:135704825-136204825`, expected target chromosome `chr3` (`config/candidate_windows.tsv:2`).
- `PAN028_chr9q_chr3q_PHR_candidate`: `PAN028#1#chr9:134380985-134880985`, expected target chromosome `chr3` (`config/candidate_windows.tsv:3`).

The package summary already reports chr3 support for both events (`summaries/candidate_window_support.tsv:2` and `:7`). I independently streamed the manifest-listed raw PAFs and reselected rows by:

1. exact query contig name;
2. positive overlap between raw PAF query interval and the candidate window;
3. target chromosome parsed from the target contig name;
4. `target_chrom == chr3`.

The independent spot-check reproduced the package summary counts:

| Event | Raw wfmash PAF | chr3 rows | Query overlap bp sum | Query overlap bp union |
| --- | --- | ---: | ---: | ---: |
| `PAN027_chr9q_chr3q_PHR_candidate` | `PAN027pat_vs_PAN011_joint.literal_p95.wfmash-v0.24.2.paf.gz` | 2 | 73,000 | 73,000 |
| `PAN028_chr9q_chr3q_PHR_candidate` | `PAN028mat_vs_PAN027_joint.literal_p95.wfmash-v0.24.2.paf.gz` | 3 | 152,860 | 150,862 |

The minimal raw PAF rows to use for any future schematic are also written to:

`paper_prep/_brainstorming/fig5_whole_genome_wfmash_p95_evidence_review/minimal_chr3_raw_paf_rows.tsv`

Those rows contain the target contig names, query and target intervals, strand, raw PAF match/block/mapping-quality fields, and the available wfmash identity/score-like optional tags (`gi`, `bi`, `md`, `ch`).

### PAN027 chr9 candidate

Raw whole-genome wfmash `-p 95` emits two chr3-target rows overlapping `PAN027#2#chr9:135704825-136204825`.

| Query interval | Target interval | Strand | Window overlap | PAF matches/block | MapQ | wfmash tags |
| --- | --- | --- | ---: | --- | ---: | --- |
| `PAN027#2#chr9:136166000-136191000` | `PAN011#joint#h1_chr3:202510328-202535329` | `+` | 25,000 bp | `24951/25001` | 26 | `gi=0.99816; bi=0.997721; md=0.9984; ch=2028.26222.496` |
| `PAN027#2#chr9:136117000-136165000` | `PAN011#joint#h1_chr3:202461289-202509327` | `+` | 48,000 bp | `47971/48038` | 28 | `gi=0.999292; bi=0.998439; md=0.9991; ch=2028.26222.495` |

These two rows are adjacent on the query and target axes and both fall inside the PAN027 candidate window.

### PAN028 chr9 candidate

Raw whole-genome wfmash `-p 95` emits three chr3-target rows overlapping `PAN028#1#chr9:134380985-134880985`.

| Query interval | Target interval | Strand | Window overlap | PAF matches/block | MapQ | wfmash tags |
| --- | --- | --- | ---: | --- | ---: | --- |
| `PAN028#1#chr9:134764000-134814998` | `PAN027#joint#h2_chr3:202431105-202482361` | `+` | 50,998 bp | `50968/51256` | 22 | `gi=0.999373; bi=0.994323; md=0.9991; ch=4165.1.1` |
| `PAN028#1#chr9:134813000-134863981` | `PAN027#joint#h2_chr3:202480363-202531372` | `+` | 50,981 bp | `50931/51009` | 28 | `gi=0.998921; bi=0.998295; md=0.9992; ch=4167.1.1` |
| `PAN028#1#chr9:134709119-134760000` | `PAN027#joint#h1_chr3:201234027-201284685` | `+` | 50,881 bp | `50356/50881` | 18 | `gi=0.999008; bi=0.984554; md=0.9981; ch=4163.1.1` |

Two PAN028 rows hit the PAN027 joint target hap2 chr3 contig and one hits the hap1 chr3 contig. The first two rows overlap each other by 1,998 bp on the query axis, explaining the 152,860 bp summed overlap versus 150,862 bp union coverage in the package summary.

## Comparison to prior sweepGA raw many:many -j 0

The prior whole-genome sweepGA/FastGA evidence review reported no chr3 raw `many:many -j 0` rows overlapping either autosomal chr9 candidate window:

- `segment_support.tsv:10`: PAN027 raw unchopped diagnostic has chr9 same-chromosome context only, 2 rows, 500,000 bp union coverage.
- `segment_support.tsv:14` and `:15`: PAN028 raw unchopped diagnostic has chr9 same-chromosome context plus one chr16 side fragment, but no chr3.

I also streamed the prior raw sweepGA PAFs at the manifest paths and got the same target-chromosome counts:

- PAN027 candidate: `{'chr9': 2}`
- PAN028 candidate: `{'chr9': 9, 'chr16': 1}`

Therefore, wfmash does change the raw evidence status for the autosomal Fig5 chr9 candidate windows: it recovers direct whole-genome chr3-target support that raw full-genome sweepGA/FastGA did not emit.

## Biological interpretation boundary

The technical result is positive for raw direct-alignment support: current upstream wfmash `v0.24.2`, run as whole-genome `-p 95`, recovers chr3-target homology overlapping both candidate windows before filtering.

This does not by itself prove the recombination mechanism or finalize Fig5 geometry. The rows are direct homology evidence suitable for future schematic consideration, but downstream biological interpretation should still remain separate from this provenance review and should account for wfmash's many-to-many raw mapping behavior, target haplotype multiplicity, same-window chr9 support, and any graph/untangle context.

## Review validation

Validation performed:

- Read the wfmash package README, configs, run script, command logs, binary manifest, job manifest, candidate-window summary, and prior sweepGA review/QA reports.
- Confirmed all three wfmash comparisons used full whole-genome FASTA inputs and the same comparison IDs as the corrected sweepGA package.
- Confirmed the primary runs used source-built upstream wfmash `v0.24.2-0-g774c01ff`, not stale Guix `0.12.5`.
- Confirmed exact primary commands used `-p 95` and `-B /dev/shm/.../tmp`.
- Recomputed SHA-256 checksums for the three primary raw PAFs.
- Independently streamed raw PAFs to reselect chr3-target rows overlapping the PAN027 and PAN028 candidate windows.
- Independently streamed prior sweepGA raw `many:many -j 0` PAFs to verify the absence of chr3 rows in those same windows.
- Checked `git status --short -- submission` and confirmed no `submission/` files were changed.
