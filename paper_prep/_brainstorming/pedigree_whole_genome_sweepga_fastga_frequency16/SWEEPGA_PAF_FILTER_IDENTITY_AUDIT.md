# sweepGA PAF filter identity scoring audit

Date: 2026-06-22

## Decision

The downstream Fig5 f16 validated chop rerun **may proceed**, but only with the
repaired `pafchop-rs` output and a sweepGA filter command that explicitly uses
`--scoring ani` for identity-first overlap ranking.

Do **not** use the existing f16 `scripts/filter_paf.py` command as-is for
identity-sensitive 1:1 / N:many filtering, because it currently omits
`--scoring ani` and therefore uses sweepGA's default `log-length-ani` score.
That default can prefer a longer, lower-identity chunk over a shorter,
higher-identity chunk.

For the f16 chopped rerun, the safe per-comparison command line is:

```bash
PACKAGE_DIR=/moosefs/erikg/phrs/paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16
CID=PAN027mat_vs_PAN010_joint
FILTER_ID=four_many_chopped
NUM_MAPPINGS=4:many
SCAFFOLD_JUMP=0
SCRATCH_BASE=/dev/shm
WORK=$(mktemp -d "$SCRATCH_BASE/sweepga_filter.${CID}.XXXXXX")
gzip -dc "$PACKAGE_DIR/chopped_paf_l10000_o0/${CID}.chopped_l10000_o0.paf.gz" > "$WORK/input.paf"
/home/erikg/.cargo/bin/sweepga "$WORK/input.paf" \
  --num-mappings "$NUM_MAPPINGS" \
  --scaffold-jump "$SCAFFOLD_JUMP" \
  --scoring ani \
  --output-file "$WORK/filtered.paf"
gzip -c "$WORK/filtered.paf" > "$PACKAGE_DIR/filtered_paf/${CID}.${FILTER_ID}.paf.gz"
gzip -t "$PACKAGE_DIR/filtered_paf/${CID}.${FILTER_ID}.paf.gz"
rm -rf "$WORK"
```

Repeat with the appropriate `CID`, `FILTER_ID`, and `NUM_MAPPINGS` from
`config/comparisons.tsv` and `config/filter_matrix.tsv`. The `many:many`
chopped layer is a copy-through control and does not need sweepGA filtering.
The current f16 matrix contains `four_many_chopped` with `NUM_MAPPINGS=4:many`;
the same safety rule applies to any future `1:1` row added to the matrix.

## Binary and source inspected

Installed executable:

- Path: `/home/erikg/.cargo/bin/sweepga`
- `ls -l`: executable, 6,700,744 bytes, timestamp `2026-06-21 09:28`
- Version: `sweepga 0.1.1`
- SHA256: `a0d7ac0c3312080d67de96d85cdcad9ce0c5a7e523897109b7f598c186ab85a6`

The installed help reports:

- `--num-mappings` default: `many:many`
- `--scoring` default: `log-length-ani`
- accepted scoring values: `ani`, `length`, `length-ani`, `log-length-ani`,
  `matches`
- `--scaffold-jump 0` disables scaffolding

Local source was available at `/home/erikg/sweepga`:

- Git commit: `d4b551d6e1fc05e3c5027371afb74481b67a1d9e`
- Worktree status at audit time: dirty (`D .cargo/config.toml`,
  `?? .cargo/config.toml.bak`, `?? build.sh`)

Relevant source observations:

- `src/main.rs` maps `--scoring ani` to `ScoringFunction::Identity` and the
  default `log-length-ani` to `ScoringFunction::LogLengthIdentity`.
- `src/paf_filter.rs` parses PAF columns 10/11, then updates identity from
  `dv:f` and `cg:Z` tags when present.
- `src/paf_filter.rs` exits after row-level plane sweep when
  `scaffold_gap == 0`; scaffold-chain merging is only reached when
  `--scaffold-jump` is nonzero.
- `src/plane_sweep_exact.rs` defines `ani` scoring as pure row identity and
  `log-length-ani` as `identity * ln(query_span)`.

The synthetic tests below match those source observations on the installed
binary.

## Synthetic fixtures

Fixtures were written under:

`paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/synthetic_paf_filter_audit/`

The harness is:

`paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/scripts/run_sweepga_paf_filter_identity_audit.sh`

The generated summary is:

`paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/summaries/sweepga_paf_filter_identity_audit.tsv`

Each PAF fixture is tiny and synthetic. No whole-genome data were processed by
the audit. The fixtures include:

- equal-length competitors with different identity;
- unequal-length competitors where default `log-length-ani` and `--scoring ani`
  choose different winners;
- a case where PAF columns 10/11 disagree with `cg:Z`;
- a case where `dv:f` disagrees with both PAF columns and `cg:Z`;
- non-overlapping rows to distinguish plane-sweep behavior from a simple
  per-name collapse;
- an explicit scaffold-chain probe with and without nonzero scaffold jump.

The fixture PAF records include recomputed PAF columns 10/11 and `cg:Z` tags in
the form sweepGA expects. Two precedence probes intentionally contain
disagreement among columns/tags to determine which fields control scoring.

## Results

All rows in `summaries/sweepga_paf_filter_identity_audit.tsv` pass after
encoding the observed semantics as expectations.

Key outcomes:

- With `--num-mappings 1:1 --scaffold-jump 0` and default scoring,
  `log-length-ani` can keep the 1000 bp / 80% identity row instead of the
  100 bp / 99% identity row.
- With `--num-mappings 1:1 --scaffold-jump 0 --scoring ani`, the 100 bp / 99%
  identity row wins, proving that `--scoring ani` ranks by per-row identity in
  the installed binary.
- With `--num-mappings many:many --scaffold-jump 0`, both overlapping
  competitors are preserved; `--scoring ani` has no effect when the mapping
  axis is unbounded.
- `1:1` is not a simple "one PAF row per query name and target name" collapse.
  It is plane-sweep filtering on query and target axes. Non-overlapping rows
  for the same query/target names can both pass.
- With `--scaffold-jump 0`, no scaffold-chain merge occurs before output.
  Filtering is per input PAF row/chunk, subject to the row-level plane sweep.
- With explicit nonzero scaffold jump and scaffold mass, nearby rows can be
  merged/evaluated as scaffold chains and emitted as chain member rows.
- When `cg:Z` is present, it affects scoring identity. In the
  `columns_vs_cg_precedence` fixture, the row with lower PAF column identity
  but higher `cg:Z` identity wins under `--scoring ani`.
- `dv:f` can override CIGAR-derived identity in the current parser order. In
  the `dv_overrides_columns` fixture, the row with lower PAF column identity
  and lower `cg:Z` identity but `dv:f:0.010000` wins under `--scoring ani`.
  Therefore stale `dv:f` tags are dangerous for identity-sensitive reruns.

## Field precedence and safety implications

For validated f16 chopped PAFs, the safe input condition is:

- PAF columns 10/11 are recomputed per chunk;
- `cg:Z` is recomputed per chunk;
- identity-sensitive optional tags such as `dv:f`, `de:f`, `NM:i`, and `df:i`
  are either recomputed per chunk or dropped.

This condition is met by the repaired `pafchop-rs` described in
`paper_prep/_brainstorming/pafchop-rs/PAF_SEMANTICS_VALIDATION.md`, not by the
old f16 chopped outputs that linearly interpolated columns 10/11.

Because sweepGA preserves original PAF rows in output, the output columns remain
whatever the input row supplied. The filtering decision, however, uses parsed
metadata from those rows/tags. That is why the rerun must start from repaired
chopped PAFs and must not trust stale identity-sensitive tags.

## Current wrapper status

`scripts/filter_paf.py` currently constructs:

```bash
/home/erikg/.cargo/bin/sweepga --num-mappings <num_mappings> --scaffold-jump <scaffold_jump> --output-file <filtered_tmp> <tmp_paf>
```

That command is not safe for the identity-first f16 filter rerun because it
omits `--scoring ani`. Use the explicit command above, or update the wrapper to
pass `--scoring ani` before running the f16 matrix.

## Validation command

The audit harness was run locally on only tiny fixtures:

```bash
paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/scripts/run_sweepga_paf_filter_identity_audit.sh
```

Validation check:

```bash
awk -F '\t' 'NR==1 || $6 != "PASS" {print}' \
  paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/summaries/sweepga_paf_filter_identity_audit.tsv
```

Only the header was printed, confirming no failing synthetic tests.
