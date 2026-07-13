# Reproducing the CHM13 physical-copy enrichment v1 package

This is the compact operational guide for the audited v1 run. For scientific
interpretation and author decisions, read [FINAL_REPORT.md](FINAL_REPORT.md).
For the normative statistical contract, read
[STATISTICAL_SPEC.md](STATISTICAL_SPEC.md).

## Frozen run identity

The run uses the six collections and family map in
[final.args](results/v1/config/final.args), with midpoint and overlap assignment,
99,999 permutations per background, a batch size of 1,000, and seeds 2026071301,
2026071302, and 2026071303 for primary, terminal, and adjacent modes,
respectively ([PREFLIGHT.json](results/v1/config/PREFLIGHT.json)). The source run
commit recorded by the completed manifests is
`58daa8dc8e37c5be2bc3e35aa3b26c0429c82496`
([RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json)).

Requirements:

- repository access to the committed `data/` and analysis inputs;
- moosefs worker environment with `guix`, Python, NumPy, and Slurm `sbatch`;
- an empty `results/v1/work/{primary,terminal,adjacent}` output state (the engine
  refuses unsafe overwrite); and
- about 18.08 GiB for the transient batch products represented by the frozen
  digest inventory
  ([RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json)).

Do not run the full permutation workload on a login node.

## Verify frozen inputs before submission

From the repository root:

```bash
cd paper_prep/_brainstorming/chm13_copy_enrichment

python3 fetch_sources.py
python3 validate_outputs.py
(cd analysis_ready && \
  awk 'NR > 1 {print $1 "  " $4}' MANIFEST.sha256 | sha256sum -c -)

python3 - <<'PY'
import hashlib
import json
from pathlib import Path

root = Path.cwd()
preflight = json.loads((root / 'results/v1/config/PREFLIGHT.json').read_text())
assert preflight['pass'] is True

def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

for record in preflight['term_maps']['observed'].values():
    # The collection name is encoded by the matching manifest record below.
    assert record['rows'] > 0

term_manifest = json.loads((root / 'engine_terms/MANIFEST.json').read_text())
for record in term_manifest['collections'].values():
    assert sha256(root / record['path']) == record['sha256']

assert sha256(root / preflight['arguments']['path']) == preflight['arguments']['sha256']
assert sha256(root / 'results/v1/config/family_map.tsv.gz') == \
    preflight['family_diagnostic_map']['sha256']
assert sha256(root / 'calibration_results/job-20260713T162803Z/combined.json') == \
    preflight['calibration']['sha256']
print('frozen input checks passed')
PY

guix shell python python-numpy -- python3 -m unittest -v \
  test_COPY_engine.py \
  test_build_inputs.py \
  tests/test_build_term_maps.py \
  test_prepare_engine_terms.py \
  test_validate_engine_run.py \
  test_calibration_suite.py
```

Expected preflight product:
[results/v1/config/PREFLIGHT.json](results/v1/config/PREFLIGHT.json) with
`"pass": true`. Its frozen SHA-256 and all expected input/term-map counts are
inside that JSON. Frozen source objects are inventoried in
[sources/SOURCE_MANIFEST.tsv](sources/SOURCE_MANIFEST.tsv), prepared inputs in
[analysis_ready/MANIFEST.sha256](analysis_ready/MANIFEST.sha256), and engine term
maps in [engine_terms/MANIFEST.json](engine_terms/MANIFEST.json).

## Exact clean rerun command

The reference engine writes the complete placement and inferential products
before entering its non-scaling reference leave-one scan. The audited v1
workflow therefore submits the three modes, waits until all three manifests and
term tables are complete, cancels the reference diagnostic stage, runs the
byte-validated cached-subtraction finalizer, and then aggregates. Run this exact
sequence from the repository root only after confirming that no prior
`results/v1/work/<mode>/run_manifest.json` exists.

```bash
set -euo pipefail

ROOT=$(git rev-parse --show-toplevel)
ANALYSIS_DIR="$ROOT/paper_prep/_brainstorming/chm13_copy_enrichment"
RUN_ROOT="$ANALYSIS_DIR/results/v1/work"

ARRAY_JOB=$(sbatch --parsable --chdir="$ROOT" \
  --export=ALL,REPO_ROOT="$ROOT",RUN_ROOT="$RUN_ROOT",PERMUTATIONS=99999,BATCH_SIZE=1000 \
  paper_prep/_brainstorming/chm13_copy_enrichment/results/v1/slurm_final_array.sbatch)

while ! python3 - "$RUN_ROOT" <<'PY'
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
for mode in ('primary', 'terminal', 'adjacent'):
    manifest_path = root / mode / 'run_manifest.json'
    terms_path = root / mode / 'term_results.tsv'
    if not manifest_path.is_file() or not terms_path.is_file():
        raise SystemExit(1)
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('completed_permutations') != 99999:
        raise SystemExit(1)
raise SystemExit(0)
PY
do
  sleep 30
done

# Only the completed non-scaling reference diagnostic stage is cancelled.
scancel "$ARRAY_JOB"
while squeue -h -j "$ARRAY_JOB" | grep -q .; do
  sleep 5
done

FINALIZE_JOB=$(sbatch --parsable --chdir="$ROOT" \
  --export=ALL,REPO_ROOT="$ROOT" \
  paper_prep/_brainstorming/chm13_copy_enrichment/results/v1/slurm_fast_finalize.sbatch)

AGGREGATE_JOB=$(sbatch --parsable --chdir="$ROOT" \
  --dependency="afterok:$FINALIZE_JOB" \
  --export=ALL,REPO_ROOT="$ROOT" \
  paper_prep/_brainstorming/chm13_copy_enrichment/results/v1/slurm_aggregate.sbatch)

printf 'array_job=%s\nfinalize_job=%s\naggregate_job=%s\n' \
  "$ARRAY_JOB" "$FINALIZE_JOB" "$AGGREGATE_JOB"
```

The generation array is expected to finish its 99,999 placements and complete
term tables before being recorded as cancelled during the reference diagnostic;
that scheduler state alone is not a failed inference run. Do not cancel it until
the loop above succeeds for all three modes. The finalizer and aggregate jobs,
in contrast, must complete with exit code zero
([SLURM_COMPLETION.json](results/v1/SLURM_COMPLETION.json)).

The historically executed submissions, including the retained retry sequence,
are frozen verbatim in
[SUBMISSION_COMMANDS.sh](results/v1/SUBMISSION_COMMANDS.sh). Scheduler states,
exit codes, retry reasons, requests, and job IDs are in
[SLURM_COMPLETION.json](results/v1/SLURM_COMPLETION.json),
[job_manifest.tsv](results/v1/job_manifest.tsv), and
[slurm_accounting.tsv](results/v1/slurm_accounting.tsv).

## Expected outputs

Each of `results/v1/work/primary`, `results/v1/work/terminal`, and
`results/v1/work/adjacent` must contain:

- `run_manifest.json`, `placement_blocks.tsv`, `candidate_spaces.tsv`, and
  `placement_qc.tsv`;
- `batches/placements.*.tsv.gz` and the assignment/collection statistic arrays;
- `burden_results.tsv`, `term_filtering.tsv`, and `term_results.tsv`;
- `driver_summary.tsv`, `driver_group_counts.tsv`, and
  `leave_one_sensitivity.tsv`.

Successful aggregation must produce or reproduce:

- [RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json) with `"pass": true`, three
  completed modes, and 99,999 effective permutations per mode;
- [all_burden_results.tsv](results/v1/final_tables/all_burden_results.tsv);
- [all_term_results.tsv.gz](results/v1/final_tables/all_term_results.tsv.gz);
- [all_term_filtering.tsv.gz](results/v1/final_tables/all_term_filtering.tsv.gz);
- [primary_evidence.tsv](results/v1/final_tables/primary_evidence.tsv);
- [composition_cross_null.tsv.gz](results/v1/final_tables/composition_cross_null.tsv.gz);
- [multiplicity_families.tsv](results/v1/final_tables/multiplicity_families.tsv);
- [primary_driver_summary.tsv](results/v1/final_tables/primary_driver_summary.tsv),
  [primary_driver_group_counts.tsv](results/v1/final_tables/primary_driver_group_counts.tsv),
  and [primary_leave_one_sensitivity.tsv](results/v1/final_tables/primary_leave_one_sensitivity.tsv);
- [COMPARATOR_deduplicated_symbol_ORA_not_copy_aware.tsv.gz](results/v1/final_tables/COMPARATOR_deduplicated_symbol_ORA_not_copy_aware.tsv.gz);
- [null_sensitivity_composition.svg](results/v1/plots/null_sensitivity_composition.svg);
- `transient_file_checksums.tsv.gz` for every large uncommitted batch product;
  and
- [SOURCE_OUTPUT_CHECKSUMS.tsv](results/v1/SOURCE_OUTPUT_CHECKSUMS.tsv) for the
  committed source-sized products.

For an internal byte-level verification of the source-sized package generated
by aggregation:

```bash
cd paper_prep/_brainstorming/chm13_copy_enrichment/results/v1
sha256sum -c <(awk 'NR>1 {print $3 "  " $1}' SOURCE_OUTPUT_CHECKSUMS.tsv)
```

The transient `work/` outputs are intentionally not committed. Their expected
11,730 paths, sizes, and SHA-256 values are frozen in
[transient_file_checksums.tsv.gz](results/v1/transient_file_checksums.tsv.gz),
with the aggregate count in
[RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json).

## Acceptance checks

Before interpreting a rerun, confirm all of the following:

```bash
python3 - <<'PY'
import json
from pathlib import Path

root = Path('paper_prep/_brainstorming/chm13_copy_enrichment/results/v1')
audit = json.loads((root / 'RESULT_AUDIT.json').read_text())
assert audit['pass'] is True
assert audit['effective_permutations'] == {
    'primary': 99999,
    'terminal': 99999,
    'adjacent': 99999,
}
assert audit['total_valid_joint_permutations'] == 299997
assert audit['output_counts']['primary_evidence_rows'] == 41
assert audit['primary_evidence_requiring_extension'] == 41
assert audit['primary_evidence_unique_terms'] == 36
assert audit['single_driver_sensitive_primary_evidence_terms'] == 36
assert audit['null_sensitive_primary_composition_evidence_terms'] == 8
print('v1 audit checks passed')
PY
```

These expected counts come directly from the frozen machine audit
([RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json)). A count mismatch is a rerun
or environment discrepancy, not a reason to edit the report manually.

The finalizer must also retain byte-exact agreement with its reference smoke
implementation:

```bash
python3 - <<'PY'
import json
from pathlib import Path

p = Path('paper_prep/_brainstorming/chm13_copy_enrichment/results/v1/FINALIZER_EQUIVALENCE.json')
r = json.loads(p.read_text())
assert r['pass'] is True
print('finalizer equivalence passed')
PY
```

See [FINALIZER_EQUIVALENCE.json](results/v1/FINALIZER_EQUIVALENCE.json) for the
machine-readable comparison.

## Interpreting the package

The provisional term screen is primary midpoint, within-collection BH q <= 0.05
and global maxT p <= 0.05. Do not treat a raw p-value, a sensitivity q-value, or
the symbol-level comparator as a discovery route. Every current provisional row
requires a stream-preserving extension, and every unique provisional term is
copy/family-sensitive
([RESULT_AUDIT.json](results/v1/RESULT_AUDIT.json)).

No command in this guide reads or modifies `submission/paper.tex`.
