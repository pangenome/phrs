# CHM13 physical-copy enrichment engine

`COPY_engine.py` implements the coordinate-randomization contract in
`STATISTICAL_SPEC.md`. It treats every CHM13 GFF3 `gene` row as a separate
physical locus, moves coherent PHR interval blocks, and never expands or
deduplicates gene symbols. This is analysis code; it does not edit or read the
manuscript.

## Reproducible environment and tests

The engine requires Python 3 and NumPy with `PCG64DXSM` support. On the moosefs
deployment environment:

```bash
cd paper_prep/_brainstorming/chm13_copy_enrichment
guix shell python python-numpy -- \
  python3 -m unittest -v test_COPY_engine.py test_build_inputs.py
```

The tests enumerate small null spaces and cover physical copy multiplicity,
duplicate term-map edges, half-open boundaries, midpoint versus overlap
assignment, masks, rigid blocks, empty collections, zero counts, degeneracy,
Clopper--Pearson intervals, collection-scoped BH, exact maxT, sampler
invariants, custom collection identity, safe output refusal, paired leave-one
diagnostics, and deterministic resume.

## Frozen term and group inputs

No functional database is silently downloaded. Every term collection is a
target-blind, frozen tab-separated file:

```text
locus_id	term_id	term_name
LOC124900618	GO:0000001	mitochondrion inheritance
```

`term_name` is optional; `locus_id` and `term_id` are required. Duplicate
`locus_id`/`term_id` edges are collapsed, but two different locus IDs with the
same gene symbol remain two copies. An unknown locus ID is a hard error.
Filtering is genome-wide and target-blind (`--min-term-loci 5` and
`--min-term-arms 2` by default). Collections remain separate even when they
contain the same term ID.

Optional family and sequence-identity maps use these schemas:

```text
locus_id	family_id
LOC124900618	FAMILY_001
```

```text
locus_id	identity_id
LOC124900618	IDCLUSTER_001
```

Missing group assignments become locus-specific singleton groups. They are
never combined into one missing-value group.

The committed annotation build uses a stricter coordinate-anchored `copy_id`,
whereas the engine's analysis table uses the GFF `locus_id`.  Generate (or
verify) the lossless one-to-one bridge before a run:

```bash
python3 prepare_engine_terms.py
```

The six deterministic files under `engine_terms/` retain `locus_id`, `copy_id`,
and `frozen_source`. `engine_terms/MANIFEST.json` records the exact coordinate
join, source/output checksums, row counts, and GO namespace split.  Use these
files as named collections; do not pass `outputs/copy_to_term.tsv.gz` directly.

## Primary run

Use a new output directory and provide the prespecified seed explicitly:

```bash
guix shell python python-numpy -- python3 COPY_engine.py \
  --terms GO_BP=/absolute/frozen/go_bp.tsv.gz \
  --terms GO_MF=/absolute/frozen/go_mf.tsv.gz \
  --terms GO_CC=/absolute/frozen/go_cc.tsv.gz \
  --family-map /absolute/frozen/families.tsv.gz \
  --identity-map /absolute/frozen/identity_clusters.tsv.gz \
  --mode primary \
  --seed 2026071301 \
  --permutations 99999 \
  --batch-size 1000 \
  --output runs/primary_2026071301
```

The default prepared inputs are the `analysis_ready/` CHM13 arm, 37-interval,
and 61,312-locus tables. There is no hardcoded target interval or locus count in
the engine: the manifest records the values read from the frozen inputs.
Midpoint and any-overlap statistics are calculated from the same placements by
default. Use a single `--assignment midpoint` only when intentionally omitting
the boundary sensitivity.

Runs below 99,999 permutations are rejected unless `--allow-pilot` is supplied.
Every pilot row is marked `pilot_nonreportable`; the flag is intended for exact
tests and runtime planning, not inference.

## Resume and extension

Batch files are committed atomically. The manifest saves the complete
`PCG64DXSM` state after each completed batch. To append unused draws without
changing any prior placement:

```bash
guix shell python python-numpy -- python3 COPY_engine.py \
  --terms GO_BP=/absolute/frozen/go_bp.tsv.gz \
  --terms GO_MF=/absolute/frozen/go_mf.tsv.gz \
  --terms GO_CC=/absolute/frozen/go_cc.tsv.gz \
  --family-map /absolute/frozen/families.tsv.gz \
  --identity-map /absolute/frozen/identity_clusters.tsv.gz \
  --mode primary --seed 2026071301 --batch-size 1000 \
  --permutations 999999 --output runs/primary_2026071301 --resume
```

All immutable options, input paths, bytes, and SHA-256 checksums must match the
original run. `--resume` refuses changed collections, thresholds, assignments,
seeds, or batch sizes. Without `--resume`, an existing output directory is a
hard error, preventing custom collections from overwriting an earlier result.
Rows marked `extension_required=1` require review and normally extension to
999,999 draws under the statistical contract.

## Prespecified background sensitivities

Run each background in a distinct directory and with its own seed. They are
never pooled with the primary null:

```bash
# Terminal-matched cross-arm sensitivity
guix shell python python-numpy -- python3 COPY_engine.py \
  --terms GO_BP=/absolute/frozen/go_bp.tsv.gz \
  --mode terminal --seed 2026071302 --permutations 99999 \
  --output runs/terminal_2026071302

# Same-arm proximal 5-Mb annulus sensitivity
guix shell python python-numpy -- python3 COPY_engine.py \
  --terms GO_BP=/absolute/frozen/go_bp.tsv.gz \
  --mode adjacent --seed 2026071303 --permutations 99999 \
  --output runs/adjacent_2026071303
```

The adjacent run also writes the deterministic immediately proximal comparator,
explicitly labeled as having no p-value. A non-empty `--mask BED` invokes exact
integer-start enumeration and enforces the one-percentage-point excluded-base
fraction match. Very large masked spaces stop with instructions to precompute
the exact candidates instead of falling back to biased rejection.

## Output contract

Each output directory is self-contained and schema-versioned:

- `run_manifest.json`: checksums, options, command, Git/Python/NumPy versions,
  target counts, batch boundaries, RNG state, and elapsed time;
- `placement_blocks.tsv`, `candidate_spaces.tsv`, and `placement_qc.tsv`: rigid
  geometry, exact candidate counts, observed-start reachability, and joint-draw
  acceptance;
- `batches/placements.*.tsv.gz`: every randomized block coordinate, sufficient
  for independent recounting and paired diagnostic removal;
- `batches/*.npy`: bounded, memory-mappable per-replicate burden, term-copy,
  collection-denominator, and breadth arrays for midpoint and overlap;
- `burden_results.tsv`: all-physical-locus burden, null interval, smoothed ratio,
  plus-one p-value, exceedances, and exact Monte Carlo interval;
- `term_filtering.tsv`: all frozen terms and target-blind retention reasons;
- `term_results.tsv`: copy burden, composition, and arm breadth, including
  within-collection BH/BY, conservative BH on Monte Carlo upper bounds,
  collection maxT, primary global maxT, effect sizes, zero-denominator rate,
  degeneracy, and boundary sensitivity;
- `driver_summary.tsv` and `driver_group_counts.tsv`: locus, symbol, arm, family,
  and identity driver concentrations without symbol deduplication; and
- `leave_one_sensitivity.tsv`: exhaustive or 80%-coverage family/identity
  removals recomputed against the identical saved placements.

BH is applied separately to each collection-by-statistic family, including all
testable terms. Global maxT includes all midpoint copy-burden, composition, and
breadth hypotheses across collections for the primary mode only. Sensitivities
receive their own collection families and never acquire a primary global maxT.

## Benchmark

Reproduce the committed GO-scale benchmark with:

```bash
guix shell python python-numpy -- python3 benchmark_COPY_engine.py \
  --terms 15000 --annotations-per-locus 12 --permutations 500 \
  --output benchmark_results.json
```

See `BENCHMARK_REPORT.md` for measured throughput, memory, disk projections,
the selected final permutation count, and the validated Slurm result in
`benchmark_results_slurm.json`.

## Slurm execution

Do not run final permutations or substantial calibration/benchmark workloads on
the login node. Create the log directory (it is already tracked), then submit a
cluster benchmark from this directory:

```bash
ROOT=$(git rev-parse --show-toplevel)
sbatch --export=ALL,REPO_ROOT="$ROOT",BENCH_PERMUTATIONS=5000,\
BENCH_TERMS=15000,ANNOTATIONS_PER_LOCUS=12,MASTER_SEED=2026071301,\
BENCH_OUTPUT="$PWD/benchmark_results.slurm.json" slurm_benchmark.sbatch
```

For a final run, construct a frozen argument file with exactly one CLI argument
per line. It can contain any number of collection pairs and optional maps:

```text
--terms
GO_BP=/absolute/frozen/go_bp.tsv.gz
--terms
GO_MF=/absolute/frozen/go_mf.tsv.gz
--family-map
/absolute/frozen/families.tsv.gz
--identity-map
/absolute/frozen/identity_clusters.tsv.gz
```

Submit the primary run (change mode/seed/output for each sensitivity):

```bash
ROOT=$(git rev-parse --show-toplevel)
sbatch --export=ALL,REPO_ROOT="$ROOT",MODE=primary,MASTER_SEED=2026071301,\
PERMUTATIONS=99999,BATCH_SIZE=1000,OUTPUT_DIR=/absolute/runs/primary,\
ARGS_FILE=/absolute/frozen/primary.args slurm_final_run.sbatch
```

Both scripts request one CPU because the reference engine is deliberately
single-threaded, and they record the Slurm job ID, host, seed, permutation count,
batch size, shell timing, stdout, and stderr. Inspect a job
with `squeue -j JOBID`; after completion, retain accounting with:

```bash
sacct -j JOBID --format=JobID,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
```

## Independent validation

`validate_engine_run.py` is a black-box validator: it does not import the
engine.  It expands every saved placement, brute-force joins physical loci,
checks all cached arrays, and independently recalculates empirical tails, BH,
collection maxT, and global maxT.  Frozen trace columns are required by default:

```bash
guix shell python python-numpy -- \
  python3 validate_engine_run.py /absolute/completed/run \
  --output /absolute/completed/run/independent_validation.json
```

The chunked synthetic calibration is submitted with
`slurm_calibration.sbatch`; `slurm_calibration_combine.sbatch` must be submitted
with an `afterok` dependency on the array.  It evaluates the prespecified exact
sampler, Type-I, BH, maxT, zero-denominator, constraint, planted-power, and
historical weighted-hypergeometric failure-control gates.  See
`VALIDATION_REPORT.md` and `calibration_results/` for the validated invocation.
