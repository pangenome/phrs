# Independent validation of CHM13 copy-aware enrichment

## Decision

The copy-aware regional permutation implementation is **calibrated for the
synthetic nulls and stress conditions tested here**, after repairing one
integration-blocking identifier defect.  The calibration suite passed every
prespecified engineering gate.  This decision validates the engine mechanics;
it does not accept any biological enrichment finding.  A downstream final run
must still complete at least 99,999 permutations, satisfy its run-specific
placement and Monte Carlo gates, and pass `validate_engine_run.py` before any
result is reportable.

No manuscript file was read or edited for this task.

## Defect found and repaired

The frozen annotation product and engine were individually valid but not
directly interoperable:

- `outputs/copy_to_term.tsv.gz` keys 532,209 edges by coordinate-anchored
  `copy_id` (`CHM13v2.0|seqid:start-end|strand|GFF_ID`).
- `COPY_engine.py` requires a `locus_id` column matching
  `analysis_ready/chm13_gene_loci.tsv.gz` and correctly rejects unknown IDs.
- The two 61,312-member identifier sets are disjoint. Passing the frozen table
  directly would therefore reject all physical IDs, rather than silently join.

`prepare_engine_terms.py` now bridges the universes using the exact tuple
`chromosome,start1,end1,strand,gff_id`. Symbols, aliases, target status, and
annotation content cannot influence the join. The bridge requires a bijection,
rejects duplicate edges and missing metadata, splits GO by its frozen namespace,
and preserves both identifiers plus `frozen_source` in every output row.

The committed `engine_terms/MANIFEST.json` proves:

| Collection | Rows | Terms |
|---|---:|---:|
| GO BP | 125,938 | 11,184 |
| GO CC | 93,496 | 1,943 |
| GO MF | 78,117 | 4,847 |
| HGNC group | 33,836 | 2,111 |
| Reactome | 139,510 | 2,835 |
| Biotype | 61,312 | 26 |
| **Total** | **532,209** | — |

The total equals the frozen source-edge count exactly. A second build in a new
directory reproduced every collection SHA-256 and row count. This repair is an
adapter and frozen output addition; the statistical engine itself did not
require a code change.

## Independent implementation checks

`validate_engine_run.py` deliberately does not import `COPY_engine`. It parses
the saved placement TSVs, expands rigid components, and uses brute-force
half-open coordinate joins to recount midpoint and overlap assignments. For
every replicate and collection it compares:

- total physical-locus burden;
- annotated-locus denominator;
- every retained term's physical-copy count; and
- every retained term's chromosome-arm breadth.

It then independently recomputes plus-one empirical tails, complete
collection-by-assignment-by-statistic BH families, collection maxT, and primary
global maxT. It rejects duplicate term edges, unknown physical loci, missing
`copy_id`/`frozen_source` traces, impossible `term count > denominator > burden`
fields, and sensitivity-assignment leakage into global maxT. A corruption test
changes a cached term count and confirms that validation fails.

A real-input smoke run used all 61,312 loci, 37 target intervals/blocks, and the
biotype collection. The independent validator recounted 74 saved block rows,
114 result rows, all caches, all empirical p-values, BH, collection maxT, and
global maxT with exact agreement. The real primary candidate audit found:

- 37/37 observed starts reachable;
- minimum candidate count 250,000 (required minimum: 100);
- 100% width, arm, component, and terminal-stratum preservation; and
- acceptance 1.0 in the two-replicate integration smoke test.

This smoke run is explicitly a nonreportable mechanics check, not inference.

## Frozen term-map audit and hand calculations

The pre-existing independent `validate_outputs.py` and source verifier were
rerun. They validated 61,312 physical loci, 532,209 unique copy/source/term
edges, 53,349 metadata terms, 74,794 hierarchy edges, and every frozen source
checksum. No duplicate edge, missing metadata key, symbol-based mapping route,
pseudogene-parent propagation, or copy outside the universe was found.

Target leakage was inspected in `build_term_maps.py`: target membership is used
only to write/audit the target flag and hand-audit examples. It does not enter
stable-ID mapping, source annotation, term filtering, or the coordinate bridge.
All collection filters remain genome-wide and target-blind.

`hand_check_selected_terms.py` independently joined the 37 raw half-open target
intervals to locus coordinates. It reproduced 402 midpoint loci and 412
any-overlap loci, exactly matching the prepared flags. Selected direct counts
are frozen in `hand_check_selected_terms.tsv`; representative calculations are:

| Collection / term | Genome copies | Midpoint copies | Arms | Overlap copies |
|---|---:|---:|---:|---:|
| HGNC_GROUP:151, olfactory receptors family 4 | 129 | 19 | 7 | 19 |
| HGNC_GROUP:3661, tubulin beta family | 11 | 2 | 2 | 2 |
| HGNC_GROUP:521, PRD homeoboxes/pseudogenes | 177 | 67 | 3 | 67 |
| GO:0004984, olfactory receptor activity | 426 | 6 | 6 | 6 |
| BIOTYPE:pseudogene | 16,018 | 204 | 28 | 204 |
| BIOTYPE:protein_coding | 20,008 | 24 | 16 | 27 |

The TSV includes every contributing `locus_id` and coordinate `copy_id`, making
these counts auditable without symbol expansion. Repeated symbols remain
separate physical rows (the analysis universe contains 853 duplicated display
symbols, with as many as 76 physical copies for one symbol). Duplicate
`locus_id,term_id` edges are rejected rather than join-expanded.

## Null calibration design

Substantial simulation ran on Slurm, not the login node. Four deterministic
chunks used `SeedSequence(2026071391).spawn(4)` with spawn keys `[0]` through
`[3]`. Each chunk ran:

- 100,000 sampler draws for each exact scenario (400,000 combined): p-arm and
  q-arm boundaries, a rigid two-component/gapped block, mask-fraction matching,
  and two-block collision rejection;
- 250 pseudo-observed regional sets against 9,999 independent reference sets
  using permutation-rank reuse (1,000 pseudo-observations combined);
- synthetic singleton loci, tandem copies, duplicated local blocks, correlated
  terms, arm effects, heavy copy multiplicity, and terminal gene-density
  gradients;
- all-locus burden, ten-term copy burden, composition, breadth, BH, and global
  maxT calibration;
- primary, terminal cross-arm, and adjacent constraint checks; and
- a 10% planted-term engineering control plus the historical instance-expanded
  hypergeometric positive failure control.

The tolerance was fixed in code and output before inspecting results:

- exact sampling: chi-square p >= 0.001, total-variation distance <= 0.01,
  every state reachable;
- Type I at alpha 0.10, 0.05, and 0.01: 95% exact-binomial lower limit <= alpha
  and point estimate <= `max(alpha + 0.02, 1.5 * alpha)`;
- BH complete-null rejection at q=0.05 and global maxT family rejection at
  0.05: point <= 0.07 and lower 95% limit <= 0.05;
- composition zero denominator rate <= 1%; and
- zero geometry/arm/stratum/annulus violations, with cross-arm sampling actually
  exercised.

## Calibration results

All gates passed. Exact sampler results were:

| Scenario | Total variation | Chi-square p | Reachable | Pass |
|---|---:|---:|---|---|
| p boundary | 0.003385 | 0.1502 | all | yes |
| q boundary | 0.002909 | 0.1460 | all | yes |
| rigid gap | 0.001985 | 0.8391 | all | yes |
| mask fraction | 0.002923 | 0.07188 | all | yes |
| collision joint state | 0.009271 | 0.1048 | all | yes |

Regional complete-null rejection rates were:

| Family | alpha=0.10 | alpha=0.05 | alpha=0.01 |
|---|---:|---:|---:|
| all-locus burden | 0.091 | 0.059 | 0.014 |
| term copy burden | 0.0648 | 0.0286 | 0.0068 |
| composition | 0.0828 | 0.0424 | 0.0072 |
| breadth | 0.0295 | 0.0112 | 0.0009 |
| any BH rejection | 0.089 | 0.048 | 0.009 |
| any global maxT rejection | 0.101 | 0.021 | 0.000 |

At alpha 0.05, the burden rate 0.059 had a 95% interval of 0.0452–0.0754 and
met the stated point ceiling 0.075. BH rejection was 48/1,000 (0.048; lower
limit 0.0356), and global maxT rejection was 21/1,000 (0.021; lower limit
0.0130). Composition had zero denominators in 0/40,996 generated sets. There
were zero primary/terminal/adjacent constraint violations, and every chunk
exercised cross-arm placements.

The planted term produced raw power 0.952 and BH power 0.786; mean false
discovery proportion among null terms was 0.0261. The deliberately mismatched
instance-expanded hypergeometric test rejected 10,536/40,000 clustered-null
draws (0.2634; 95% interval 0.2591–0.2677). Its lower limit exceeds 0.05, so the
required positive failure control passed and demonstrates that the suite can
detect the historical anti-conservative method.

Full machine-readable values are in
`calibration_results/job-20260713T162803Z/combined.json` and its four chunk
files.

## Slurm provenance

Before submission, `sinfo` and `scontrol show partition workers` were inspected.
The `workers` partition was UP with seven nodes, 336 CPUs, no configured time
limit, and no partition-level memory cap; the scripts nevertheless request
bounded resources.

| Job | Request | Host | State / exit | Elapsed |
|---|---|---|---|---:|
| 1754105 array tasks 0–3 | each 1 CPU, 4 GiB, 30 min | octopus08 | COMPLETED / 0:0 | 20–21 s |
| 1754106 dependent merge | 1 CPU, 1 GiB, 5 min | octopus08 | COMPLETED / 0:0 | 1 s |

The merge was submitted with `--dependency=afterok:1754105`, so it could not
validate partial/failed chunks. Slurm accounting did not populate MaxRSS for
these short jobs; requested resources, host, parameters, seeds, output paths,
shell timing, and successful exits are retained in the committed stdout/stderr
logs. Chunks are combined by summing rejection/trial and exact-state counts;
rates, exact binomial intervals, chi-square p-values, and pass/fail gates are
recomputed from the aggregate rather than averaging chunk rates.

## Adversarial checklist

| Risk | Result |
|---|---|
| Silent symbol deduplication | absent; physical IDs are primary, duplicate-symbol tests and hand loci pass |
| Term-map join inflation | absent; duplicate edges rejected, source total equals collection total |
| Target/background leakage | absent in mapping/filter code; target used only for audits |
| Multi-copy / duplicated blocks | retained and calibrated; weighted failure control detects fragmentation error |
| Interval boundaries / masks / rigid gaps | exact scenarios pass |
| Arm and sensitivity restrictions | zero violations; primary, terminal, adjacent exercised separately |
| Impossible contingency tables | validator rejects them; none in real smoke output |
| Multiple-testing family leakage | independent BH/maxT recomputation agrees; overlap excluded from global maxT |
| Output overwrite | existing-output refusal and immutable resume tests pass |
| Seed/resume reproducibility | repeated seed test, byte-exact resume test, spawn-key chunks, repeated term build pass |
| Result-to-source traceability | every engine term edge contains locus_id, copy_id, frozen_source; checksums frozen |

## Commands and tests

The clean validation invocation is:

```bash
cd paper_prep/_brainstorming/chm13_copy_enrichment
python3 fetch_sources.py
python3 validate_outputs.py
python3 prepare_engine_terms.py
python3 hand_check_selected_terms.py
python3 -m py_compile COPY_engine.py prepare_engine_terms.py \
  validate_engine_run.py calibration_suite.py combine_calibration.py
guix shell python python-numpy -- python3 -m unittest -v \
  test_COPY_engine.py test_build_inputs.py tests/test_build_term_maps.py \
  test_prepare_engine_terms.py test_validate_engine_run.py \
  test_calibration_suite.py
```

The Slurm commands and all parameters are encoded in
`slurm_calibration.sbatch` and `slurm_calibration_combine.sbatch`. Unit tests
cover coordinate mismatch, duplicate-edge rejection, duplicated-symbol copy
retention, black-box recount, cached-array corruption, deterministic chunk
seeds, and the required weighted-method failure control, in addition to the
upstream engine and term-map suites.
