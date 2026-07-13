# CHM13 copy-enrichment engine benchmark

## Decision

Use **99,999 valid joint permutations** for each final primary, terminal, and
adjacent run, in batches of 1,000. This meets the statistical contract and is
feasible on the moosefs analysis host. Extend an unchanged run to 999,999 when
the engine flags low exceedance counts, a near-threshold result, or an
observed-positive/null-always-zero statistic. A 19,999-draw pilot is useful for
planning but is not reportable inference.

## Measured workload

Benchmark date: 2026-07-13 UTC. The committed `benchmark_results.json` was
generated with NumPy 1.23.2 and the command in `README.md`. It used the complete
frozen CHM13 universe and target geometry:

- 61,312 physical CHM13 gene loci;
- 37 PHR interval templates and 37 placement blocks;
- 15,000 synthetic GO-scale retained terms;
- 12 correlated annotations per physical locus (735,744 locus-term edges);
- both midpoint and any-overlap assignment; and
- 500 primary same-arm joint permutations.

The synthetic map is deterministic and deliberately gives nearby/copy-related
loci overlapping term memberships. It measures the coordinate recount and
many-to-many term burden rather than an unrealistically sparse one-term map.

## Results

| Measurement | Observed |
|---|---:|
| Input and synthetic-map setup | 2.886 s |
| 500 joint placements, both assignments | 6.526 s |
| Throughput | 76.61 joint permutations/s |
| MaxT kernel, all 15,000 terms | 0.588 s |
| Peak resident memory | 173.1 MiB |
| Increase over process start | 132.0 MiB |
| Projected generation time at 99,999 | 0.363 h (21.8 min) |
| Projected uncompressed statistic arrays at 99,999 | 8.383 GiB |

The disk projection includes copy counts, arm breadth, and collection
denominators for both assignments. It excludes compressed placement TSVs,
small manifests/results, filesystem metadata, and any additional annotation
collections beyond a combined total of 15,000 terms. Counts use `uint16` for
this 61,312-locus universe, breadth uses `uint8`, and the engine automatically
widens either type if a supplied universe cannot fit.

## Scaling and memory assessment

Generation scales approximately linearly in permutations and in the number of
selected locus-term edges. Extrapolating the measured 15,000-term workload gives
about 3.6 hours for 999,999 draws before filesystem contention and diagnostic
rescans. That extension is expensive but operationally realistic when limited
to decision-relevant runs.

The engine writes fixed-size `.npy` batches atomically and reads them by memory
map. Final inference processes 64 term columns at a time. At 99,999 draws, one
float64 64-term working matrix is about 48.8 MiB; simultaneous count,
composition, breadth, standardization, and summary temporaries remain bounded
well below 1 GiB rather than materializing the full 100,000 by 15,000 cube.
Leave-one diagnostics rescan compressed placement coordinates only for terms
that cross the prespecified q/maxT diagnostic threshold.

## Limitations and run planning

This benchmark measures one local process without competing moosefs load. It is
not a biological enrichment run because frozen GO/family/identity maps are a
separate input product. Actual term prevalence and annotation density can alter
runtime. Plan at least 12 GiB of free output space per 15,000-term primary run,
plus comparable space for each full background sensitivity retained on disk.
Run backgrounds sequentially when storage is constrained.

The benchmark establishes computational feasibility only. Reportable inference
still requires the input checks, minimum candidate counts, placement/calibration
gates, Monte Carlo decision intervals, and independent recount specified in
`STATISTICAL_SPEC.md`.
