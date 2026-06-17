# F3 Mouse Shape Contrast Report

Date: 2026-06-17  
Task: `manuscript-revision-f3`  
Script: `paper_prep/manuscript_revision/F3_mouse_shape/F3_mouse_shape_analysis.py`

## Question

The manuscript draft claimed that the mouse meiosis sequence-to-contact
correlation peaks at zygotene. The requested F-3 check is narrower than asking
whether each stage has rho different from zero: it tests whether zygotene rho is
higher than the flanking stages.

## Inputs

The analysis used only cached repo TSVs. No upstream matrix rebuild was run.

Primary 1 Mb, 50 kb resolution inputs:

| Analysis family | Cached inputs | Rows |
|---|---|---:|
| Sequence-level exact PHR | `data/mouse_meiosis_sweep/seqlevel/1Mb/mouse_<stage>_phr_50000bp_seqlevel.tsv` | 999-1,135 rows/stage |
| Arm-pair collapsed / community-analysis | `data/mouse_meiosis_sweep/zuo/1Mb/50000bp/zuo2021_<stage>_phr_pair_correlation.tsv` | 344 rows/stage |

The full file inventory is in `F3_input_inventory.tsv`.

## Method

For each stage I recomputed Spearman rho between sequence similarity and Hi-C
contact:

- Sequence-level exact PHR: `jaccard` versus `hic_contact_norm`.
- Arm-pair collapsed: `mean_jaccard` versus `hic_contact`.

For the shape test, I computed paired contrasts:

`delta = rho(zygotene on common rows) - rho(comparator stage on the same common rows)`

Contrasts were bootstrapped with the same paired row universe for each
comparator:

- Sequence-level exact PHR: clustered by unordered arm pair, because many PHR
  rows share the same arms.
- Arm-pair collapsed: arm-block bootstrap, because the 344 arm-pair entries are
  not independent of shared chromosome arms.

The script uses 2,000 bootstrap replicates with deterministic seed `2463`.
Two-sided bootstrap p-values are computed from the fraction of bootstrap
contrasts crossing zero. The p-values are decision-support values, not exact
permutation p-values.

## Slurm Record

No Slurm job was submitted. There is no Slurm job ID.

Reason: all cached inputs are small (`~16-198 KB` per file; 8 files total), and
the 2,000-replicate clustered bootstrap completed locally in about 30 seconds.
This did not require matrix scans, large RDS loads, or heavy permutation loops.
A 20,000-replicate pure-Python pilot was interrupted as unnecessarily slow for
this decision record; the committed script uses the lighter reproducible
bootstrap count.

## Stage Series

Output table: `F3_stage_series.tsv`.

| Analysis | Stage | rho | Rows | Resampling unit |
|---|---:|---:|---:|---|
| Sequence-level exact PHR | leptotene | 0.372 | 1,088 | arm pair |
| Sequence-level exact PHR | zygotene | 0.425 | 1,135 | arm pair |
| Sequence-level exact PHR | pachytene | 0.428 | 1,135 | arm pair |
| Sequence-level exact PHR | diplotene | 0.416 | 999 | arm pair |
| Arm-pair collapsed | leptotene | 0.680 | 344 | arm |
| Arm-pair collapsed | zygotene | 0.715 | 344 | arm |
| Arm-pair collapsed | pachytene | 0.677 | 344 | arm |
| Arm-pair collapsed | diplotene | 0.574 | 344 | arm |

These reproduce the apparent tension in existing manuscript-era notes:

- The exact-coordinate sequence-level PHR series is lower (`rho ~0.37-0.43`)
  and does not have a clean zygotene maximum; pachytene is slightly higher than
  zygotene in the all-row stage series.
- The arm-pair collapsed series is higher (`rho ~0.57-0.71`) and has the
  largest point estimate at zygotene.

## Direct Zygotene Contrasts

Output table: `F3_zygotene_contrasts.tsv`.

| Analysis | Contrast | rho zygotene | rho comparator | delta rho | 95% bootstrap CI | p(two-sided) |
|---|---|---:|---:|---:|---:|---:|
| Sequence-level exact PHR | zygotene - leptotene | 0.434 | 0.372 | +0.062 | 0.004 to 0.120 | 0.033 |
| Sequence-level exact PHR | zygotene - pachytene | 0.425 | 0.428 | -0.002 | -0.051 to 0.049 | 0.944 |
| Sequence-level exact PHR | zygotene - diplotene | 0.408 | 0.416 | -0.008 | -0.079 to 0.059 | 0.812 |
| Arm-pair collapsed | zygotene - leptotene | 0.715 | 0.680 | +0.035 | -0.148 to 0.253 | 0.810 |
| Arm-pair collapsed | zygotene - pachytene | 0.715 | 0.677 | +0.038 | -0.153 to 0.266 | 0.772 |
| Arm-pair collapsed | zygotene - diplotene | 0.715 | 0.574 | +0.141 | -0.143 to 0.452 | 0.440 |

Interpretation:

- The only contrast whose bootstrap interval excludes zero is the
  sequence-level exact-PHR zygotene-minus-leptotene contrast.
- Zygotene is not higher than pachytene in the sequence-level analysis.
- Zygotene is not higher than diplotene in the sequence-level analysis on the
  common-row paired set.
- The arm-pair collapsed point estimates are highest at zygotene, but the
  arm-block bootstrap intervals are wide and all include zero.

Therefore, the direct shape test does not support an unqualified "zygotene peak"
claim. At most, it supports a descriptive statement that the arm-pair collapsed
point estimate is numerically largest at zygotene, while clustered contrasts do
not establish a statistically resolved zygotene-specific peak.

## Per-PHR-Pair Versus Arm-Level Tension

The two analysis families answer related but different questions.

The sequence-level exact-PHR table treats individual PHR sequence pairs at their
exact coordinates. This is the most local test, but rows are highly
non-independent because many sequence pairs share the same arms and the same
stage-specific contact structure. It also includes multiple PHR pairs per arm
pair, so within-arm-pair heterogeneity and multi-mapping/contact sparsity can
lower the rank correlation.

The arm-pair collapsed table averages sequence similarity into one mean
Jaccard/contact observation per unordered arm pair. This suppresses noisy
within-arm-pair variation and makes the sequence-to-contact relationship look
much stronger (`rho ~0.57-0.71`). But the nominal 344 arm-pair rows are still
not 344 independent observations because there are only 27 arms. Once the
contrast is resampled at the arm level, the zygotene-minus-flanking intervals
are wide.

Reconciliation: both views support a broad positive mouse meiosis
sequence-to-contact relationship across all stages. They do not jointly support
a precise zygotene-specific curve-shape claim. The per-PHR view weakens the
shape claim because pachytene/diplotene are comparable to zygotene; the
arm-level view weakens it because the apparent zygotene maximum is not stable
under arm-block uncertainty.

## Manuscript Recommendation

Do not write that the mouse correlation "peaks at zygotene" as an inferential
result.

Exact Results replacement:

> In mouse meiotic Hi-C, subtelomeric sequence similarity was positively
> correlated with inter-chromosomal contact across leptotene, zygotene,
> pachytene and diplotene. At 1 Mb and 50 kb resolution, exact PHR-pair
> correlations were similar across stages (Spearman rho 0.37-0.43), while
> arm-pair-collapsed correlations were higher (rho 0.57-0.71) with the largest
> point estimate at zygotene. Direct clustered contrasts did not resolve a
> zygotene-specific peak, so we treat the mouse data as evidence for a broad
> prophase-I sequence-to-contact association rather than a stage-specific
> maximum.

Exact caption replacement for any mouse panel:

> Mouse meiotic Hi-C shows a positive subtelomeric sequence-to-contact
> association at all four prophase-I stages. The arm-pair-collapsed point
> estimate is largest at zygotene, but clustered zygotene-versus-flanking-stage
> contrasts do not establish a statistically resolved zygotene-specific peak.

If space is constrained, use this shorter caption sentence:

> Mouse meiotic Hi-C supports a broad prophase-I sequence-to-contact
> association; zygotene is the largest arm-collapsed point estimate but not a
> resolved stage-specific peak by clustered contrast.

## Validation Checklist

- Artifact report exists: this file.
- Zygotene-vs-flanking contrast tested directly: `F3_zygotene_contrasts.tsv`.
- Per-pair vs arm-level tension explained: see "Per-PHR-Pair Versus Arm-Level
  Tension".
- Slurm job ID or not-run reason recorded: see "Slurm Record".
- Exact manuscript/caption recommendation provided: see "Manuscript
  Recommendation".
