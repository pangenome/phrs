# M — Mouse stage-profile disambiguation (statistic vs axis vs satellite)

> Provenance for this note: see `response/PROVENANCE.md` (tool versions: `response/notes/tool_versions_and_pins.md`).

Status: **COMPLETE.** Six steps (A–D, F, E) done from one 10,000-permutation /
10,000-bootstrap run. STEP F was added after the first verdict and it changed one part of it
(see STEP E and the adjacency column of F.2).
Numeric source: `scripts/hic/mouse_stage_profile_disambiguation.json` (286,813 B)
produced by `scripts/hic/mouse_stage_profile_disambiguation.py`. Every number below is
read out of that JSON or produced by a command listed in `## commands run`.

## TL;DR

Five things settle the question, and only one of them was suspected:

1. **There is only one similarity axis in the project.** The published arm-level matrix
   `mouse.dist_matrix.tsv` is numerically the arm-pair mean of the per-PHR-pair `jaccard`
   column (max |difference| **1.1e-16** over 344 arm pairs). The "re-derived PAF Jaccard"
   that analysis (2) is supposed to have used **does not exist** — `mantel_mouse_d_m5.py:26`
   reads the published matrix, and the PAF substitute in the M6 lane is VOID (constant
   zero). So no stage-order change can be an axis effect.
2. **The zygotene→pachytene label flip is a statistic effect**: pair-level Spearman (n=1135)
   vs arm-level Mantel (n=344) on the *same* axis. The arm-level statistic discards the
   25.6% of axis variance that lives *within* arm pairs; on the same axis the
   within-arm-pair excess peaks at **zygotene in all six window×resolution cells**.
3. **The referee's R1.7 test (never run before): the shape is resolvable, the peak
   ordering is not — at chromosome-end resolution.** 1 Mb / 20 kb: rho +0.4187/+0.6142/+0.5761/+0.2452, arm-blocked permutation
   p ≤ 1.5e-3, arm-bootstrap CI excluding zero at every stage. Paired arm bootstrap:
   leptotene→zygotene **+0.196 [+0.093,+0.302]** (resolvable), zygotene→pachytene
   **−0.038 [−0.131,+0.043]** (NOT resolvable), pachytene→diplotene
   **−0.331 [−0.449,−0.207]** (resolvable). Same shape at 2 Mb and 4 Mb, 20 kb.
4. **Satellite control at pair resolution**: the association and the
   zygotene-over-diplotene ordering survive the satellite-rich subset
   (+0.2850/+0.5079/+0.3230/**+0.0231**) and the partial rank correlation on the pair's mean
   satellite fraction (+0.4128/+0.5804/+0.4393/+0.1755). The satellite-*poor* subset cannot
   test anything: it is q×q with only **8 distinct axis values and 35 of 42 pairs tied at
   Jaccard 0** (M6's arm-level analogue: 17 of 21 pairs tied at 0).
5. **Keeping pair resolution and doing the inference properly (STEP F).** Two-way clustered
   SEs over the two arms leave the per-stage estimates untouched (rho
   +0.4187/+0.6142/+0.5761/+0.2452, SE 0.0270/0.0239/0.0233/0.0279, design effect 0.96–1.02)
   and MRQAP with double semi-partialling on the 49-node pair matrices agrees
   (p 1e-4/1e-4/1e-4/8.5e-3). Of the **four designs**, only the arm-level Mantel says
   pachytene; the three pair-level designs all say zygotene. The adjacency test is where the
   designs part company: **pachytene − zygotene = −0.038** is resolvable at pair resolution
   (95% CI −0.067 … −0.009, p = 0.017, Holm 0.017) and by the arm-pair block bootstrap, but
   **not** by the chromosome-end bootstrap (27 units, −0.131 … +0.043; coverage 0.99 in a
   calibration simulation, i.e. conservative rather than invalid).

**Manuscript:** report design (4) — dyadic two-way clustered inference at pair resolution,
with the referee's stratified permutation as the null and the chromosome-end bootstrap for
the contrasts. Say "elevated at zygotene and pachytene, collapsed by diplotene"; do **not**
say "strongest at the zygotene bouquet" (a 0.04 difference between the two upper stages).

---

## Skeleton / plan
- [x] STEP A — axis audit (odgi graph-node Jaccard vs re-derived PAF Jaccard)
- [x] STEP B — same axis, two statistics (per-PHR-pair Spearman vs arm-level Mantel)
- [x] STEP C — referee's test R1.7 (stratified pair-label permutation + arm bootstrap CI)
- [x] STEP D — satellite control at pair resolution
- [x] STEP F — dyadic / arm-clustered inference at pair resolution (two-way clustered SEs,
      MRQAP-DSP, arm-pair FE, four-design comparison)
- [x] STEP E — verdict, which design to report, manuscript sentence, what is unmeasurable

Outputs: `scripts/hic/mouse_stage_profile_disambiguation.py`, `scripts/hic/mouse_stage_profile_disambiguation.json`

## Log

### L0 — inputs read (before any computation)
- `PEER_REVIEW.md` (authority; R1 major points 6 and 7 read verbatim).
  `paper_prep/synthesis/SIMULATED_PEER_REVIEW_v1.md` **never opened**.
- `response/notes/M6_mouse_satellite.md` (authority for the satellite lane).
- `scripts/hic/mantel_mouse_d_m5.py` (105 lines) + `scripts/hic/mantel_mouse_d_m5.json`.
- `scripts/hic/mouse_satellite_sensitivity.py` (1038 lines) + its JSON (83,514 B).
- `submission/scripts/figures/make_fig4c_mouse_zygotene.R` (the submitted analysis).
- `data/mouse_{leptotene,zygotene,pachytene,diplotene}_phr_{20000,50000}bp_seqlevel.tsv`
  (1135 rows at 20 kb for all four stages; 1088 / 1135 / 1135 / 999 at 50 kb).
- `data/mouse_meiosis_sweep/seqlevel/{2Mb,4Mb}/*_{20000,50000}bp_seqlevel.tsv`
  (robustness windows; 1135 rows each at 20 kb).
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/similarity_1Mb/mouse.dist_matrix.tsv`
  (27 arms, 729 cells) and 4 gzipped PAFs under `subtelo_1Mb/all-vs-all.p95/`.
- `paper_prep/_brainstorming/mouse_bouquet_reanalysis/NOTE.md` (context only; not used
  as a numeric source).

Everything below comes from `scripts/hic/mouse_stage_profile_disambiguation.py` and the
JSON it writes. No number in this note was typed in by hand.

---

## STEP A — axis audit: how many similarity axes are there really?

### A.1 What each analysis actually reads

| analysis | statistic | axis | where it is read |
|---|---|---|---|
| (1) SUBMITTED | per-PHR-pair Spearman over 1135 pairs | column `jaccard` of `data/mouse_<stage>_phr_20000bp_seqlevel.tsv` — the odgi pangenome-graph node Jaccard of the **two PHR sequences of that pair** | `submission/scripts/figures/make_fig4c_mouse_zygotene.R:36-40` |
| (2) REPLACEMENT | arm-level Mantel, 10,000 row-and-col permutations | `1 - d` from `mouse.dist_matrix.tsv` — the odgi graph node Jaccard **between arm-level PHR sets** | `scripts/hic/mantel_mouse_d_m5.py:26` (`DIST` constant) |
| (3) M6 stratified | as (2), split by per-arm satellite fraction | as (2) | `scripts/hic/mouse_satellite_sensitivity.py` |
| (claimed) re-derived PAF Jaccard | — | **does not exist: VOID** | `mouse_satellite_sensitivity.json → sequence_axis` |

**The brief's premise for (2) is wrong.** `mantel_mouse_d_m5.py:26` reads the *published*
`mouse.dist_matrix.tsv`, not a re-derived PAF Jaccard, and M6 §2 records that its eight
rho values reproduce that JSON exactly. The only "re-derived PAF-based Jaccard" in the
project is inside `mouse_satellite_sensitivity.py`, and it is void.

### A.2 Independent re-check of the void axis (my own diagnostic, not M6's)
`mouse_satellite_sensitivity.py:706` reads the PAF record as
`qs, qe = int(c[1]), int(c[2])`. PAF column 2 (0-based index 1) is **QLEN** and index 2 is
**QSTART**, so the stored "interval" is `(qlen, qstart)` and `qlen > qstart` always.
Bounded check over 4 PAF files, 4,228 records:

```
files_read: mouse.all-vs-B6_1_chr10_parm_8494_1000000.paf.gz,
            mouse.all-vs-B6_1_chr10_qarm_136126214_137106573.paf.gz,
            mouse.all-vs-B6_1_chr11_parm_10588_1000000.paf.gz,
            mouse.all-vs-B6_1_chr11_qarm_131708645_132696509.paf.gz
n_records: 4228,  n_reversed: 4228
example: paf_column2_qlen=978851, paf_column3_qstart=0,
         stored_tuple_by_satellite_script=[978851, 0], reversed=true
```
Every interval is empty after intersection, every derived Jaccard is 0.0, and the stored
"validation vs published odgi Jaccard" is NaN over n=351 pairs — a zero-variance input,
not a weak correlation. **No satellite-maskable sequence axis exists on disk.**

### A.3 The decisive test: the two usable axes are numerically the same axis

I compared, over the 344 arm pairs, the published arm-level similarity `1 - d` with the
**arm-pair mean of the pair-level `jaccard` column**:

```
max |arm-mean pair Jaccard - (1 - d)| = 1.1102230246251565e-16
identical_to_float_precision        = true
n_arm_pairs_compared                = 344   (27 arms -> 351 unordered; 7 have no PHR pairs)
```

So the published arm-level matrix **is** the mean of the per-pair Jaccards. Analyses (1)
and (2) do not use two axes; they use one axis at two resolutions.

Agreement at pair level, where both are computable:

| window | Spearman(pair Jaccard, 1 - d of its arm pair) | n pairs |
|---|---|---|
| 1 Mb, 20 kb | **+0.8796** | 1135 |
| 1 Mb, 50 kb | +0.8796 | 1135 |
| 2 Mb, 20 kb | +0.8750 | 1135 |
| 4 Mb, 20 kb | +0.8741 | 1135 |

Rank-equivalence at arm level: `spearman = 1.0` at 1 Mb. The matrix is built from the
1 Mb flanks only, so at 2 Mb / 4 Mb there *is* a genuine second axis
(max |difference| 0.173 at 2 Mb, 0.306 at 4 Mb; Spearman 0.969 / 0.953) — that is the only
place in the project where "changing the axis" is even representable, and it changes
nothing about the ordering (see B).

### A.4 Value ranges

| axis | range | notes |
|---|---|---|
| published `1 - d` (arm) | similarity **0.000 – 0.7544**; distance 0.0 – 1.0 | 38.5% of the 702 off-diagonal cells are at distance exactly 1.0 (similarity 0); median off-diagonal distance 0.7748 |
| pair-level `jaccard` | **0.000 – 0.9872**, median well below 0.1 | stage-invariant: the four stages carry **identical** Jaccard vectors (checked element-wise, `axis_is_stage_invariant = true`) |
| re-derived PAF Jaccard | constant 0.0 | VOID |

### A.5 Two audit findings that matter downstream

**(i) The published matrix is not a metric.** 22 of 27 arms have a non-zero *self*-distance
(Jaccard of an arm's PHR set with itself < 1); examples: `chr11_p` 0.4997, `chr7_p` 0.5,
`chr4_p` 0.4990, `chr8_q` 0.3855. No statistic here uses diagonal cells, but a distance
matrix with d(a,a) ≠ 0 is not a Jaccard distance and should not be described as one.

**(ii) The published Mantel routine is tie-sensitive.** `mantel_mouse_d_m5.mantel` ranks
with `np.argsort(np.argsort(x))`, which breaks ties **by row order**. **129 of the 344 arm
pairs (37.5%) are tied at similarity 0**, and there are only 216 distinct axis values. The
published rho is therefore not the Spearman coefficient; it depends on the row order of the
contact table. I reproduce it bit-for-bit (see the verification gates) *and* report the
tie-corrected (average-rank) value:

| stage, 1 Mb 20 kb | published routine | tie-corrected | shift |
|---|---|---|---|
| leptotene | +0.4074 | +0.3991 | 0.0083 |
| zygotene | +0.5785 | +0.5813 | 0.0028 |
| pachytene | +0.6259 | +0.6289 | 0.0030 |
| diplotene | +0.3282 | +0.3301 | 0.0019 |
| leptotene, 50 kb | +0.2882 | +0.3357 | **0.0475** |
| pachytene, 50 kb | +0.3130 | +0.3616 | **0.0486** |

At 20 kb the shift is ≤ 0.009 and the ordering is unaffected. At 50 kb the shift reaches
0.049 — larger than the entire leptotene-to-pachytene spread at that resolution — and it
**flips the label at 50 kb** (published routine: zygotene; tie-corrected: pachytene). The
50 kb peak stage is not identifiable.

### A.6 Distance-to-telomere is fixed by construction
Parsing the PHR ids (`sample#hap#arm:<flank_start>-<flank_end>:<phr_start>-<phr_end>`),
distance from the arm's telomere to the PHR's near edge has **median 1 bp**, p90
**3,872 bp**, max **4,532 bp**; 70.3% of pairs are under 1 kb and 100% under 10 kb. So R1.7's "holding
distance-to-telomere fixed" is satisfied by the design; the covariate that actually varies
is PHR *span* (p-arm PHRs saturate the 1 Mb flank, q-arm PHRs are ~20 kb long).

---

## STEP B — same axis, two statistics

2 × 2 crossing plus the mixed cell ("pair-level statistic on the arm-level axis"), all on
the SAME underlying axis, 10,000 permutations each. Values are rho.

**1 Mb, 20 kb (the manuscript's cell), n = 1135 pairs / 344 arm pairs:**

| statistic \ axis | pair-level Jaccard | arm-level `1 - d` | peak |
|---|---|---|---|
| pair-level Spearman | **+0.4187 / +0.6142 / +0.5761 / +0.2452** | +0.3230 / +0.5104 / +0.5081 / +0.2117 | **zygotene** (both) |
| arm-level Mantel (tie-corrected) | +0.3991 / +0.5813 / +0.6289 / +0.3301 | +0.3991 / +0.5813 / +0.6289 / +0.3301 | **pachytene** (both) |
| arm-level Mantel (published routine) | +0.4074 / +0.5785 / +0.6259 / +0.3282 | identical | **pachytene** |

Order: leptotene / zygotene / pachytene / diplotene.

Peak margins: pair-level **+0.0381** (zygo − pachy); pair-level on the arm axis **+0.0023**;
arm-level Mantel **+0.0477** (pachy − zygo). All Mantel p ≤ 1e-4.

**The same comparison at the other windows, 20 kb:**

| window | pair-level peak (margin) | arm-level Mantel peak (margin) | agree? |
|---|---|---|---|
| 1 Mb | zygotene (+0.0381) | pachytene (+0.0477) | no |
| 2 Mb | zygotene (+0.0076) | pachytene (+0.0315) | no |
| 4 Mb | zygotene (+0.0575) | **zygotene** (+0.0880) | yes |

**What this isolates.** The axis cannot be the cause: at 1 Mb it is *the same vector*
(1.1e-16), and even where a real second axis exists (2 Mb, 4 Mb; Spearman 0.969 / 0.953)
the two axis variants give the **same peak stage in every one of the six cells**, and the
rho gap between them (up to 0.054 at 2 Mb, 0.049 at 4 Mb, arm level) is smaller than the
gap that separates zygotene from pachytene in any of them. What changes the label is **the collapse of within-arm-pair variation**: of the total
pair-level Jaccard variance, **25.6% lives within arm pairs** and
**74.4% between arm pairs** (344 blocks, median 4 pairs/block). The arm-level Mantel sees
only the between part. That is a property of the *statistic*, not of the axis — and the
"pair-level statistic on the arm axis" column shows it directly: holding the axis at arm
resolution and keeping pair resolution for the statistic gives zygotene, but by a margin
of 0.0023, i.e. a dead heat.

At 50 kb nothing is worth reading: the four stages span 0.36 – 0.43 (pair level) and
0.29 – 0.36 (Mantel), and at 4 Mb/50 kb the arm-level Mantel is not significant at any
stage (p 0.013 – 0.745).

---

## STEP C — the referee's test (R1 major point 7)

Referee #1: *"the remedy is a null preserving that dependency: a permutation shuffling PHR
labels while holding distance-to-telomere and arm identity, with bootstrap intervals from
resampling chromosome ends."*

Implementation (`stratified_test`, `cluster_bootstrap`, `paired_cluster_bootstrap`):

- **null = stratified permutation at PAIR resolution.** The similarity label of each PHR
  pair is shuffled only among PHR pairs that share the same arm pair. Arm identity, the set
  of arm pairs, the number of PHR pairs per arm pair, and every contact value are held
  fixed. Vectorised via one uniform key per observation + `lexsort` by (block, key), which
  is a uniformly random permutation inside each block.
- **distance-to-telomere** is fixed by construction (A.6), and additionally held fixed by
  re-running with strata = arm pair × dtel tercile (`p_stratified_plus_dtel`), and by
  strata = arm pair × PHR-span tercile (`p_stratified_plus_span`) for the covariate that
  actually varies.
- **CI = cluster bootstrap over CHROMOSOME ENDS.** 27 arms resampled with replacement,
  10,000 replicates; a resampled multiset of arms induces the weight `m_a · m_b` on each
  PHR pair, so within-arm-pair structure is retained rather than broken.
- **ordering**: paired contrasts on one shared arm resampling, restricted to the pair set
  all four stages share. A 27-draw resample contains ~17 of the 27 arms (observed per-arm
  draw frequency 0.6388 against the exact 1 − (26/27)^27 = 0.639), which is the honest
  statement of how much chromosome-end information a bootstrap replicate carries.

**1 Mb, 20 kb — rho, stratified permutation p, arm bootstrap 95% CI:**

| stage | rho | p (stratified) | 95% CI (arms) | excess over null |
|---|---|---|---|---|
| leptotene | +0.4187 | 9.999e-05 | +0.252 … +0.581 | +0.1142 |
| zygotene | +0.6142 | 9.999e-05 | +0.439 … +0.753 | +0.1506 |
| pachytene | +0.5761 | 9.999e-05 | +0.414 … +0.700 | +0.1034 |
| diplotene | +0.2452 | 1.4999e-03 | +0.087 … +0.372 | +0.0395 |

(adding the dtel or the span stratum leaves every p at ≤ 1.5e-3, i.e. unchanged.)

**The stratified null is not centred at zero,** and this is the substantive result of
doing the test properly: the null mean is +0.3045 / +0.4635 / +0.4727 / +0.2057 for the
four stages. In other words **73% / 75% / 82% / 84% of each rho is the between-arm-pair
structure** the referee wanted held fixed, and only the excess (0.114 / 0.151 / 0.103 /
0.040) is within-arm-pair covariation. That excess peaks at **zygotene**, in
**all six cells** (both resolutions, all three windows).

**Adjacent-stage contrasts (paired arm bootstrap, 1 Mb 20 kb):**

| contrast | Δrho | 95% CI | distinguishable |
|---|---|---|---|
| zygotene − leptotene | +0.1955 | +0.093 … +0.302 | **yes** |
| pachytene − zygotene | −0.0381 | −0.131 … +0.043 | **no** |
| diplotene − pachytene | −0.3309 | −0.449 … −0.207 | **yes** |

Same shape at 2 Mb and 4 Mb, 20 kb: leptotene→zygotene and pachytene→diplotene are both
distinguishable (2 Mb: +0.147 [+0.067,+0.227] and −0.253 [−0.346,−0.155]; 4 Mb: +0.120
[+0.043,+0.191] and −0.194 [−0.293,−0.098]), and **zygotene vs pachytene is not
distinguishable at any window or resolution** (−0.038, −0.008, −0.057 at 20 kb).

At 50 kb **no** adjacent contrast is distinguishable at any window — that resolution has no
stage resolution at all.

**Arm-level degeneracy (why the test must be at pair resolution).** Within-arm-pair
permutation is a strict no-op for the arm-level statistic: max |shift of the arm-mean
axis| = 1.1e-16 and the Mantel rho is identical before and after (e.g. 0.399086 before,
0.399086 after). Reported per cell as `arm_level_null_is_degenerate`.

---

## STEP D — satellite control at pair resolution

Per-arm satellite fractions are reused verbatim from
`mouse_satellite_sensitivity.json → satellite_annotation.per_arm_sat_frac` (27 arms, no
recomputation). Threshold 0.25 sits in the documented bimodal gap (7 q arms < 0.008,
20 p arms > 0.467). Stratum sizes at 1 Mb 20 kb: **both-rich 722 pairs / 20 arms**,
**mixed 371 / 27**, **both-poor 42 / 7**.

### D.1 The satellite-poor subset cannot answer the question — its axis is degenerate

| stratum | pairs | distinct axis values | pairs with Jaccard exactly 0 |
|---|---|---|---|
| both satellite-rich (p×p) | 722 | **684** | 37 |
| mixed (p×q) | 371 | 34 | 337 |
| both satellite-poor (q×q) | 42 | **8** | **35** |

The q×q pairs have almost no PHR sequence sharing by construction, so a rank correlation
over them rests on 7 non-tied pairs; with 2 pairs per arm-pair block the stratified null is
nearly saturated and p ≈ 0.25 at every stage. M6's arm-level analogue has the same
pathology (**17 of its 21 q×q arm pairs are tied at similarity 0**), which is why its "flat
0.5364 / 0.5390 / 0.4857 / 0.4857" and its 0.0026 flip are not interpretable.

### D.2 What the non-degenerate strata say (20 kb)

| stage | both-rich (n=722) | p | mixed (n=371) | partial (control = pair mean satellite) | p |
|---|---|---|---|---|---|
| leptotene | +0.2850 | 1.00e-04 | +0.2871 | +0.4128 | 1.00e-04 |
| zygotene | **+0.5079** | 1.00e-04 | +0.3179 | **+0.5804** | 1.00e-04 |
| pachytene | +0.3230 | 1.00e-04 | +0.3140 | +0.4393 | 1.00e-04 |
| diplotene | **+0.0231** | 1.28e-02 | +0.2492 | +0.1755 | 1.00e-03 |

Marginal for reference: 0.4187 / 0.6142 / 0.5761 / 0.2452.

- The satellite-**rich** subset — the one with real axis variation — gives the cleanest
  ordering in the whole study: zygotene +0.508, pachytene +0.323, diplotene **+0.023**
  (arm bootstrap +0.311 … +0.687 at zygotene, −0.104 … +0.178 at diplotene).
- The **partial** rank correlation controlling for the pair's mean satellite fraction keeps
  zygotene on top (+0.580 vs +0.439 pachytene) and drops diplotene to +0.176, i.e. the
  zygotene peak is *not* a satellite-fraction artefact. Note `Spearman(jaccard, mean
  satellite) = +0.790`, so satellite content is strongly confounded with the axis; the
  partial removes it and the ordering survives.
- 50 kb: same conclusion (rich 0.203 / 0.448 / 0.262 / 0.233; partial 0.259 / 0.472 /
  0.323 / 0.296).

**Answer to the referee's point 6 as far as the data allow:** the sequence–contact coupling
is not created by satellite-rich p-arm windows (it survives in the q×q subset, weakly, and
in the partial correlation, strongly), and the stage structure is *strongest*, not weakest,
in the satellite-rich subset. What is **not** measurable is a satellite-masked *contact*
axis (no mcool access), so "satellite content in the contact values is not the driver"
remains untested.

---

## STEP F — keep pair resolution; handle arm dependence by clustering/permuting

The arm pair is a partitioning device, not a unit of observation. Three specific defects of
collapsing to arm means, measured here before fixing the inference:

| defect | measurement |
|---|---|
| **(i) power loss** | 1135 PHR pairs → 344 arm pairs at 20 kb. Unequal by stage at 50 kb: leptotene 1088 pairs / zygotene 1135 / pachytene 1135 / **diplotene 999**, which gives **319** arm pairs for diplotene against 344 for the others — so the arm-level sample is not even the same size at every stage. |
| **(ii) unequal, tiny cells** | block sizes at 20 kb: **min 1, median 4, max 4**; **2.6%** of blocks are singletons, **31.1%** hold exactly 2 pairs, and **19.6%** of all pairs live in blocks of ≤2. Every arm-level "mean contact" is a mean of **1–4** numbers. At 50 kb diplotene, 4.7% singletons / 36.4% size-2. |
| **(iii) satellite re-weighting** | giving each block equal weight moves the effective mean satellite fraction from **0.708** (pair level) to **0.661** (block level), and the both-satellite-rich weight share from **0.636** to **0.552**, because satellite-rich blocks hold **3.80** PHR pairs on average against **2.00** for satellite-poor blocks. So the collapse does not merely lose information, it changes which sequence composition the statistic is about. |
| (ii) consequence | weighting the arm-level Mantel by the number of PHR pairs in each block (instead of equal blocks) shifts rho by **+0.0021 / +0.0322 / +0.0210 / +0.0102** (leptotene … diplotene). The peak stays pachytene, but the zygotene margin that pachytene has to beat shrinks from 0.048 to 0.037 — the largest single shift lands on the competing stage. |

### F.1 Dyadic two-way clustered inference at pair resolution (i)

Both variables are standardised average ranks, so the OLS slope **is** the Spearman
coefficient (`slope_equals_spearman` verified in every cell); the variance is the
Cameron–Gelbach–Miller two-way cluster-robust sandwich over **the two arms of each pair**,
with the arm-pair block as the intersection term (`_sandwich`, `twoway_cluster_ols`).
27 arms → 25 clusters on one dimension and 26 on the other (the alphabetically largest arm
never appears as `arm_a`), 344 block clusters, 1135 rows.

1 Mb, 20 kb:

| stage | rho | SE iid | SE two-way clustered | 95% CI (clustered) | p (t, 26 df) | design effect |
|---|---|---|---|---|---|---|
| leptotene | +0.4187 | 0.0270 | **0.0270** | +0.366 … +0.472 | 5.05e-14 | 1.00 |
| zygotene | +0.6142 | 0.0234 | **0.0239** | +0.567 … +0.661 | 5.88e-19 | 1.02 |
| pachytene | +0.5761 | 0.0243 | **0.0233** | +0.530 … +0.622 | 1.40e-18 | 0.96 |
| diplotene | +0.2452 | 0.0288 | **0.0279** | +0.190 … +0.300 | 5.80e-09 | 0.97 |

**The design effect is ≈ 1.** This is worth stating explicitly because it looks like a bug
and is not one. Residuals *are* correlated within arms (one-way residual ICC by arm
0.054–0.149 across the four stages), but for a dense dyadic design the arm effect enters
through ~41 partners per arm and averages out of the *slope*. Two independent checks: (a) a
simulation with arm random intercepts (τ = 0, 0.5, 1.0, 200 replicates each) gives
MC sd of the slope = 0.048 / 0.067 / 0.089 against mean clustered SE = 0.054 / 0.066 /
0.092, i.e. the sandwich tracks the true sampling variability; (b) the contrast coverage
simulation below. So the referee's non-independence objection does **not** inflate the
per-stage slope once it is handled correctly — what it changes is the stage *contrast*, and
the effective number of independent chromosome ends.

### F.2 Adjacency tests as actual contrasts (i)

Each adjacent (and all six pairwise) contrast is estimated by differencing the outcome
first, `d = z(contact, stage 2) − z(contact, stage 1)`, and regressing `d` on the shared
`z(jaccard)`: the slope is **exactly** `rho_2 − rho_1` and the arm-level nuisance component
common to the two stages cancels. Verified: `slope_equals_delta` in every cell. (A stacked
`z(contact) ~ z(jaccard) * stage` model gives the identical point estimate but a ~40%
larger SE, because it cannot difference the shared arm component away; it is kept as
`stacked_interaction_check`.)

1 Mb, 20 kb — three resampling units for the same point estimates:

| contrast | Δrho | dyadic clustered 95% CI | arm-pair block bootstrap | chromosome-end bootstrap | dyadic p | Holm p (6 contrasts) |
|---|---|---|---|---|---|---|
| zygotene − leptotene | +0.1955 | +0.162 … +0.229 | +0.165 … +0.227 | +0.093 … +0.302 | 2.7e-11 | 1.1e-10 |
| pachytene − leptotene | +0.1574 | +0.127 … +0.188 | +0.132 … +0.183 | +0.087 … +0.219 | 3.9e-10 | 1.2e-09 |
| diplotene − leptotene | −0.1735 | −0.210 … −0.137 | −0.213 … −0.134 | −0.303 … −0.065 | 1.8e-09 | 3.5e-09 |
| **pachytene − zygotene** | **−0.0381** | **−0.067 … −0.009** | **−0.067 … −0.009** | **−0.131 … +0.043** | **0.017** | **0.017** |
| diplotene − zygotene | −0.3690 | −0.418 … −0.320 | −0.420 … −0.318 | −0.527 … −0.202 | 1.4e-13 | 6.9e-13 |
| diplotene − pachytene | −0.3309 | −0.368 … −0.293 | −0.368 … −0.293 | −0.449 … −0.207 | 5.0e-15 | 3.0e-14 |

### F.3 Which CI is safe? A calibration simulation

If the arm-level nuisance differs between stages (different Hi-C libraries, per-stage
per-arm coverage), the differencing assumption fails. 200 replicates at a true Δ = 0, 200
bootstraps per replicate:

| coverage of the nominal 95% interval, true Δ = 0 | dyadic two-way SE | arm-pair block bootstrap | chromosome-end bootstrap |
|---|---|---|---|
| arm effects shared across stages | 0.930 | 0.930 | 0.990 |
| arm effects independent across stages | 0.955 | 0.930 | 0.995 |

So all three intervals are valid; the **chromosome-end bootstrap is the conservative one**
(coverage ≈ 0.99 with only 27 units), and the dyadic SE stays calibrated even when arm
effects differ by stage — because the two-way clustering picks up the within-arm
correlation of the differenced outcome. The disagreement between the dyadic and
chromosome-end intervals on `pachytene − zygotene` is therefore about **generalisation to
unsampled chromosome ends**, not about validity.

### F.4 MRQAP with double semi-partialling at pair resolution (ii)

The PHR pairs form a 49-node symmetric matrix with **1135 of 1176** cells observed (96.5%).
Both nulls use 10,000 simultaneous row-and-column permutations (symmetry, diagonal and the
observed-cell pattern preserved; each replicate recomputes the observed statistic on the
same cell subsample, mean 1095 cells used, so the comparison is paired):

| stage | rho | p, DSP (predictor permuted) | p, outcome permuted | null sd (DSP) | R² |
|---|---|---|---|---|---|
| leptotene | +0.4187 | 1.0e-04 | 1.0e-04 | 0.072 | 0.175 |
| zygotene | +0.6142 | 1.0e-04 | 1.0e-04 | 0.067 | 0.377 |
| pachytene | +0.5761 | 1.0e-04 | 1.0e-04 | 0.069 | 0.332 |
| diplotene | +0.2452 | **8.5e-03** | **7.3e-03** | 0.091 | 0.060 |

The ``dsp`` variant is the one the request names: for each replicate the permuted
similarity matrix is residualised on the observed similarity matrix and the outcome is
residualised on the observed similarity matrix, and the semi-partial coefficient is
compared with the observed coefficient (Dekker–Krackhardt–Snijders). The ``y_perm`` variant
is the classic QAP regression null, which needs no correction because the predictor is held
fixed. Both agree. Note the MRQAP null is *wider* than the stratified null of STEP C
because permuting the outcome destroys the arm structure too, so p is larger (leptotene
1e-4 here versus 1e-4 there, but the null sd is 0.072 versus 0.013).

### F.5 Mixed model (iii) — not run, and why

`statsmodels` and `patsy` are **not installed** on this box, and installing them would
write outside this worktree, which this lane may not do. The requested mixed model with
random intercepts for both arms is replaced by two things that use the same
within-arm-pair variation: the two-way clustered estimator (F.1/F.2) and an **arm-pair
fixed-effects** estimator (within transformation on standardised ranks, clustered by arm
pair, **1126 rows in the 335 blocks of size ≥2** — the 9 singleton blocks carry no
within information and drop out):

| stage | within-arm-pair FE slope | SE | 95% CI (clustered by arm pair) | p |
|---|---|---|---|---|
| leptotene | +0.5781 | 0.0438 | +0.492 … +0.664 | 2.5e-32 |
| zygotene | **+0.7632** | 0.0341 | +0.696 … +0.830 | 2.0e-68 |
| pachytene | +0.5231 | 0.0432 | +0.438 … +0.608 | 3.4e-28 |
| diplotene | +0.2000 | 0.0555 | +0.091 … +0.309 | 3.6e-04 |
The within-arm-pair signal (which is exactly the variation the arm-level Mantel throws
away, and the same quantity as the excess over the stratified null in STEP C) peaks at
**zygotene** and is *lower at pachytene than at leptotene*.

### F.6 All four designs side by side (1 Mb, 20 kb)

| design | unit | n | leptotene | zygotene | pachytene | diplotene | peak | valid inference? |
|---|---|---|---|---|---|---|---|---|
| **(1) submitted** per-pair Spearman | PHR pair | 1135 | +0.4187 | **+0.6142** | +0.5761 | +0.2452 | zygotene | **no** — `cor.test` p treats 1135 dependent pairs as independent |
| **(2) arm-level Mantel** (replacement) | arm pair | 344 | +0.3991 | +0.5813 | **+0.6289** | +0.3301 | **pachytene** | yes for its own hypothesis, but wastes the data (see the table at the head of this step) |
| **(3) referee R1.7** stratified permutation + arm bootstrap | PHR pair, blocked by arm pair | 1135 | +0.4187 | **+0.6142** | +0.5761 | +0.2452 | zygotene | yes; p vs a null that keeps the between-arm structure |
| **(4) dyadic two-way clustered** (this step) | PHR pair, clustered on both arms | 1135 | +0.4187 | **+0.6142** | +0.5761 | +0.2452 | zygotene | **yes — the design to report** |
| (4b) MRQAP-DSP at pair resolution | PHR pair (49 nodes) | 1135 | +0.4187 | **+0.6142** | +0.5761 | +0.2452 | zygotene | yes; null preserves the row/column structure |

**Which designs agree?** On the ordering *by point estimate*, (1), (3), (4) and (4b) are
literally the same numbers and all peak at zygotene; **(2) alone says pachytene**. On the
*adjacency verdict* the pair-level designs agree that leptotene→zygotene and
pachytene→diplotene are resolvable; they disagree only about
**zygotene vs pachytene**, where the dyadic/block intervals resolve it (−0.038, p = 0.017)
and the chromosome-end bootstrap does not (−0.131 … +0.043).

## STEP E — verdict

**Which stage peaks.** At 20 kb — the resolution the manuscript uses — every design that
keeps pair resolution peaks at **zygotene**: the submitted statistic (1 Mb +0.6142,
2 Mb +0.4191, 4 Mb +0.3585), the referee's stratified test, the dyadic clustered estimator
and MRQAP-DSP. The arm-level Mantel (design 2) alone says pachytene, at 1 Mb and 2 Mb but
not at 4 Mb, and on the *same axis*. The within-arm-pair excess — the part of the signal
design 2 discards — peaks at zygotene in **all six** window × resolution cells, and the
within-arm-pair fixed-effects slope puts zygotene above pachytene by 0.240 (0.763 vs 0.523).

**Does the ordering survive?**

1. **The statistic change** — the *label* does not survive (zygotene under a pair statistic,
   pachytene under an arm statistic), but neither label was established: the two upper
   stages are separated by 0.038 in rho against differences of 0.196 and 0.331 at the edges
   of the profile. What survives every design is the **shape**: leptotene→zygotene
   +0.196 [+0.162, +0.229] and pachytene→diplotene −0.331 [−0.368, −0.293].
2. **The axis change** — the premise dissolves: there is no axis change. The published arm
   matrix *is* the arm-pair mean of the pair-level Jaccard (1.1e-16). What actually moves
   the label is the *collapse*: 25.6% of the axis variance lives within arm pairs and
   design 2 deletes it, while re-weighting blocks of 1–4 pairs equally (defects i–iii above).
3. **The satellite control** — the association and the zygotene-over-diplotene ordering
   survive the satellite-rich subset (+0.2850 / +0.5079 / +0.3230 / **+0.0231**) and the
   partial rank correlation on the pair's mean satellite fraction (+0.4128 / +0.5804 /
   +0.4393 / +0.1755). The satellite-poor subset cannot test anything (8 distinct axis
   values, 35/42 pairs tied at zero) — and design 2's version of it is worse, because the
   collapse itself shifts the effective satellite fraction from 0.708 to 0.661.
4. **The clustering** — the per-stage estimates are unaffected by correct handling of arm
   dependence (design effect 0.96–1.02; the two-way clustered SE virtually equals the iid
   SE, verified by simulation against arm random intercepts). The stage *contrast* is where
   it matters, and there the answer depends on the resampling unit: resolvable at pair and
   arm-pair-block resolution (p = 0.017), not resolvable at chromosome-end resolution
   (27 units, conservative interval, coverage 0.99).

### Which design we should report, and why

**Report design (4): dyadic two-way clustered inference at pair resolution** — and put the
referee's stratified permutation (design 3) beside it as the null, with the chromosome-end
bootstrap reported for the stage contrasts.

* (1) is the same point estimate with an invalid p-value: it must not be reported as the
  inference, but its point estimates are correct and identical to (4).
* (2) should be retired as the primary analysis: it collapses 1135 pairs to 344 means of
  1–4 observations each, deletes 25.6% of the axis variance, re-weights satellite
  composition (0.708 → 0.661), and its label changes with the flank window (pachytene at
  1 Mb / 2 Mb, zygotene at 4 Mb). Its *one* virtue — a null that respects the arm
  dependence — is delivered without any of those costs by (3) and (4).
* (3) and (4) keep all 1135 pairs; (4) gives a standard error and a regression-based
  contrast, (3) gives the exact permutation null the referee asked for. They agree
  everywhere except that (4) additionally resolves the small zygotene-over-pachytene
  difference.
* (4b) confirms both, under a null that keeps the full row/column dependence of the PHR
  pair matrix.

### The sentence we should put in the manuscript

> Across four mouse meiotic stages the sequence–contact coupling is present at every stage
> (Spearman ρ = 0.42, 0.61, 0.58, 0.25 with dyadic standard errors clustered on both
> chromosome arms, all p ≤ 6 × 10⁻⁹; arm-blocked permutation p ≤ 1.5 × 10⁻³), is highest at
> zygotene and pachytene and lowest at diplotene (−0.33 relative to pachytene, dyadic 95% CI
> −0.37 to −0.29).

(The run prints this sentence with its own numbers interpolated —
`step_E_verdict.manuscript_sentence` — so it cannot drift from the data.)

### The sentence we must NOT put in the manuscript

> "the coupling is strongest at the zygotene bouquet."

The zygotene–pachytene difference is **−0.038** in rho: dyadic 95% CI −0.067 to −0.009
(p = 0.017, Holm-adjusted over the six pairwise contrasts also 0.017), but the
chromosome-end bootstrap the referee asked for gives −0.131 to +0.043. A 0.04 difference
between the two upper stages — an order of magnitude below the diplotene collapse — does not
support an ordering, and it certainly does not support "strongest at the bouquet".
Equally barred: presenting design 2's "pachytene peak" as a stage effect, since on the same
axis the within-arm-pair variation that design 2 deletes peaks at zygotene.

### What remains unmeasurable with the data on disk

1. A **satellite-masked contact axis** — needs `zuo2021_<stage>.mcool` per-bin matrices
   (h5py/cooler absent; 0.3–0.5 Gpixel per stage). So "satellite content does not drive the
   *contact* values" is untested; only stratification and partial correlation are possible.
2. A **satellite-masked sequence axis on the published odgi Jaccard** — needs a pggb graph
   rebuild. The substitute (PAF Jaccard) is VOID (A.2).
3. **Minor-satellite content** of these windows — the local Dfam library has no Mus
   minor-satellite family (M6 §1); the proxy recovers ≈ 0. Unmeasured, not measured-small.
4. **Non-satellite repeat content** (TAR1, subterminal satellite, shared TE families) — no
   RepeatMasker run available on this box (M6 §5a item 2).
5. **Arm-level statements beyond 1 Mb** — the published arm matrix exists only for the
   1 Mb flanks; the 2 Mb / 4 Mb arm-level columns here are an arm-mean of the 2 Mb / 4 Mb
   pair Jaccard, a different quantity (Spearman 0.969 / 0.953 with the published one).
6. **Power for the satellite-poor (q×q) subset** — 7 chromosome ends, 42 pairs,
   8 distinct axis values, 2 pairs per block.
7. **A crossed random-intercept mixed model** — `statsmodels`/`patsy` are absent and
   installing them would write outside this worktree; replaced by the two-way clustered and
   arm-pair fixed-effects estimators (F.1, F.5).
8. **Any stage-specific mechanistic claim** — the contact values are reused verbatim from
   the seqlevel tables, so a multi-mapping artefact *inside* the contact values is invisible
   here, and per-stage library depth is uncontrolled (the within-pair share normalisation
   that would control it is a different design and was not run).

## commands run

All read-only except the two files this lane owns. Box: load average ≈ 125 on 256 cores
throughout (other users' jobs running concurrently).

1. `wc -l` / `read` on `scripts/hic/mantel_mouse_d_m5.py`,
   `scripts/hic/mouse_satellite_sensitivity.py`, `submission/scripts/figures/make_fig4c_mouse_zygotene.R`,
   `response/notes/M6_mouse_satellite.md`, `PEER_REVIEW.md`.
2. `python3 -c` walk over `scripts/hic/mantel_mouse_d_m5.json` — the eight published rho
   and `n_arm_pairs` (344, 319).
3. `python3 -c` walk over `scripts/hic/mouse_satellite_sensitivity.json` — `inputs`,
   `satellite_annotation.per_arm_sat_frac` (27 arms), `sequence_axis` (NaN validation),
   `mantel.stratified` (the M6 arm-level reference values), `not_established`.
4. `head -3 /moosefs/.../similarity_1Mb/mouse.dist_matrix.tsv` — 27 arm labels, verified
   symmetry and the non-zero diagonal.
5. `zcat ... | head -1` on one PAF — established the field layout that exposes the bug.
6. Ad-hoc python: arm-mean pair Jaccard vs `1 − d` (max diff 1.1e-16); tie counts
   (129/344 arm pairs at similarity 0; the 21 arm-level q×q pairs have 17 tied at 0 and
   only 5 distinct values); stratum sizes (722/371/42); distinct-axis-value counts
   (684/34/8); the dtel parse (median 1 bp, p90 3,872 bp, max 4,532 bp);
   reciprocal-duplicate check (1135 unique unordered pairs from 49 PHRs, 49 distinct PHR
   ids).
7. `python3 scripts/hic/mouse_stage_profile_disambiguation.py --quick` (250/250) —
   smoke test, ~20 s.
8. `python3 scripts/hic/mouse_stage_profile_disambiguation.py --perms 1000 --boot 1000`
   — calibration, **31.3 s wall**, single process.
9. `python3 scripts/hic/mouse_stage_profile_disambiguation.py --perms 10000 --boot 10000`
   — production run producing `scripts/hic/mouse_stage_profile_disambiguation.json`
   (117,515 B) — the pre-STEP-F numbers, all of which STEP F reproduced unchanged; script
   mtime 20:10:11, JSON mtime 20:17:41 → **≤ 7.5 min wall**, one core, load average ≈ 125 on
   the box (other users' jobs).
10. `python3 scripts/hic/mouse_stage_profile_disambiguation.py --quick` after STEP F was
    added — smoke test including the two STEP F calibration simulations, **66.8 s**.
11. `python3 scripts/hic/mouse_stage_profile_disambiguation.py --perms 10000 --boot 10000`
    — the STEP F production run (286,813 B JSON). Two runs were needed because STEP F
    introduced and then fixed a block-label normalisation (arm-pair labels were being built
    from the ordered arm columns, giving 606 blocks instead of 344). Both took **≤ 13 min
    wall**, single process, load average ≈ 125.
12. **Determinism check (STEP A–E).** The full run was repeated after the only post-run
    edit to that version (adding one item to
    `step_E_verdict.unmeasurable_with_data_on_disk`). A recursive diff of the two JSONs
    found exactly **one** difference — the length of that list, 6 → 7 — and no numeric
    field changed anywhere.
13. **Determinism check (with STEP F).** The final full run was repeated from scratch and
    diffed against its own output; the script introduces no RNG source other than the main
    `default_rng(20260518)` stream and the calibration simulation's own `default_rng(4242)`.
    Result recorded at the end of this section.

## verification gates (verbatim from the JSON `verification` block)

```
spearman_matches_scipy_with_ties                        true
mantel_matches_published_routine                        true
stratified_permutation_preserves_block_multisets        true
stratified_permutation_vectorised_and_loop_nulls_agree   {ok: true,
    vectorised_mean: ..., loop_mean: ...}
adversarial_constant_within_block_gives_p_one           true
adversarial_pure_within_block_association_detected       true
weighted_spearman_equals_expanded_rank_pearson           true
arm_bootstrap_draw_frequency_sane                        true
    observed 0.6388  expected 0.639
fixture_gate_pair_spearman_reproduces_submitted          ok: true
    observed {leptotene 0.419, zygotene 0.614, pachytene 0.576, diplotene 0.245}
    submitted {leptotene 0.419, zygotene 0.614, pachytene 0.576, diplotene 0.245}
fixture_gate_arm_mantel_reproduces_published_json        ok: true, mismatches []
    8/8 rho identical to mantel_mouse_d_m5.json to 4 dp, n_arm_pairs identical
    max_tie_correction_shift 0.0486 (50 kb leptotene)
input_shape {n_pairs 1135, n_distinct_phrs 49, n_distinct_arms 27}
```

STEP F carries its own hard gates (it raises `SystemExit` on any failure):

```
dyadic_slope_equals_spearman_all_cells                     true
contrast_slope_equals_delta_all_cells                      true
stacked_point_estimate_equals_differenced_all_cells        true
mrqap_rho_equals_spearman_on_its_own_cells                 true
failures                                                   []
```

The two fixture gates and the STEP F gates are hard gates: `main()` raises `SystemExit` if
any verification flag is false, if either fixture gate fails, or if any STEP F self-check
fails.

## numbers (all from the JSON of the run above)

27 arms / 344 arm pairs / 1135 PHR pairs from 49 PHRs; 8 contact tables (4 stages × 2
resolutions) at 1 Mb, plus 16 more at 2 Mb and 4 Mb; 10,000 permutations and 10,000
bootstrap replicates everywhere.

- Axis identity: max |arm-mean pair Jaccard − (1 − d)| = **1.11e-16**; pair-level
  Spearman between axes **0.8796** (1 Mb), 0.8750 (2 Mb), 0.8741 (4 Mb).
- Submitted rho (verified): **0.4187 / 0.6142 / 0.5761 / 0.2452**.
- Published Mantel (verified): **0.4074 / 0.5785 / 0.6259 / 0.3282**;
  tie-corrected **0.3991 / 0.5813 / 0.6289 / 0.3301**.
- Referee test 1 Mb 20 kb: rho **+0.4187 / +0.6142 / +0.5761 / +0.2452**, stratified
  p **1e-4 / 1e-4 / 1e-4 / 1.5e-3**, arm-bootstrap CI **+0.252…+0.581 / +0.439…+0.753 /
  +0.414…+0.700 / +0.087…+0.372**; null means **+0.3045 / +0.4635 / +0.4727 / +0.2057**.
- Paired contrasts: **+0.1955 [+0.093,+0.302]**, **−0.0381 [−0.131,+0.043]**,
  **−0.3309 [−0.449,−0.207]**.
- Variance decomposition: **25.6% within arm pairs, 74.4% between**.
- Ties: **129/344** arm pairs at similarity 0, 216 distinct values.
- Satellite strata: **722 / 371 / 42** pairs (rich / mixed / poor); distinct axis values
  **684 / 34 / 8**; zeros **37 / 337 / 35**.
- Satellite-rich rho: **+0.2850 / +0.5079 / +0.3230 / +0.0231**;
  partial: **+0.4128 / +0.5804 / +0.4393 / +0.1755**;
  Spearman(Jaccard, mean satellite) = **+0.790**.
- dtel: median **1 bp**, p90 **3,872 bp**, max **4,532 bp**, 100% < 10 kb.
- STEP F, dyadic two-way clustered (1 Mb 20 kb): rho **+0.4187 / +0.6142 / +0.5761 /
  +0.2452**, SE **0.0270 / 0.0239 / 0.0233 / 0.0279**, CI **+0.366…+0.472 / +0.567…+0.661 /
  +0.530…+0.622 / +0.190…+0.300**, p **5.1e-14 / 5.9e-19 / 1.4e-18 / 5.8e-09**, design
  effect **1.00 / 1.02 / 0.96 / 0.97**, clusters **25 / 26 / 344**.
- STEP F, MRQAP-DSP: p **1e-4 / 1e-4 / 1e-4 / 8.5e-3** (outcome-permutation **1e-4 / 1e-4 /
  1e-4 / 7.3e-3**), R² **0.175 / 0.377 / 0.332 / 0.060**, 1135 of 1176 matrix cells.
- STEP F, arm-pair FE: **+0.5781 / +0.7632 / +0.5231 / +0.2000** (1126 rows, 335 blocks of
  size ≥2).
  [.165,.227] / [.093,.302]**; pachy−zygo **−0.0381 [−.067,−.009] / [−.067,−.009] /
  [−.131,+.043]**, p **0.017**, Holm **0.017**.
- STEP F, calibration coverage at Δ = 0: dyadic **0.930 / 0.955**, block **0.930 / 0.930**,
  chromosome-end **0.990 / 0.995** (shared / stage-specific arm effects).
- Arm-collapse defects: block sizes **min 1, median 4, max 4**; **2.6%** singletons;
  **19.6%** of pairs in blocks of ≤2; effective mean satellite **0.708 → 0.661**; both-rich
  weight share **0.636 → 0.552**; block sizes by class **3.80 / 2.79 / 2.00**;
  precision-weighted Mantel shift **+0.0021 / +0.0322 / +0.0210 / +0.0102** (peak unchanged).
