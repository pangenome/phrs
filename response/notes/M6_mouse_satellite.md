# M6 — Mouse p-arm satellite confound (referee #1, major point 6)

> Provenance for this note: see `response/PROVENANCE.md` (tool versions: `response/notes/tool_versions_and_pins.md`).

Status: **annotation-level answer complete; contact-axis masking NOT done; and the
"masked sequence axis" numbers already present in the JSON are VOID (a coordinate-order bug,
diagnosed below — do not quote them).**

Lane E write-up. Numeric source: `scripts/hic/mouse_satellite_sensitivity.json` (83,514 B,
mtime Oct 23 18:57) produced by `scripts/hic/mouse_satellite_sensitivity.py` (44,837 B,
mtime Oct 23 18:54). **The analysis was not re-run for this note.** Every number below is
either read out of that JSON or produced by a command recorded in `## commands run`.

Referee #1, major point 6, verbatim (field `referee_verbatim` of the JSON):

> Mouse chromosomes are acro/telocentric with essentially no short arm, so a 1 Mb
> telomere-anchored "p-arm" flank is largely composed of minor and major satellite, which
> are near-identical between chromosomes, placing these windows at once at the top of the
> Jaccard axis and at the maximum of the multi-mapping problem.

---

## 1. What was real annotation and what was a proxy

This is stated plainly because the distinction matters for how far the answer can be pushed.

| element | status | evidence in JSON |
|---|---|---|
| major satellite | **REAL curated annotation** — Dfam 3.8 family `GSAT_MM`, 471 consensus entries, header `GSAT_MM#Satellite @Mus_genus [S:15,55]`, extracted from RepeatMasker 4.1.7p1 `RepeatMasker.lib` and searched with blastn | `satellite_annotation.families.GSAT_MM = 471`, `family_headers.GSAT_MM` |
| minor satellite | **PROXY** — Dfam 3.8 in this local library has **no** Mus family named "minor satellite"; `MurSAT1` (120) + `MMSAT4` (168) + `CENSAT_MC` (59) were substituted | `satellite_annotation.label`, `families`, `family_headers` |
| tandem periodicity | **PROXY** — composition-only detector, no annotation involved: a 5 kb bin is called satellite-like if the fraction of positions with `s[i] == s[i+lag]` exceeds 0.5 for some lag in 100–260 bp (major 234 bp, minor 120 bp tandem) while the 6 bp match fraction is below 0.8 | script docstring § (a) SOURCE 2; `per_flank.*.proxy_periodic_bin_frac`, `proxy_best_lag_median` |

Hit filters: `blastn length >= 100 bp, identity >= 75%; bin covered if >= 50%`
(`satellite_annotation.hit_filters`). Bins are 5 kb (`inputs.bin = 5000`), 78 haplotype
flanks (39 p-end, 39 q-end), 27 arms.

Two things this table has to say out loud:

* **The minor-satellite half of the referee's sentence cannot be answered from this
  analysis.** The proxy families yield almost nothing on the very windows at issue: mean
  minor-proxy fraction is **0.00004** on the 39 p-end flanks and **0.00112** on the 39 q-end
  flanks (computed from `per_flank`). Whether *minor* satellite contributes to the p-arm
  signal is therefore **unmeasured**, not measured-and-small.
* The two independent sources agree at bin level: mean `agree_bins` = **0.937** over 78
  flanks (min 0.5226, max 1.0000); 0.8993 on p-end flanks, 0.9747 on q-end flanks. So the
  periodicity proxy does track the annotation where the annotation is dense.

**Which mouse ends carry the satellite (referee asks this explicitly).** The 1 Mb
telomere-anchored flanks are strongly bimodal:

* p-end flanks: `sat_frac_union` **min 0.2597, max 0.9928**; **38 of 39** have ≥ 0.50
  satellite.
* q-end flanks: `sat_frac_union` **min 0.0021, max 0.0077** — effectively satellite-free.
* Per-arm (`satellite_annotation.per_arm_sat_frac`, 27 arms): 20 p arms range
  **0.4673 (chr7_p) – 0.9853 (chr1_p)**; the 7 q arms present in this dataset
  (chr8_q 0.0049, chr9_q 0.0077, chr10_q 0.0053, chr11_q 0.0057, chr17_q 0.0064,
  chr19_q 0.0062, chrX_q 0.0037) range **0.0037–0.0077**.
* Major-satellite-only (`per_arm_major`): p arms 0.0001 (chrX_p) – 0.6960 (chr18_p); every
  q arm is 0.0000.

So the referee's premise is quantitatively correct for this dataset: the p-arm flank *is*
mostly satellite (mean union fraction 0.8848 across the 39 p-end flanks, mean
`frac_bins_sat50` 0.8607), and the q-arm flank is not (mean 0.0044).

---

## 2. Baseline reproduction (does the re-run match the manuscript numbers?)

Arm-level Mantel, sequence axis = the published `mouse.dist_matrix.tsv` odgi Jaccard,
contact axis = `data/mouse_<stage>_phr_<res>bp_seqlevel.tsv` `hic_contact_norm`
(`mantel.baseline_reproduced`, `mantel.perms = 200`). `null95` is the 2.5–97.5 percentile
band of the 200-draw row-and-column permutation null; `boot` is the bootstrap 95% CI of rho.

| stage / res | rho | p_perm | null mean | null 95% band | boot 95% rho | n arm pairs |
|---|---|---|---|---|---|---|
| leptotene 20 kb | +0.4074 | 0.004975 | +0.0043 | −0.0895 … +0.0969 | 0.2817 … 0.4818 | 344 |
| leptotene 50 kb | +0.2882 | 0.004975 | +0.0031 | −0.0951 … +0.1183 | 0.2255 … 0.4146 | 344 |
| zygotene 20 kb | +0.5785 | 0.004975 | +0.0007 | −0.1116 … +0.0997 | 0.4714 … 0.6488 | 344 |
| zygotene 50 kb | +0.3369 | 0.004975 | −0.0011 | −0.1053 … +0.0976 | 0.1992 … 0.4157 | 344 |
| **pachytene 20 kb** | **+0.6259** | 0.004975 | +0.0043 | −0.1011 … +0.1149 | 0.5297 … 0.6846 | 344 |
| pachytene 50 kb | +0.3130 | 0.004975 | −0.0006 | −0.0933 … +0.1114 | 0.2471 … 0.4481 | 344 |
| diplotene 20 kb | +0.3282 | 0.004975 | +0.0035 | −0.1086 … +0.1052 | 0.2249 … 0.4078 | 344 |
| diplotene 50 kb | +0.3043 | 0.004975 | −0.0095 | −0.1163 … +0.1155 | 0.2180 … 0.4218 | 319 |

Reproduction is **exact on rho**: all eight values are identical to
`scripts/hic/mantel_mouse_d_m5.json` (0.4074, 0.2882, 0.5785, 0.3369, 0.6259, 0.3130,
0.3282, 0.3043) with identical `n_arm_pairs` (344, and 319 for diplotene 50 kb). The shape
the manuscript reports is reproduced: **a maximum at 20 kb in pachytene (+0.6259, next
highest zygotene +0.5785) and an essentially flat 50 kb profile (+0.2882 … +0.3369 across
all four stages)**.

Two honest caveats on the p column: this run used **200 permutations**, so the smallest
attainable p is 1/201 = **0.004975** — every baseline row sits on that floor and the
baseline p-values are uninformative beyond "p < 0.005". (D-M5 itself used 10,000
permutations and reports p = 9.999e-05.) The permutation null is centred at ≈0 with a
±0.10–0.12 band, i.e. the observed rho is far outside the null.

---

## 3. Stratified result: does the association survive in the low-satellite subset?

Stratification variable = satellite fraction of an arm pair = max over its two arms
(so a pair is "high" if **either** flank is satellite-rich). `p_perm` from 200 permutations.

**The partition is degenerate**, and this has to be said: arm-level satellite fraction is
bimodal with **nothing between 0.10 and 0.25** (arms < 0.10: exactly the 7 q arms; arms in
[0.10, 0.25): none; arms ≥ 0.25: the 20 p arms). Thresholds 0.10 and 0.25 therefore select
**the same** partition, and the threshold-0.50 low set additionally admits chr7_p
(0.4673) — 8 arms, C(8,2) = 28 pairs. Consequently:

> **The "low-satellite" subset *is* the q–q subset: 21 pairs among 7 arms
> (chr8_q, chr9_q, chr10_q, chr11_q, chr17_q, chr19_q, chrX_q).** It is not a large
> independent sample; it is 7 chromosome ends, and pairs within it are not independent.

### 3a. All-satellite split (major + minor proxy + any Dfam satellite family)

20 kb:

| stage | low (thr 0.10 / 0.25, n=21, mean sat 0.0065) | high (thr 0.10 / 0.25, n=323, mean sat 0.9139) |
|---|---|---|
| leptotene | rho **+0.5364**, p 0.0249 / 0.0050, boot −0.1256 … +0.7028 | rho +0.4143, p 0.004975, boot 0.3139 … 0.5153 |
| zygotene | rho **+0.5390**, p 0.004975 / 0.0199, boot −0.1157 … +0.7208 | rho +0.6515, p 0.004975, boot 0.5669 … 0.7280 |
| pachytene | rho **+0.4857**, p 0.004975 / 0.0199, boot −0.0849 … +0.7352 | rho +0.6394, p 0.004975, boot 0.5635 … 0.7233 |
| diplotene | rho **+0.4857**, p 0.0149 / 0.004975, boot −0.1144 … +0.7356 | rho +0.2775, p 0.004975, boot 0.1839 … 0.3831 |

| stage | low (thr 0.50, n=28, mean sat 0.1217) | high (thr 0.50, n=316, mean sat 0.9238) |
|---|---|---|
| leptotene | rho +0.5884, p 0.004975, boot 0.1563 … 0.7505 | rho +0.4421, p 0.004975 |
| zygotene | rho +0.6409, p 0.004975, boot 0.1319 … 0.7531 | rho +0.6713, p 0.004975 |
| pachytene | rho +0.6141, p 0.004975, boot 0.1508 … 0.7554 | rho +0.6595, p 0.004975 |
| diplotene | rho +0.6426, p 0.004975, boot 0.1499 … 0.7537 | rho +0.2822, p 0.004975 |

50 kb:

| stage | low (n=21, thr 0.10 / 0.25) | high (n=323) |
|---|---|---|
| leptotene | rho +0.0857, p 0.363 / 0.413 | rho +0.2585, p 0.004975 |
| zygotene | rho +0.3351, p 0.0697 / 0.0846 | rho +0.3593, p 0.004975 |
| pachytene | rho +0.2429, p 0.119 / 0.119 | rho +0.2947, p 0.004975 |
| diplotene | rho +0.1870, p 0.264 / 0.209 | rho +0.2632, p 0.004975 |

(thr 0.50, n=28: leptotene +0.2868 p 0.0746; zygotene +0.4756 p 0.004975;
pachytene +0.3848 p 0.0299; diplotene +0.3552 p 0.0249.)

### 3b. Major-satellite-only split (`GSAT_MM` alone, threshold 0.05)

| stage / res | major-low (n=27) | major-high (n=317) |
|---|---|---|
| leptotene 20 kb | rho +0.3242, p 0.0398, boot −0.1137 … 0.6442 | rho +0.3994, p 0.004975 |
| zygotene 20 kb | rho +0.3504, p 0.0299, boot −0.1502 … 0.6150 | rho +0.6349, p 0.004975 |
| pachytene 20 kb | rho +0.3559, p 0.0299, boot −0.0478 … 0.6564 | rho +0.6316, p 0.004975 |
| diplotene 20 kb | rho +0.4304, p 0.0149, boot −0.0676 … 0.6068 | rho +0.2886, p 0.004975 |
| leptotene 50 kb | rho +0.1300, p 0.289 | rho +0.2935, p 0.004975 |
| zygotene 50 kb | rho +0.2625, p 0.0647 | rho +0.3674, p 0.004975 |
| pachytene 50 kb | rho +0.2814, p 0.0846 | rho +0.3442, p 0.004975 |
| diplotene 50 kb | rho +0.1870, p 0.184 (n=21) | rho +0.2632, p 0.004975 (n=298) |

### 3c. What this answers

* **Yes — a positive sequence↔contact association survives in the satellite-poor subset.**
  At 20 kb the 21 q–q pairs give rho between **+0.4857 and +0.6426** (p between 0.004975
  and 0.0249 across stages and thresholds). So the arm-level association is *not* created
  solely by satellite-rich p-arm windows.
* **But the stage structure does not survive.** At 20 kb the low-satellite subset is flat:
  leptotene +0.5364, zygotene +0.5390, pachytene +0.4857, diplotene +0.4857. The
  **pachytene maximum of the baseline (+0.6259) is gone**, and the "peak stage" label flips
  to zygotene by a margin of 0.0026 (JSON `verdicts.peak_stage_20kb_low_sat(both<0.25) =
  "zygotene"`). Note that the same flip happens in the *high* satellite subset
  (zygotene +0.6515 vs pachytene +0.6394), so the pachytene maximum is fragile to
  stratification in general, not only to masking — the baseline gap between pachytene and
  zygotene is only 0.047 while the bootstrap intervals overlap heavily.
* At 50 kb the low-satellite subset is not significant at 0.10/0.25 (p 0.0697–0.413) and
  only marginal at 0.50 — consistent with the baseline being flat there anyway.
* The bootstrap intervals for the low-satellite rows nearly always include 0
  (e.g. −0.1256 … +0.7028). With 21 pairs drawn from 7 arms, the low-satellite estimate is
  directionally clear but not precisely known.

---

## 4. Verdict summary and the p-vs-q comparison

`verdicts` (verbatim from the JSON):

```
peak_stage_20kb_baseline            : pachytene
peak_stage_50kb_baseline            : zygotene
peak_stage_20kb_low_sat(both<0.25)  : zygotene
peak_stage_20kb_high_sat(both>=0.25): zygotene
peak_stage_20kb_masked_axis         : zygotene     <-- VOID, see §5b
```

`p_vs_q`: class sizes pp 190, pq 133, qq 21 (20 p arms, 7 q arms); the contrast tested is
mean contact(pp) − mean contact(qq), with a permutation null that shuffles the p/q arm
label over the 27 arms (200 permutations, so smallest p = 0.004975).

| stage / res | mean contact pp (n=190) | pq (n=133) | qq (n=21) | pp − qq | p_perm (arm label) |
|---|---|---|---|---|---|
| leptotene 20 kb | 0.000661 | 0.000566 | 0.000910 | −2.48e-04 | 0.2985 |
| leptotene 50 kb | 0.003181 | 0.002507 | 0.001926 | +1.26e-03 | 0.3532 |
| **zygotene 20 kb** | 0.000505 | 0.000288 | 0.000881 | **−3.75e-04** | **0.009950** |
| **zygotene 50 kb** | 0.001316 | 0.001337 | 0.002550 | **−1.24e-03** | **0.019900** |
| pachytene 20 kb | 0.000618 | 0.000362 | 0.000783 | −1.65e-04 | 0.2239 |
| pachytene 50 kb | 0.002451 | 0.002066 | 0.001812 | +6.39e-04 | 0.4975 |
| diplotene 20 kb | 0.001709 | 0.001052 | 0.000796 | +9.13e-04 | 0.1642 |
| diplotene 50 kb | 0.003843 (n=171) | 0.002872 (n=127) | 0.001872 | +1.97e-03 | 0.2438 |

The informative part: the p–p pairs do **not** show higher inter-arm contact than q–q pairs
at any stage, and the only significant arm-type contrast (zygotene, both resolutions) has
**q–q above p–p**, i.e. the opposite of what a satellite-driven, multi-mapping-inflated
p-arm signal would produce. This is a point in the paper's favour on this referee item, but
it is a mean-contact contrast, not the Mantel association.

**The `rho_masked_axis` / `p_masked_axis` columns of every `p_vs_q` row are void** for the
same reason as §5b and must not be quoted.

---

## 5. Not established

### 5a. Verbatim from the JSON (`not_established`, three items)

1. > contact axis is NOT satellite-masked here: re-masking contacts needs the
   > zuo2021_<stage>.mcool per-bin matrices (h5py/cooler absent from system python;
   > 0.3-0.5 Gpixel per stage) plus the paternal/maternal haplotype mapping
2. > no full RepeatMasker run was possible (RepeatMasker 4.1.7p1 is missing
   > PubRef.pm/FastaDB.pm on this box); only the curated Dfam satellite consensus
   > sequences were used, so non-satellite repeat content is unmeasured
3. > the published arm-level Jaccard is an odgi pangenome-graph node Jaccard and cannot be
   > recomputed with satellite bins excluded without rebuilding the pggb graph; the
   > masked-axis result uses the re-derived PAF-based Jaccard instead

Item 3 is exactly why everything in §3 is a *stratification* of arm pairs and not a
*masking* of the sequence axis: the published axis cannot be masked without rebuilding the
pggb graph, and the substitute axis the script built instead **does not work** (5b).

### 5b. NEW — the "masked sequence axis" results in the JSON are void (identity bug)

This was found while reading the JSON for this note, and it must be recorded as a
limitation, not silently dropped. Three independent signatures agree:

1. `sequence_axis.validation_vs_published_odgi_jaccard = {n: 351, spearman: NaN,
   pearson: NaN}` — the derived axis failed to validate against the published odgi Jaccard,
   and the failure is a **zero-variance** input (scipy returns NaN), not a weak correlation.
2. For all eight stage×resolution rows, `rho` (masked) is **numerically identical** to
   `rho_derived_axis_unmasked` (e.g. leptotene 20 kb −0.0542 / −0.0542; pachytene 20 kb
   +0.0001 / +0.0001), and the masked-axis rhos are ≈0 everywhere (pachytene 0 kb
   +0.0001 p=0.4328; zygotene +0.0553 p=0.1592; leptotene −0.0542 p=0.8557; diplotene
   +0.0149 p=0.3483). Identical statistics before and after masking, with values pinned
   at zero, is the signature of a constant vector — `mantel_mouse_d_m5.mantel` ranks with
   `argsort(argsort(x))`, which turns a constant x into an index ramp, so the reported rho
   is an artefact of row order in the contact table.
3. Cause: `scripts/hic/mouse_satellite_sensitivity.py:705` stores PAF alignments as
   `(int(c[1]), int(c[2]))` = `(qs, qe)` without normalising, and **all 60,766
   gi≥0.95 records in these PAFs have qs > qe**. `intersect_intervals`
   (`...py:244`) sorts the tuples but keeps the reversed pair, so
   `s = max(qs, ·)`, `e = min(qe, ·)` gives `s < e == False` and every intersection is
   empty (`...py:738`). Every derived Jaccard is therefore 0.0, the axis is constant, and
   `...py:780` yields NaN.

Bounded diagnostic (not a re-run of the analysis; 78 gzipped PAFs read, ~93k alignments):

```
records with qs>qe (reversed): 60766 / 60766
```

So: **no satellite-masked sequence axis exists yet.** §3's evidence is stratification only,
which is nevertheless the direct answer to the referee's question ("re-run the stage
analysis with satellite excluded" → we excluded the satellite-rich *windows*, not the
satellite *bases*).

### 5c. Also not established

* **Minor satellite composition** of these windows (§1) — the proxy recovers ~0 everywhere.
* **Non-satellite repeat content** of the p-arm flanks (no RepeatMasker run; 5a item 2).
  TAR1/subterminal/TE families are unmeasured here.
* **Statistical power of the low-satellite subset**: 7 arms. The p-values in §3a for that
  subset are computed as if 21 pairs were 21 independent observations; they are not.
* The `p_perm` values for identical pair sets differ between the thr 0.10 and thr 0.25 rows
  (e.g. leptotene 20 kb 0.0249 vs 0.004975 with identical rho +0.5364) purely because the
  200-permutation null is redrawn from a shared RNG stream. Monte-Carlo noise at this
  permutation count is of order 0.005–0.02; the peak-stage verdicts in §4 rest on margins
  (0.0026, 0.012, 0.047) that are of that order.

---

## 6. One-paragraph verdict

Annotation-level masking does **not** remove the mouse association, but it also does not
rescue the specific stage claim. With the satellite-rich windows set aside, the arm-level
sequence↔contact association persists in the only satellite-poor subset available — the
21 q–q arm pairs, rho +0.486 to +0.643 at 20 kb with permutation p between 0.005 and 0.025
— so the association is not an artefact of satellite-rich p-arm flanks alone, and
separately the p–p arm pairs do not show higher contact than q–q pairs (the only
significant arm-type contrast, at zygotene, runs the other way). What does not survive is
the **stage ordering**: in the satellite-poor subset the 20 kb profile flattens
(0.5364 / 0.5390 / 0.4857 / 0.4857 for leptotene / zygotene / pachytene / diplotene) and
the pachytene maximum of the baseline (+0.6259) disappears; the same flip occurs in the
satellite-rich subset, so the "pachytene peak" is fragile to stratification generally and
should not be stated as a satellite-independent result. Three things would settle the
question and none were available here: (i) contacts with satellite bins masked from the
`zuo2021_<stage>.mcool` per-bin matrices together with the paternal/maternal haplotype
mapping, which requires h5py/cooler and 0.3–0.5 Gpixel per stage; (ii) a real RepeatMasker
run (this box's 4.1.7p1 is missing `PubRef.pm`/`FastaDB.pm`), so that non-satellite repeats
are excluded too; (iii) rebuilding the pggb graph so the published odgi graph-node Jaccard
can itself be recomputed with satellite nodes removed — until then the masked-axis numbers
would have to come from a re-derived PAF Jaccard, and the one attempt at that in this JSON
is void because of the qs>qe coordinate bug in §5b. Until (i)–(iii) are done, the honest
claim is: the association is not explained away by satellite content, but the stage
attribution is.

---

## commands run

All read-only. Nothing was re-run; `/moosefs` reads succeeded (no I/O errors encountered).

1. `ls -la scripts/hic/mouse_satellite_sensitivity.*` — sizes/mtimes quoted in the header.
2. `grep -n -i -A18 "satellite" PEER_REVIEW.md` — retrieved the referee's point 6 text;
   confirmed authority file used (the `paper_prep/synthesis/SIMULATED_PEER_REVIEW_v1.md`
   file was never opened).
3. `python3 -c "import json; ..."` walks over `mouse_satellite_sensitivity.json` — dumped
   `inputs`, `satellite_annotation`, `sequence_axis`, `mantel`, `p_vs_q`, `verdicts`,
   `not_established` and the per-arm aggregates (all tables above).
4. `head -3 /moosefs/.../similarity_1Mb/mouse.dist_matrix.tsv` and
   `head -3 .../subtelo_1Mb/mouse.all-vs-all.p95.id95.len.tsv` — confirmed the 27-arm
   published axis and the PHR region columns; 29 of 78 flanks are PHR-free (`.`).
5. `ls /moosefs/.../subtelo_1Mb/all-vs-all.p95/ | wc -l` → 78 PAF files.
6. `zcat .../mouse.all-vs-B6_1_chr10_parm_8494_1000000.paf.gz | head -1` — PAF field/name
   layout (`gi:f:` tag present, names identical to FASTA/PHR keys).
7. Bounded PAF diagnostic (heredoc, stdlib only, no file written, 78 gz files, 93,210
   alignments): counted `qs > qe` reversed records → `60766 / 60766`; recomputed the
   derived Jaccard for all 1,499 directed pairs → `n_zero: 1499, n_distinct: 1`.
8. `python3 - <<EOF` aggregation over `per_flank` — p-end vs q-end means
   (sat 0.8848 vs 0.0044; minor proxy 0.00004 vs 0.00112; `agree_bins` 0.8993 vs 0.9747),
   flank ranges, and the arm-threshold partition check
   (`arms<0.10: 7 q arms`, `arms in [0.10,0.25): none`, `arms>=0.25: 20`).
9. `ls -la data/mouse_*_phr_*bp_seqlevel.tsv` + `head -2` of the pachytene 20 kb table —
   confirmed the 8 contact tables exist and the `hic_contact_norm` column the Mantel uses.
10. `grep -n` on `scripts/hic/mouse_satellite_sensitivity.py` and
    `scripts/hic/mantel_mouse_d_m5.py` — located `...:705`, `...:738`, `...:244`, `...:780`
    and the `argsort` ranking in `mantel()`.

## numbers

See the tables in §1–§4. Key figures (all from the JSON or the commands above):
27 arms (20 p, 7 q); 344 arm pairs (319 for diplotene 50 kb); 78 flanks (39 p / 39 q);
5 kb bins; 200 permutations (p floor 0.004975); baseline rho pachytene 20 kb +0.6259,
50 kb range +0.2882…+0.3369 (exact match to `mantel_mouse_d_m5.json`); low-satellite subset
rho +0.4857…+0.6426 at 20 kb (n=21 pairs / 7 arms); p-end satellite fraction 0.2597–0.9928
per flank (mean 0.8848), q-end 0.0021–0.0077 (mean 0.0044); Dfam `GSAT_MM` = 471 consensus
entries, `MurSAT1` 120 / `MMSAT4` 168 / `CENSAT_MC` 59 for the minor proxy; blastn
filters ≥100 bp / ≥75% identity / ≥50% bin coverage; `agree_bins` mean 0.937;
60,766/60,766 PAF records with qs>qe; validation spearman/pearson NaN.

## what this does NOT establish

* No satellite-**masked** contact axis (no mcool access) and no satellite-masked published
  Jaccard (would need a pggb rebuild); §3 is a window *stratification*.
* The re-derived PAF-based Jaccard that was supposed to carry the masked axis is void
  (constant zero vector; qs>qe bug). Nothing should be quoted from
  `mantel.stratified.*_axis_maskedSat`, `p_vs_q.*.rho_masked_axis`, or
  `verdicts.peak_stage_20kb_masked_axis`.
* Minor-satellite content is a proxy that recovers ≈0 on these windows; "minor satellite
  does not contribute" is **not** supported — it is unmeasured.
* Non-satellite repeat content (TAR1, subterminal, shared TE families) is unmeasured.
* The low-satellite subset is 7 chromosome ends; its p-values treat 21 non-independent
  pairs as independent, and its bootstrap CIs mostly straddle zero.
* Nothing here re-tests the contact data itself; the contact axis is reused verbatim from
  `data/mouse_<stage>_phr_<res>bp_seqlevel.tsv`, so a multi-mapping artefact **inside** the
  contact values would be invisible to this analysis.
* No conclusion is offered about the human maps, which are somatic and not addressed here.
