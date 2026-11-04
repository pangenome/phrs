# Paralog / allelic divergence as a function of distance from the telomere

Script `scripts/ci/telomere_distance_profile.py`; figure `telomere_profile.pdf`/`.png`
(`scripts/ci/make_telomere_profile_figure.py`). Full cohort: 18,827 window PAFs, 465 assemblies,
2 kb bins near the junction.

**Zero point.** The telomere-subtelomere junction, i.e. the end of each window at the chromosome
terminus: p arms at sequence position 0, q arms at position L (window intervals run to the contig
end; windows are stored in original contig orientation).

**Rules.** A record is assigned only when the midpoint of *both* the query and target spans fall in
the same 2 kb distance bin, so "10-12 kb" means both sides lie 10-12 kb from their own telomeres.
Paralogs are inter-chromosomal only (same-chromosome pairs excluded). Divergence is substitution-only
(`X/(=+X)` from CIGAR), allelic and paralog computed identically.

| distance (kb) | allelic median | paralog median | ratio | allelic bp (Gbp) | paralog bp (Gbp) |
|---|---|---|---|---|---|
| 0-2 | 2.260% | 2.940% | **1.30** | 0.4 | 4.0 |
| 2-4 | 0.660% | 1.960% | **2.97** | 2.4 | 11.5 |
| 4-6 | 0.850% | 1.600% | **1.88** | 6.6 | 158.5 |
| 6-8 | 0.260% | 1.290% | **4.96** | 5.2 | 34.0 |
| 8-10 | 0.580% | 1.150% | **1.98** | 5.5 | 76.3 |
| 10-12 | 0.220% | 1.270% | **5.77** | 6.5 | 42.9 |
| 12-14 | 0.290% | 1.520% | **5.24** | 5.7 | 66.5 |
| 14-16 | 0.180% | 1.000% | **5.56** | 9.5 | 75.1 |
| 16-18 | 0.180% | 0.810% | **4.50** | 10.1 | 20.2 |
| 18-20 | 0.220% | 0.730% | **3.32** | 13.2 | 22.4 |
| 20-22 | 0.120% | 0.930% | **7.75** | 8.1 | 14.7 |
| 22-24 | 0.160% | 0.940% | **5.88** | 7.4 | 13.5 |
| 24-26 | 0.100% | 0.490% | **4.90** | 60.6 | 50.9 |
| 26-28 | 0.170% | 0.600% | **3.53** | 9.1 | 15.1 |
| 28-30 | 0.140% | 0.620% | **4.43** | 6.6 | 7.9 |
| 30-32 | 0.110% | 0.520% | **4.73** | 5.8 | 11.5 |
| 32-34 | 0.120% | 0.540% | **4.50** | 6.1 | 7.7 |
| 34-36 | 0.090% | 0.590% | **6.56** | 6.4 | 4.6 |
| 36-38 | 0.100% | 0.570% | **5.70** | 3.5 | 5.8 |
| 38-40 | 0.090% | 0.620% | **6.89** | 5.1 | 9.5 |
| 40-42 | 0.090% | 0.490% | **5.44** | 4.5 | 4.2 |
| 42-44 | 0.100% | 0.600% | **6.00** | 5.9 | 4.3 |
| 44-46 | 0.110% | 0.470% | **4.27** | 5.2 | 6.3 |
| 46-48 | 0.110% | 0.610% | **5.55** | 4.8 | 4.9 |
| 48-50 | 0.070% | 0.540% | **7.71** | 6.1 | 3.6 |
| 50-75 | 0.130% | 0.620% | **4.77** | 125.6 | 206.8 |
| 75-100 | 0.110% | 0.560% | **5.09** | 165.9 | 217.2 |
| 100-150 | 0.140% | 0.600% | **4.29** | 319.3 | 456.6 |
| 150-200 | 0.130% | 0.680% | **5.23** | 312.2 | 176.0 |
| 200-500 | 0.150% | 1.050% | **7.00** | 2219.7 | 137.4 |

## What this says

**The relative homogenisation does peak at the junction, but as a step rather than a gradient.**
In the first 10 kb the paralog/allelic ratio is **1.3-3.0**; it then jumps to ~5 across 10-16 kb,
falls back to **3.32 at 18-20 kb**, and only settles into the 4-8 band (median ~5) beyond 20 kb.
The 12-20 kb band is where the allelic estimate first becomes stable (0.18-0.29% versus 0.26-2.26%
closer in), so the transition is best described as two breakpoints at roughly 10 kb and 20 kb, not a
monotonic decay. Note also that 12-14 kb carries the **highest paralog divergence in the whole
0-50 kb range (1.52%)**, which is not what a simple monotonic homogenisation model predicts. At the extreme
tip, non-homologous copies are within 30% of the divergence measured between homologous copies -
the strongest form of the inter-arm exchange statement, and the shape a crossover/conversion model
predicts.

**But the absolute values move the other way, and that matters.** Both classes are *most* divergent
at the junction (allelic 2.26%, paralog 2.94%) and *least* in the interior (0.09-0.15% and
0.5-0.6%). So the tip is not "cleaner"; it is a regime where alleles and paralogs converge because
both are inflated. Two candidate explanations, not yet separated:

1. **Genuine exchange** - terminal repeats are the most actively exchanged sequence in the genome, so
   every copy there is a recent mosaic of every other copy, and there is little information left to
   distinguish a homolog from a paralog.
2. **Repeat-array convergence / alignment ambiguity** - terminal arrays are shared repeat families
   with poor anchorability, so allelic comparisons there partly align paralogous copies, inflating
   the denominator and driving the ratio toward 1 without any exchange.

**The control that separates them** (not yet run): recompute the profile using only long
(>= 5 kb), high-identity (>= 95%) alignment records, so the allelic estimate cannot be built from
short repeat matches; and additionally restrict to records where the query and target alignments span
comparable lengths. If the tip ratio stays near 1.3-2 with those records, explanation 1 holds; if it
rises to the interior value, explanation 2 does.

**Second caveat.** The 200-500 kb bin reverses (ratio 7.0, allelic 0.15%, paralog 1.05%) on a much
larger allelic base (2.2 Tbp) - that region is interstitial-adjacent sequence dominated by ordinary
segmental duplications, not the subtelomere tip, so it should not be read as part of the gradient.

---

## CORRECTION (supersedes the "step structure" reading above)

The per-arm, equal-weight analysis (`telomere_end_ratio_curves.py`, ED Fig 6) changes the shape, and the
earlier description in this note was a **pooling artefact**.

**What the two aggregations do:**

- **Pooled** (this note, all arms weighted by aligned base pairs): ratio 4.96 at 6-8 kb, dipping to 1.98
  at 8-10 kb, 5.77 at 10-12 kb, 3.32 at 18-20 kb, then 4-8 in the interior. From that I inferred steps
  with breakpoints near 10 and 20 kb.
- **Per-arm, one vote per chromosome** (ED Fig 6a, histogram medians within each arm, then median across
  48 arms): ratio **~1.3 in the first 8 kb**, rising continuously through ~2-3 at 20-60 kb, ~3 at
  100-200 kb, and **~6-7 beyond 300 kb**. No steps, no dip.

**The per-arm version is the one to report**, because a chromosome end is the biological unit and pooling
lets a handful of high-record arms drive the curve. So the correct description is a **monotone increase of
the paralog/allelic ratio with distance from the telomere, from ~1.3 within the first 8 kb to ~6-7 in the
far interior** - not a step function, and not a plateau.

**Consequence.** The original prediction - homogenisation strongest nearest the telomere - is **supported
in the measurable range**, and more strongly than the pooled analysis suggested: within the first ~8 kb
(immediately outside the unmeasurable zone) paralogs are only ~1.3x more divergent than alleles, i.e.
nearly interchangeable with them. My earlier "step, not gradient" correction was itself the artefact.

**Unchanged:** the terminal ~6 kb remains unmeasurable (allelic denominator fails), so the gradient cannot
be extrapolated to the telomere itself, and "ratio approaches 1 at the junction" still cannot be claimed
from data - only "ratio is ~1.3 at 8 kb and rises inward".

## Why the terminal ~6 kb is unmeasurable: identified mechanism (not repeats, not the array)

An earlier revision attributed the excluded zone loosely to "repeat arrays". Two direct measurements rule
both obvious explanations out, and a third supplies the real one.

**(i) The telomeric arrays are not in these windows.** Median (TTAGGG)+(CCCTAA) unit count per kb over 40
sampled windows, oriented telomere-first:

| 0-2 kb | 2-6 kb | 6-10 kb | 10-20 kb | 20-40 kb |
|---|---|---|---|---|
| 1.0 | 0.5 | 0.2 | 0.5 | 0.6 |

A genuine array runs ~150-170 units per kb. The windows were trimmed at the array (hence
`telo_500kb_trimmed`), so tip divergence is subtelomeric sequence, not telomere.

**(ii) It is not a low-complexity zone - the opposite.** Over 30 sampled windows:

| distance | N% | low-complexity % (6-mers repeated within segment) | GC% |
|---|---|---|---|
| 0-2 kb | 0.00 | 32.8 | **62.3** |
| 2-6 kb | 0.00 | 28.4 | 41.6 |
| 6-10 kb | 0.00 | 31.8 | 39.5 |
| 10-20 kb | 0.00 | 51.0 | 41.0 |
| 20-50 kb | 0.00 | **82.4** | 44.7 |

Repetitiveness rises *inward* (82% at 20-50 kb, where satellite and duplicon arrays sit). The first 6 kb is
the least repetitive part of the subtelomere, with an unusually GC-rich first 2 kb.

**(iii) What actually fails is allelic alignability.** For two CHM13-target windows, same-arm (allelic)
alignments by distance on the target:

| distance (kb) | chr9_q records | median block (bp) | chr9_p records | median block (bp) |
|---|---|---|---|---|
| 2-4 | **60** | **2,728** | **41** | **3,007** |
| 8-10 | 92 | 14,780 | 29 | 13,004 |
| 24-26 | 82 | 12,417 | 245 | 50,033 |
| 28-30 | **4,933** | **49,988** | **5,283** | **49,983** |

Within the first few kb, **two orders of magnitude fewer individuals align, in blocks ~18x shorter**. So the
elevated "allelic" divergence in that zone is estimated from a small, haplotype-biased set of short
fragments - it is not orthologous divergence, and the ratio built on it has no interpretation.

**Correct statement:** the terminal few kb carry haplotype-specific structural variation to the point that
same-arm alignments are sparse and fragmentary, which invalidates the allelic denominator. This matches
Tan et al.'s independent finding that subtelomeric haplotype differences lie largely at the telomeric end.
The exclusion is therefore a consequence of structural polymorphism, not of repetitiveness and not of the
telomere array.
