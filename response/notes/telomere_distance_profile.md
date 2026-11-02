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

| distance (kb) | allelic | paralog | **ratio** |
|---|---|---|---|
| 0-2 | 2.260% | 2.940% | **1.30** |
| 2-4 | 0.660% | 1.960% | 2.97 |
| 4-6 | 0.850% | 1.600% | 1.88 |
| 6-8 | 0.260% | 1.290% | 4.96 |
| 8-10 | 0.580% | 1.150% | 1.98 |
| 10-12 | 0.220% | 1.270% | 5.77 |
| 12-14 | 0.290% | 1.520% | 5.24 |
| 14-16 | 0.180% | 1.000% | 5.56 |
| 16-20 | 0.18-0.22% | 0.73-0.81% | 3.3-4.5 |
| 20-50 | 0.07-0.17% | 0.47-0.94% | 3.5-7.8 |
| 50-100 | 0.11-0.13% | 0.56-0.62% | 4.8-5.1 |
| 100-200 | 0.13-0.14% | 0.60-0.68% | 4.3-5.2 |
| 200-500 | 0.150% | 1.050% | 7.00 |

## What this says

**The relative homogenisation does peak at the junction.** In the first 10 kb the
paralog/allelic ratio is **1.3-3.0**; from 20 kb inward it is **4-8** (median ~5). At the extreme
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
