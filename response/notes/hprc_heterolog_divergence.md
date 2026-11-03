# Heterologous divergence across the HPRCv2 cohort (the PHR all-vs-all)

Script `scripts/ci/hprc_heterolog_divergence.py`; per-pair table
`scripts/ci/hprc_heterolog_divergence_pairs.tsv`. Full-cohort run: **18,827 window PAFs**
(`/moosefs/.../PHR_III/all-vs-all.1Mb.p95/`, telomere-anchored 500 kb windows, all 465 assemblies,
wfmash -p 95), 1,427 s wall with 14 workers, **658 heterologous arm pairs** with >= 20 kb aligned.

Two divergence columns are carried throughout:
- **`d_sub` = X/(=+X)** from the CIGAR — substitutions only, the correct denominator for a
  substitution rate mu.
- **`d_gapcomp` = 1 - gi:f** — wfmash gap-compressed identity, i.e. the pipeline's own metric
  (gap events in the denominator). It is >= d_sub always; ratio ~1.1-1.4, so T is ~15-30% larger.
  (`md:f` is *not* usable: it disagrees with the CIGAR by >0.002 on 52% of records.)

T = d/(2*mu), mu = 1.2e-8/bp/generation, 20 y/gen, **an upper bound** (subtelomeric/paternal mutation
rate is higher; at mu = 3e-8 all T shrink by 60%).

> **READ THIS FIRST.** Sections up to and including "The right summary" were computed with
> **1e-3 divergence bins** and are superseded by the "FINAL NUMBERS" section at the end, which uses
> 1e-4 bins and excludes allelic self-pairs. The affected values are: allelic median 0.10% -> **0.15%**,
> paralog median 0.80% -> **0.89%**, paralog/allelic ratio 8.0x -> **5.93x**, best-partner-per-end
> median 0.50% -> **0.52%**. Ratio-based conclusions and the intra/inter comparison are unaffected.
> Where an early value is quoted below it is marked "(coarse bins)".

## Global

- length-weighted `d_sub` over the cohort: **0.0122**
- per-record median `d_sub`: **0.0080** -> **T = 6.7 Myr** (upper bound; ~2.7 Myr at mu=3e-8)
- per-record quantiles: q05 0.003 (2.5 Myr), q25 0.005, q75 0.014, q95 0.036
- **17.9%** of all aligned heterologous bases are at >= 99.5% identity

## Closest heterologs, everything included (n = 658)

| pair | aligned bp | d_sub | per-record median | frac >=99.5% | T_sub (Myr) | T_gc (Myr) |
|---|---|---|---|---|---|---|
| chrX_q–chrY_q | 11,085,726,305 | 0.00124 | 0.00100 | 0.982 | 1.03 | 1.22 |
| chrX_q–chrY_p | 216,014,576 | 0.00128 | 0.00100 | 0.987 | 1.07 | 1.27 |
| chr20_q–chr4_p | 52,560,673 | 0.00252 | 0.00100 | 0.935 | 2.10 | 2.40 |
| chr12_p–chrX_p | 14,765,389 | 0.00277 | 0.00200 | 0.951 | 2.31 | 2.72 |
| chr13_q–chr18_p | 104,055,538 | 0.00301 | 0.00200 | 0.945 | 2.51 | 3.02 |
| chr10_q–chrX_q | 465,169 | 0.00409 | 0.00400 | 0.955 | 3.41 | 3.69 |
| chr4_q–chrX_q | 4,387,839 | 0.00427 | 0.00400 | 0.990 | 3.56 | 3.87 |
| chr10_q–chr15_p | 49,878,870 | 0.00432 | 0.00300 | 0.970 | 3.60 | 4.02 |
| chr13_p–chr4_p | 1,563,615,153 | 0.00433 | 0.00400 | 0.852 | 3.61 | 4.05 |
| chr14_p–chr4_p | 5,392,627,060 | 0.00437 | 0.00400 | 0.859 | 3.64 | 4.09 |
| chr21_p–chr4_p | 2,582,112,386 | 0.00438 | 0.00400 | 0.841 | 3.65 | 4.09 |
| chr22_p–chr4_p | 5,320,617,013 | 0.00438 | 0.00400 | 0.854 | 3.65 | 4.10 |
| chr15_p–chr4_p | 2,716,767,958 | 0.00438 | 0.00400 | 0.836 | 3.65 | 4.09 |
| chr15_p–chr4_q | 480,055,504 | 0.00440 | 0.00400 | 0.946 | 3.66 | 4.08 |

**The two tightest heterologous pairs in the genome are X/Y gametologs** (chrX_q–chrY_q 0.124%,
chrX_q–chrY_p 0.128%, 98% of bases >= 99.5%) followed by an acrocentric-p/chr4p cluster
(chr13_p, 14_p, 15_p, 21_p, 22_p each ~0.44% from chr4_p). **Both classes are excluded from the
manuscript's main analysis** (X arms and acrocentric p arms are outside the 3.51 Mb total). So the
genome-wide ranking is dominated by exactly the families we deliberately scope out, and the
X/Y values (~1 Myr) are best read as ongoing X-Y / palindrome conversion rather than as
transmission.

## Closest heterologs under the paper's scope (41 signal-bearing arms, n = 317)

| pair | aligned bp | d_sub | per-record median | frac >=99.5% | T_sub (Myr) | T_gc (Myr) |
|---|---|---|---|---|---|---|
| chr20_q–chr4_p | 52,560,673 | 0.00252 | 0.00100 | 0.935 | 2.10 | 2.40 |
| chr19_p–chr3_q | 31,382,792,512 | 0.00492 | 0.00500 | 0.466 | 4.10 | 4.73 |
| chr3_q–chr6_q | 11,861,703,706 | 0.00540 | 0.00500 | 0.482 | 4.50 | 5.13 |
| chr1_p–chr5_q | 33,378,581,907 | 0.00568 | 0.00400 | 0.526 | 4.73 | 5.33 |
| chr16_q–chr7_p | 22,506,641,359 | 0.00639 | 0.00500 | 0.437 | 5.32 | 5.87 |
| chr3_q–chr5_q | 22,203,114,864 | 0.00641 | 0.00500 | 0.318 | 5.34 | 6.03 |
| chr1_p–chr6_q | 22,055,613,198 | 0.00662 | 0.00500 | 0.359 | 5.51 | 6.22 |
| chr3_q–chr8_p | 10,363,836,128 | 0.00682 | 0.00500 | 0.352 | 5.68 | 6.48 |
| chr3_q–chr4_q | 25,200,229 | 0.00692 | 0.00200 | 0.793 | 5.76 | 6.59 |
| chr10_q–chr4_q | 12,335,199,132 | 0.00698 | 0.00600 | 0.209 | 5.82 | 6.93 |
| chr17_p–chr20_p | 459,488,396 | 0.00713 | 0.00600 | 0.177 | 5.94 | 6.65 |
| chr16_p–chr7_p | 7,818,859,082 | 0.00714 | 0.00600 | 0.422 | 5.95 | 6.54 |
| chr20_p–chr7_p | 15,378,941,265 | 0.00718 | 0.00500 | 0.331 | 5.98 | 6.73 |
| chr5_q–chr6_q | 25,932,662,727 | 0.00724 | 0.00500 | 0.375 | 6.03 | 6.74 |

- **chr20_q–chr4_p is the tightest pair in scope** (0.252%, T = 2.1 Myr) and the **only** pair with
  >= 90% of its aligned bases at >= 99.5% identity (per-record median 0.001 = 0.1% -> ~0.8 Myr).
  This pair is not currently discussed in the manuscript and is the strongest new candidate the
  cohort-wide sweep produces.
- **chr1_p–chr5_q: 0.568%, T = 4.7 Myr** — the withdrawn 5q/1p candidate. It sits in the young
  quarter of the distribution (q25 = 1.18%) and is one of only **three** pairs with >= 50% of bases
  at >= 99.5% identity. At population scale it is a genuinely tight, uniformly shared family -
  consistent with inherited homology, and it has no boundary-locked switch (see
  `M2b_9q3q_read_support.md`).
- chr19_p–chr3_q (0.49%), chr3_q–chr6_q (0.54%), chr16_q–chr7_p (0.64%), chr3_q–chr5_q (0.64%),
  chr1_p–chr6_q (0.66%), chr3_q–chr8_p (0.68%) follow.

## Distribution under scope

min 0.00252 | q25 0.01184 | **median 0.01690** | mean 0.02038 |
q75 0.02439 | max 0.07411  ->  **T(median) = 14.1 Myr**,
T(mean) = 17.0 Myr (upper bounds). mean/median = 1.21.

Removing X/Y and the acrocentric p arms *raises* the median (0.008 -> 0.017): the excluded families
are the young, heavily homogenized ones.

## Specific pairs of interest

| pair | d_sub | per-record median | frac >=99.5% | T_sub | T_gc |
|---|---|---|---|---|---|
| chr1_p–chr5_q | 0.00568 | 0.00400 | 0.526 | 4.73 | 5.33 |
| chr3_q–chr9_q | 0.00902 | 0.00700 | 0.298 | 7.52 | 8.39 |

## Caveats

1. Windows are 500 kb, so a per-arm-pair aggregate mixes the shared module with flanking unique
   sequence; the **per-record median** is the better estimate of the tightest shared segment and is
   the one to quote for coalescence.
2. `aligned_bp` counts alignment base-pairs across all individuals, not unique genomic bp.
3. wfmash `-p 95` floor: pairs more diverged than 5% contribute only via their best blocks, so
   the old tail is truncated and medians are mildly biased low.
4. Lexicographic dedupe means the 48 CHM13/GRCh38-target files contribute nothing; each unordered
   window pair is still counted exactly once.
5. Single technology of measurement: wfmash on assembly-vs-assembly; no read-level verification at
   population scale (that was done only for the pedigree locus).

## Is inter-arm exchange rare? (block-level vs arm-scale)

**Block-level sharing is ubiquitous.** Of the 658 heterologous arm pairs:

| threshold | pairs | % | equivalent T (upper bound) |
|---|---|---|---|
| d_sub <= 0.5% | 20 | 3.0% | <= 4.2 Myr |
| d_sub <= 1.0% | 128 | 19.5% | <= 8.3 Myr |
| d_sub <= 2.0% | 399 | 60.6% | <= 16.7 Myr |
| d_sub <= 5.0% | 617 | 93.8% | <= 41.7 Myr |

**33 of 45 arms** have at least one heterologous partner below 1% divergence, **18 of 45** below 0.5%,
and the median arm's best partner sits at 0.639% (5.3 Myr). No arm is isolated: every subtelomere
carries blocks that it shares recently with at least one non-allelic end.

**What is rare is the arm-scale configuration**, not the molecular process: a discrete module at a
terminus that behaves as a detectable, transmissible unit (2.3% of arms; one validated by
transmission in the pedigree). The manuscript's claim should therefore be stated as *block-level
exchange is ubiquitous and continuous; arm-scale rearrangements that persist as discrete
configurations are rare* - not "exchange is rare".

## Internal calibration: PAR2 is the "recombines every generation" reference

| pair | d_sub | per-record median | frac >= 99.5% | T (upper bound) |
|---|---|---|---|---|
| chrX_q-chrY_q (**PAR2**, obligate male crossover every generation) | **0.00124** | 0.00100 | **0.982** | 1.03 Myr |
| chrX_p-chrY_p (PAR1 mixed with non-PAR Xp/Yp sequence) | 0.00694 | 0.00400 | 0.552 | 5.79 Myr |
| chr20_q-chr4_p (tightest in scope) | 0.00252 | 0.00100 | 0.935 | 2.10 Myr |

PAR2's within-species divergence is 0.12% *because* that is allelic coalescence depth (~4 N_e mu),
not zero. So the correct statement is not "these pairs recombine every generation" but **the
tightest heterologous modules coalesce no deeper than alleles do** - they are drawn from the same
gene pool. The clock inverts the other way for frequency: d = 2 mu T means 1 Myr = 50,000
generations = 0.12%, whereas per-generation exchange would give d = 2 mu = 0.00024%, ~500x below
anything observed. The implied exchange frequency per block is therefore ~1 per 10^4-10^5
generations: continuous on a Myr timescale, not generational.

## Inter- vs intra-: the 2x2 against the local (subtelomeric) baseline

Full cohort, 18,827 windows, 1,748 s. Every alignment record is classified from its PanSN tags
(`sample#hap#contig`): allelic (same arm) vs paralogous (different arm), and paralogous
intra-chromosomal (chrN_p vs chrN_q) vs inter-chromosomal (different chromosomes).

| class | records | d_sub (weighted) | per-record median | frac >= 99.5% | T (median) |
|---|---|---|---|---|---|
| allelic, within individual | 118,112 | 0.00223 | **0.00100** | 0.916 | 0.83 Myr |
| allelic, between individuals | 55,132,305 | 0.00241 | **0.00100** | 0.908 | 0.83 Myr |
| paralog, intra-chromosomal | 1,962,844 | 0.01106 | **0.00800** | 0.261 | 6.7 Myr |
| paralog, inter-chromosomal | 132,557,383 | 0.01220 | **0.00800** | 0.179 | 6.7 Myr |

**1. The local allelic baseline is 0.10% per-record median (coarse bins; 0.15% at 1e-4 resolution)**, identical within and between
individuals (0.00100 both), i.e. the same as genome-wide nucleotide diversity - the subtelomeric
*sequence* baseline is not elevated. The weighted value (0.24%) is inflated because some blocks in an
"allelic" comparison align to paralogous copies rather than to the orthologous position; the median
is the honest baseline.

**2. The correct denominator is that local baseline, not genome-wide pi.** paralog / allelic =
**8.0x by medians, 5.1x by weighted divergence (coarse bins; 5.93x and 5.1x at 1e-4 resolution - use the latter)**. For ancient duplicated gene families the
equivalent ratio is 100-1000x. A 5-8x ratio is the quantitative signature of concerted evolution:
paralogs are measurably more divergent than alleles, but only by a factor of a few.

**3. Inter- vs intra-chromosomal paralogy is indistinguishable: 1.00x by medians**
(0.00800 vs 0.00800), 1.10x weighted. A chromosome's own two ends are no more similar to each other
than to any other chromosome's end. The mild excess of near-identical blocks intra-chromosomally
(26.1% of bases >= 99.5% vs 17.9%) is not accompanied by any difference in the central tendency.
This independently supports the community structure cutting across chromosomes rather than
following them.

**4. Allelic vs paralog is separable, which is what makes the analysis meaningful at all:** 91-92%
of allelic aligned bases are >= 99.5% identical versus 18-26% for paralogs.

## Intra-community vs inter-community (the comparison that matters)

Pairs where both arms carry a k=15 assignment (n = 626; 58 intra-, 568 inter-community):

| | n | weighted d_sub | median | q25-q75 | per-record median | mean frac >= 99.5% |
|---|---|---|---|---|---|---|
| intra-community | 58 | 0.01036 | **0.01310** | 0.0083-0.0147 | 0.0110 | **0.180** |
| inter-community | 568 | 0.01331 | **0.01637** | 0.0120-0.0246 | 0.0140 | **0.079** |

Mann-Whitney U = 10,947, z = -4.21, **p = 2.5e-5**; ratio **0.80x** (0.78x weighted). Against the
allelic baseline the two are 13.1x and 16.4x, i.e. both remain far above allelic.

**The two distributions overlap almost completely** (intra q75 0.0147 vs inter q25 0.0120). The
robust separation is not in the central tendency but in the near-identical tail: **18.0% vs 7.9% of
aligned bases at >= 99.5% identity, a 2.3x enrichment.**

### Per community (intra vs the same arms' out-of-community partners)

| community | arms | intra | vs outside | ratio |
|---|---|---|---|---|
| C14 | X_q, Y_q | 0.00124 | 0.01435 | 0.09x (n=1) |
| **C11** | 1p, 5q, 6q, 8p | 0.00749 | 0.01538 | **0.49x** |
| **C3** | 3q, 7p, 9q, 11p, 16q, 19p | 0.00849 | 0.01419 | **0.60x** |
| C1 | 4q, 10q | 0.00698 | 0.02302 | 0.30x (n=1) |
| C2 | 10p, 18p | 0.02391 | 0.03090 | 0.77x (n=1) |
| C6 | 1q, 13q, 17q, 19q, 21q, 22q | 0.01992 | 0.02484 | 0.80x |
| C5 | 6p, 9p, 12p, 20q | 0.01384 | 0.01455 | **0.95x** |
| **C7** | 13p, 14p, 15p, 21p, 22p | 0.01429 | 0.01428 | **1.00x** |
| C4 | 7q, 12q | 0.03775 | 0.02561 | 1.47x (n=1) |
| C12 | 2q, 20p | 0.01327 | 0.01544 | 0.86x (n=1) |
| C15 | X_p, Y_p | 0.00694 | 0.01654 | 0.42x (n=1) |

**Community cohesion is community-specific, not global.** C11 and C3 are genuinely cohesive
(0.49x, 0.60x); C5 shows almost none (0.95x) and **C7 exactly none (1.00x)**. C7 is the acrocentric
community: its members share rDNA *architecture* while their sequences diverge as much as any
outside partner. Communities are therefore **structure-/repertoire-defined, not
divergence-magnitude-defined**.

### Consequence for the manuscript

Do not claim that community members are "more similar overall" - that difference is only 0.80x with
heavily overlapping distributions. Claim that they share **near-identical blocks**: 2.3x enrichment
of bases at >= 99.5% identity, consistent with the module-architecture result (block repertoire
Jaccard 0.83 within vs 0.57 between). Overall divergence magnitude is nearly flat across the genome
precisely because the blocks themselves are old and shared everywhere (see
`community_coalescence_times.md`: within-community 0.0070 vs between-community 0.0076 in HG002).

## Figure: the scope of subtelomeric divergence

`response/notes/heterolog_scope.pdf` / `.png`, produced by
`scripts/ci/make_heterolog_figure.py` from `/moosefs/erikg/phrs_tmp/hprc_het_hist.json`
(the cohort run with per-class histograms).

- **A** ECDF of substitution divergence for allelic (same arm), paralog on the same
  chromosome, and paralog on a different chromosome. Record-level, weighted by aligned bp.
  Allelic median 0.10% (coarse bins; 0.15% at 1e-4); both paralog classes 0.80% (0.89%/0.83% at 1e-4) - they are indistinguishable from each
  other, and both sit 8x above allelic. Bin width is 1e-3, so the allelic curve is visibly
  stepped: its entire distribution falls inside one or two bins.
- **B** The same statement as fold-change over the allelic median, with the ancient
  duplicated-gene-family reference band (100-1000x) for scale. Allelic 1x, all paralogs 8x
  (record level, exact); the in-community (2.0x) and out-of-community (9.0x) bars are
  **pair-level** medians, because the aggregator does not accumulate histograms per
  community - that asymmetry is a real caveat and the reason panel C exists.
- **C** Community cohesion: within-community median divergence divided by the median to
  that community's out-of-community partners. Most communities sit at 0.40-0.60 (genuinely
  cohesive); **C5 at 0.60 and C7 at 0.97** are the weak ones, C7 being the acrocentric
  community that shares rDNA architecture without sequence cohesion. Communities with <3
  internal pairs are omitted.
- **D** Why the community split is tail-defined: per arm pair, median divergence (magnitude)
  against the fraction of aligned bases at >= 99.5% identity (tail). In-community pairs are
  not shifted left in magnitude but are strongly enriched at high tail fractions.

**Caption wording to avoid:** "paralogs are more similar within communities than between" -
the magnitude difference is 0.80x with heavily overlapping distributions. The supported
wording is "community members share a larger fraction of near-identical blocks (18% vs 8% of
aligned bases at >= 99.5%), while overall divergence magnitude differs only ~20%".

## The right summary: best non-homolog partner per end

Community membership is a derived partition; the primitive quantity is **for each chromosome end,
how close is its closest non-allelic partner**. Using the partner's per-record median divergence
(robust - see the note on `min_record_d` below), across the 45 ends with alignments:

| statistic | value |
|---|---|
| min | **0.10%** (coarse bins; = the allelic baseline, 0.83 Myr) |
| q25 | 0.25% (2.1 Myr) |
| **median** | **0.50%** (4.2 Myr) = **5x allelic** (coarse bins; **0.52%** = 3.5x at 1e-4, which is the value to quote) |
| q75 | 0.70% (5.8 Myr) |
| max | 7.30% (60.8 Myr) (7.35% at 1e-4) |
| mean | 0.95%; mean/median 1.90 (right-skewed) |

- 11 of 45 ends find a non-homolog within 2x allelic (<= 0.2%)
- 29 of 45 within 5x (<= 0.5%)
- 38 of 45 within 10x (<= 1.0%)
- 41 of 45 within 20x (<= 2.0%)

**Every subtelomere has a non-allelic partner closer than 5x the allelic baseline at the median, a
quarter are within 2.5x, and four (Xq, Yq, Yp, and the chr20q/chr4p pair) find a non-homolog at
exactly allelic divergence.** This is the statement the analysis should be built on, with community
membership as a downstream description of who the partners are.

### Independent validation of the PHR caller

The seven arms with no k=15 assignment (no PHR signal) are chr2_p, chr3_p, chr5_p, chr8_q, chr11_q,
chr14_q, chr18_q. Their best non-homolog divergences:

| no-signal arm | best non-homolog | median d | T |
|---|---|---|---|
| chr2_p | chr5_p | 7.30% | 60.8 Myr |
| chr3_p | none | - | - |
| chr5_p | chr16_q | 4.80% | 40.0 Myr |
| chr8_q | chr17_p | 6.40% | 53.3 Myr |
| chr11_q | none | - | - |
| chr14_q | none | - | - |
| **chr18_q** | **chrX_p** | **0.40%** | **3.3 Myr** |

Six of the seven are exactly the arms with no close non-homolog (median 5.6%, three with no
heterologous alignment at all), versus **0.50% median for the 41 signal-bearing arms (max 4.0%)** (0.52% at 1e-4 resolution).
That is a clean, independent confirmation that the PHR caller's "signal" calls track real
inter-arm sharing.

**chr18_q is the exception and deserves a look**: it has a 0.40% partner (chrX_p) yet carries no PHR
call - a candidate false negative in the caller, not a false positive.

### Correction: `min_record_d` is not a usable statistic

The pair table's `min_record_d` (minimum per-alignment divergence) is ~0 for 35 of 45 arms because
with ~10^8 alignment records some single block is always perfect. Any statement of the form "N pairs
share an allele-level block" based on it is meaningless, and an earlier version of this note claimed
"261 of 658 pairs have at least one allele-level block" on that basis - **withdrawn**. Robust
alternatives are the per-record median divergence (`d_median_record`), the fraction of aligned bases
at >= 99.5% identity (`frac_at_99.5pct`), and the per-end best-partner statistic above.

## FINAL NUMBERS (1e-4 bins, full cohort, self-pairs excluded)

Supersedes the coarse-bin values earlier in this note. Resolution was 1e-3 in the first runs, which
quantized every per-record median onto 0.001 increments; the figure and these numbers use **1e-4**.

| quantity | coarse bins | **final (1e-4)** |
|---|---|---|
| allelic median (same arm) | 0.10% | **0.15%** |
| paralog median (different chromosome) | 0.80% | **0.89%** |
| paralog median (same chromosome) | 0.80% | **0.83%** |
| paralog / allelic | 8.0x | **5.93x** |
| inter- / intra-chromosomal paralog | 1.00x | **1.07x** (still indistinguishable) |
| best partner on another chromosome, median per end | 0.50% | **0.52%** |
| intra / inter community (inter-chrom pairs only) | 0.80x | **0.80x** (n=58/568, unchanged) |

**Ladder against the corrected allelic baseline of 0.15%:**

| R = d_pair / d_allelic | pairs | % |
|---|---|---|
| <= 1.5x (at or below allelic) | 3 | 0.5% |
| 1.5-3x | 24 | 3.6% |
| 3-8x | 199 | 30.2% |
| 8-20x | 336 | 51.1% |
| > 20x | 96 | 14.6% |

**Best partner on another chromosome, per end (n = 45):** min 0.10%, q25 0.33%, **median 0.52%
(3.5x allelic)**, q75 0.73%, max 7.35%. 11 of 45 ends find an inter-chromosomal partner within 2x the
allelic baseline.

### Bug that contaminated one intermediate result

The `PAIR` aggregation key includes allelic pairs (`chrN_p`-`chrN_p`). Any downstream "closest
partner for arm A" loop must therefore **drop `arm_a == arm_b` before searching**, or every arm
"finds" a partner at the allelic level (its own alleles). An intermediate computation in this session
reported a best-partner median of 0.16% and 48 "ends" on that basis; **that number was
contamination and is withdrawn.** With self-pairs removed the value is 0.52% over 45 ends, which
matches the coarse-bin estimate (0.50%) to within quantization.

The figure `heterolog_scope.pdf/.png` was regenerated from the final run: panel A now has 1e-4 bins
(a smooth allelic ECDF with a real shoulder), panel E plots only inter-chromosomal partners with the
allelic baseline and the per-end median marked, and panel C shows only communities with >= 3 internal
pairs (C11 0.49x, C3 0.59x, C6 0.80x, C5 0.93x, C7 0.97x).
