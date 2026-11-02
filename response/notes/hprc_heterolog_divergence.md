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
