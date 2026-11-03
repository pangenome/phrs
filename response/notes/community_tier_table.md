# Community tiers: allelic, near-allelic, and residual

Derived from the cohort run (`hprc_het_final.json`): per community, the member arms' own allelic
divergence, their intra-community paralogous divergence, and their divergence to out-of-community
partners. Inter-chromosomal pairs only; per-record medians.

| community | arms | allelic | intra-paralog | ratio | out-paralog | out-ratio | verdict |
|---|---|---|---|---|---|---|---|
| **C7** | 13p, 14p, 15p, 21p, 22p | **1.130%** | **1.135%** | **1.00** | 1.085% | 0.96 | allelic-like (acrocentric rDNA/satellite arrays) |
| **C15** | Xp, Yp | 0.455% | 0.450% | **0.99** | 1.605% | 3.53 | allelic-like (PAR1, obligate X-Y recombination) |
| **C14** | Xq, Yq | 0.070% | 0.100% | **1.43** | 1.415% | 20.21 | near-allelic (PAR2) |
| C1 | 4q, 10q | 0.190% | 0.620% | 3.26 | 1.990% | 10.47 | near-allelic, cohesive |
| **C11** | 1p, 5q, 6q, 8p | 0.170% | 0.585% | **3.44** | 1.330% | 7.82 | near-allelic, cohesive |
| **C3** | 3q, 7p, 9q, 11p, 16q, 19p | 0.190% | 0.670% | **3.53** | 1.290% | 6.79 | near-allelic, cohesive |
| C12 | 20p, 2q | 0.185% | 1.340% | 7.24 | 1.400% | 7.57 | weakly cohesive |
| C5 | 6p, 9p, 12p, 20q | 0.140% | 1.365% | 9.75 | 1.310% | 9.36 | weakly cohesive |
| C2 | 10p, 18p | 0.105% | 1.280% | 12.19 | 2.800% | 26.67 | weakly cohesive |
| C6 | 1q, 13q, 17q, 19q, 21q, 22q | 0.155% | 1.930% | 12.45 | 2.240% | 14.45 | weakly cohesive |
| C4 | 7q, 12q | 0.150% | 4.070% | **27.13** | 2.390% | 15.93 | **residual, not a community** |

## The three tiers

1. **Allelic-like (ratio ~1): C7, C15, C14.** In the acrocentric p arms, alleles and non-homologs are
   *equally* divergent (1.13% vs 1.14%, ratio 1.00; 0.96 against outside partners). Homologs and
   paralogs are indistinguishable - complete homogenisation of a single repeat pool, with known
   mechanisms in all three cases (rDNA/satellite arrays; PAR1 obligate X-Y crossover; PAR2).
2. **Near-allelic and cohesive (ratio 3.3-3.5): C1, C11, C3.** Members sit 3.5x their allelic baseline
   but 7-10x against outside partners - genuine shared-block communities. **The two pedigree exchange
   candidates fall in this tier**: chr3q/chr9q in C3, chr1p/chr5q in C11.
3. **Weak or residual: C5, C2, C6, C12, and C4.**

## C4 is not a community

chr7_q-chr12_q is the only pair in C4: **3.775% weighted, 4.0% median, fraction of aligned bases at
>= 99.5% identity = 0.001**. There is essentially no near-identical sequence shared.

- chr7_q's best partners are chr17_q, chr16_q, chr1_p and chr2_q at ~2.2% - **chr12_q is not among
  chr7_q's best partners at all**.
- chr12_q's best partner is chr7_q at 3.78%, with everything else at 6.8-7.7% - it has no good partner
  anywhere.

So C4 is a residual pairing of two weakly connected arms, produced because k=15 forces every arm into
some community. It should be reported as such (or merged), not presented as a community with internal
cohesion.

## Why this table matters for the manuscript

It replaces a single global "paralog/allelic = 5.9x" with a structured result: the subtelomere
contains at least three qualitatively different regimes - fully homogenised arrays (ratio 1), cohesive
shared-block communities (ratio ~3.5), and residual one-off pairings (ratio >7, no identical blocks) -
and the pedigree exchange candidates arise specifically in the cohesive middle tier.
