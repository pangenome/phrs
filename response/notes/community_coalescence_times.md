# Best-guess coalescence times per arm-level community

Script: `scripts/cladistics/community_coalescence.py`; results:
`scripts/cladistics/community_coalescence_communities.tsv`,
`..._pairs.tsv`, `...json`. Method is crude and stated as a bound — see §4.

## 1. Method

For two contemporaneous copies of a subtelomeric block **in one genome**, the substitution-only
divergence over aligned homologous blocks is

    d_sub = X / (= + X)          (from wfmash CIGARs; indels excluded by construction)

and, with `d_sub ≈ 2 μ T` for copies that last shared an ancestor T generations ago:

    T = d_sub / (2 μ),      μ = 1.2e-8 substitutions/bp/generation   (human pedigree average)

Reported in generations and Myr at 20 y/gen. Data: HG002 hap1, 41 telomere-anchored 500 kb windows
(`combined_with_hprcv2.1Mb.telo_500kb_trimmed.fa.gz`), aligned all-vs-all with
`wfmash -p 75 -s 1000 -w 3000 -o`; pairs required ≥ 20 kb of aligned block.

## 2. Per community

| community | arms | pairs | median d_sub | T (Myr) | T (Myr) if μ = 3e-8 |
|---|---|---|---|---|---|
| C1 | 4q, 10q | 1 | 0.00651 | **5.4** | 2.2 |
| C3 | 3q, 7p, 9q, 11p, 16q, 19p | 7 | 0.00599 | **5.0** | 2.0 |
| C5 | 9p, 12p, 20q (+6p) | 2 | 0.01472 | **12.3** | 4.9 |
| C6 | 1q, 13q, 17q, 19q, 21q, 22q | 1 | 0.01419 | **11.8** | 4.7 |
| C7 | 13p, 14p, 15p, 21p, 22p | 1 | 0.01377 | **11.5** | 4.6 |
| C11 | 1p, 5q, 6q, 8p | 5 | 0.00604 | **5.0** | 2.0 |
| C2, C4, C8, C9, C10, C12, C13, C14, C15 | — | 0 | no pair aligned ≥ 20 kb at p75 | — | — |

Median over the six communities with data: **d_sub = 0.0101 → 8.5 Myr upper bound (3.4 Myr if
μ = 3×10⁻⁸)**.

## 3. Selected pair divergences (aligned bp, d_sub, T upper bound)

| pair | aligned | d_sub | T (Myr) | note |
|---|---|---|---|---|
| chr5q–chr6q | 285 kb | 0.00498 | 4.2 | most extensive alignment of any pair |
| chr5q–chr8p | 145 kb | 0.00552 | 4.6 | same community |
| chr1p–chr5q | 91 kb | 0.00604 | 5.0 | the withdrawn candidate: copies are *least* divergent |
| chr20p–chr7p | 52 kb | 0.00670 | 5.6 | between communities |
| chr5q–chr9q | 81 kb | 0.00987 | 8.2 | between communities |
| **chr3q–chr9q** | 43 kb | **0.01228** | **10.2** | the retained candidate: arm-scale copies are the *most* divergent |
| chr1p–chr8p | 88 kb | 0.02084 | 17.4 | same community, different block |
| chr2q–chr8p | 65 kb | 0.02682 | 22.4 | between communities |

## 4. What this all means — and why every number is an upper bound

**It is not one generation and it is not the age of the chromosomes.** The copies we are dating
have diverged by 0.5–1.5% in substitutions, which is two orders of magnitude more than a single
transmission could produce and one to two orders less than the age of the duplications that created
these families. The process has been running for millions of years, and it is still running.

**Why upper bound.** Subtelomeres sit in the paternal germline and in repeat-rich sequence, where
the mutation rate is plausibly 2–4× the genome average; μ is also not constant across the block.
With μ = 3×10⁻⁸ every T shrinks by 60% (last column). The published rate is therefore the
conservative choice, and **the true coalescence times are probably 2–5 Myr, not 5–12 Myr**.

**Within-community ≈ between-community.** Median d_sub for the 17 within-community pairs is 0.00702
(5.9 Myr) and for the 23 between-community pairs 0.00765 (6.4 Myr) — indistinguishable. So the
*blocks* are shared across communities; communities differ in **which blocks they combine**, not in
the age of the blocks. This is the same finding as the module-architecture analysis (block repertoire
Jaccard 0.83 within vs 0.57 between: high everywhere) and it is what concerted evolution predicts —
old blocks, continually reshuffled.

**The mosaic at 9q/3q is visible in the numbers.** chr3q–chr9q is the *most* divergent
within-community pair at arm scale (1.23%, 10.2 Myr), while the 54 kb module inside it is the *least*
divergent sequence we have measured anywhere (0.1–0.5%, 0.8–4 Myr). That combination — a young
module in an old background — is exactly the signature of a recent transfer into a long-diverged
pair of arms, and it is what the pedigree independently showed by transmission.
By contrast chr1p–chr5q copies are 0.60% divergent at arm scale with no such young insert, i.e. an
old, uniformly shared family.

**What is dated and what is not.** What these numbers date is *the last exchange affecting each
block*, which differs block by block — the distribution across pairs, not any single value, is the
object of interest. They do not date the founding duplication of the families (which is older and
only lower-bounded, since the most diverged pairs fall below the alignment floor), and they do not
date any population process: a single genome cannot.

## 5. Caveats

1. One genome (HG002 hap1). Copies within an individual are contemporaneous, which is what the
   calculation requires, but it is one draw from the population.
2. `-p 75` alignment floor: more diverged pairs are measured on their most similar blocks only, so
   d_sub is biased **low** there.
3. Blocks shorter than 20 kb of aligned sequence are dropped, which removes many communities from
   the table (9 of 15).
4. Indels are excluded from d_sub (correct for converting with a substitution rate) but they mean
   these numbers describe substitutional change only; the indel-rich reality is far more divergent
   in length-penalised terms (d_lenpen is 3–20× larger, see the pairs table).
5. Communities are not clades (character bootstrap: TIGHT_q and 4q/10q are not monophyletic), so
   "community coalescence time" summarizes member pairs, not one ancestor.
