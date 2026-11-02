# Is the telomere-proximal homogenisation a property of the genes there?

Script `scripts/ci/telomere_gene_profile.py`. Full cohort, 18,827 window PAFs, 2 kb bins, every bin
split by whether the alignment midpoint falls inside a RefSeq gene (CHM13 Liftoff v5.2, projected to
distance-from-telomere per arm). Genes are 36.8% of the telomere-proximal 20 kb and are the canonical
subtelomeric families: **WASH** (WASH8P chr12p 2.5-20 kb, WASH3P chr15q 2-20 kb), **DDX11L**
(L5/L8/L17), **DUX4**, **IGH** (chr14q), **PGBD2**, **RPL23AP** paralogs, **LOC107986338**.

## Answer: no. The tip effect is INTERGENIC, and the genes are not homogenised.

| distance (kb) | gene allelic | gene paralog | **gene ratio** | intergenic allelic | intergenic paralog | **intergenic ratio** |
|---|---|---|---|---|---|---|
| **0-2** | 0.300% | 2.960% | **9.87** | **2.530%** | 2.940% | **1.16** |
| 2-4 | 1.170% | 1.850% | 1.58 | 0.520% | 2.040% | 3.92 |
| 4-6 | 1.450% | 1.560% | 1.08 | 0.640% | 1.680% | 2.62 |
| 6-8 | 0.310% | 1.410% | 4.55 | 0.240% | 1.200% | 5.00 |
| 8-10 | 0.780% | 1.250% | 1.60 | 0.430% | 1.100% | 2.56 |
| 10-12 | 0.250% | 1.170% | 4.68 | 0.200% | 1.420% | 7.10 |
| 12-16 | 0.27-0.32% | 0.97-1.56% | 3.6-4.9 | 0.14-0.22% | 1.07-1.41% | 6.4-7.6 |
| 20-40 | 0.09-0.14% | 0.51-1.08% | 3.3-9.8 | 0.18-0.21% | 0.52-0.80% | 2.5-4.4 |
| 50-200 | 0.11-0.13% | 0.47-0.69% | 4.3-6.3 | 0.11-0.17% | 0.59-0.69% | 3.5-5.5 |

## Three findings

**1. The unmasked tip signal is entirely intergenic.** The previously reported "ratio 1.30 at 0-2 kb"
decomposes into **intergenic ratio 1.16** and **genic ratio 9.87**. The terminal 2 kb is where the
TTAGGG-adjacent repeat arrays live, and those arrays have allelic divergence of **2.53%** - 17x the
interior allelic level - while their paralog divergence (2.94%) is barely higher. So the two classes
converge there because *both* are inflated by repeat-array alignment ambiguity, not because exchange
homogenises them. This is the artefact explanation, and it is now the supported one.

**2. Genes at the very tip are the LEAST homogenised sequence in the subtelomere.** At 0-2 kb the
genic ratio is 9.87, driven by normal allelic divergence (0.30%) against strongly divergent paralogs
(2.96%) - higher paralog divergence than the tip intergenic value. There is **no** evidence of
concerted evolution operating on WASH, DDX11L, IGH or PGBD2 at the telomere junction.

**3. The low genic ratios at 2-6 kb (1.08-1.58) are an orthology artefact, not exchange.** In those
bins *allelic* divergence is 1.17-1.45%, roughly 10x the interior allelic level, which means same-arm
alignments there are not aligning orthologous copies: these are multicopy pseudogene families
(WASH arrays, RPL23AP) where "allelic" and "paralogous" comparisons sample the same copy set. A ratio
near 1 produced by inflating the denominator is not evidence of homogenisation.

## Consequence

The claim "subtelomeric homogenisation is strongest at the telomere, and the subtelomeric genes are
under concerted evolution" is **not supported**. What the data support is narrower and cleaner:
**allelic and paralogous divergence both rise steeply in the terminal few kb, converge there, and
that convergence is a property of the terminal repeat arrays, not of the genes.** The gene-level
behaviour (paralog divergence ~1-3% against allelic ~0.1-0.3%, ratio ~5-10) is the ordinary
subtelomeric-paralog regime measured elsewhere in this analysis.

Whether WASH/DDX11L specifically are homogenised between arms remains open, but it cannot be answered
with arm-level allelic comparisons in these assemblies - it needs copy-resolved (per-locus) comparison
of gene copies, which the multi-mapping assemblies do not provide.
