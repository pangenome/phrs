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

## Correction: there is no mask bug - the allelic denominator itself fails at the tip

An earlier revision of this note claimed the genic 2-6 kb values were "internally inconsistent with the
per-arm allelic data, therefore a bug". **That reasoning was wrong and is withdrawn.** The per-arm
allelic median is an average over the whole 500 kb window, dominated by the interior, whereas a
distance bin is a positional subset. Comparing the two is not a valid consistency check.

Verified directly instead, on a single PAF whose target is gene-rich
(`all-vs-HG002_1_chr15_95754595_96254594_chr15_qarm`), recomputing every class/bin by hand:

| distance | class | genic | n | median divergence |
|---|---|---|---|---|
| 0-2 kb | paralog | no | 81 | 1.42% |
| 2-4 kb | paralog | yes | 497 | 2.10% |
| **4-6 kb** | **allelic** | **yes** | **357** | **1.43%** |
| 4-6 kb | paralog | yes | 2396 | 1.66% |
| 8-10 kb | paralog | yes | 465 | 1.33% |
| 100-150 kb | allelic | no | 580 | 0.16% |
| 150-200 kb | allelic | no | 370 | 0.12% |
| 200-500 kb | allelic | yes | 1159 | 0.10% |

The bin reproduces the aggregate exactly. So the mask is doing what it was written to do.

**The real problem is the denominator.** Allelic divergence within ~10 kb of the telomere is
0.30-2.53% (genic 0.30-1.45%, intergenic 0.52-2.53%), against 0.09-0.16% in the interior. Human
intraspecific divergence is ~0.1-0.3% even in hypervariable regions, so **an "allelic" value of
0.6-2.5% cannot be orthologous divergence** - same-arm alignments in the terminal few kb are not
aligning homologous positions. This applies to genes and intergenic sequence alike.

**Consequence: the telomere-distance profile is uninterpretable inside roughly the first 10 kb**, and
the "ratio approaches 1 at the junction" result is an artefact of that failure, not a measurement of
homogenisation. Both the genic and intergenic tip values must be withdrawn as evidence about
exchange. What remains valid from the profile is the interior (10-500 kb), where allelic divergence
is 0.09-0.16% and the paralog/allelic ratio sits at 4-8.

A defensible operational rule for any future version of this figure: **discard bins whose allelic
divergence exceeds ~0.4%** (4x the interior baseline), because the denominator has failed there.
By that rule every bin at 0-6 kb is discarded, and the first usable bin is 6-8 kb.

Testing terminal exchange specifically needs copy-resolved methods - per-locus assignment of reads or
assemblies to individual repeat units - not arm-level window comparison. That is out of reach for the
multi-mapping HPRCv2 assemblies.

## The first ~2 kb is a CpG island, and on the ends that have one it is a DDX11L promoter

Two measurements, from 30 sampled windows oriented telomere-first:

| distance | CpG observed/expected | GC% | island (O/E>0.6, GC>55%) |
|---|---|---|---|
| 0-2 kb | **1.02** | **62.3** | **YES** |
| 2-5 kb | 0.87 | 40.9 | no |
| 5-10 kb | 0.93 | 40.2 | no |

Genome-wide CpG O/E is ~0.2-0.3; values near 1.0 indicate a CpG island, i.e. a promoter. So the first
~2 kb of the subtelomere carries a promoter signature, and it is the same zone in which allelic
alignments collapse (see `telomere_distance_profile.md`).

**Annotation agrees on the ends that have a gene there.** Gene TSS distance from the telomere end in
CHM13, from RefSeq Liftoff v5.2:

- **7 of 46 ends** have a TSS within 2 kb of the tip; **14 of 46** within 5 kb.
- The recurrent gene family is **DDX11L** (DEAD/H-box helicase-like pseudogene): DDX11L1, L5 (x2), L8
  (x2), L9, L10, L16, L17 - **8 of the 15 tip genes**, with TSS 1.7-2.5 kb from the tract.
- Others at the tip: LOC102723376 (x2), LINC01986 (TSS 749 bp), LINC02684 (859 bp), IGH/IGHVIII-82
  (chr14q, 4,978 bp).

**Caveats.** (i) It is not every chromosome end - 14/46, i.e. ~30%. (ii) That count should be treated as
a LOWER bound: RefSeq Liftoff is least reliable in these terminal repeat regions, and the CpG signature
appears broader than the annotation. (iii) CpG island does not prove an expressed promoter in every
individual.

**Why this strengthens rather than weakens the exclusion of the terminal zone.** The zone is not
unreadable sequence - it is the zone where gene and promoter copy number, and their position relative to
the tract, differ between haplotypes. That is the same copy-number variation Tan et al. report for
WASH, olfactory receptors and DDX11L. Allelic alignments fragment there because the gene complement
differs between people, which is a biological finding and the reason a divergence ratio cannot be
computed - not a technical nuisance.

**Still open.** (a) Whether the GC-rich tip on the other 32 ends is an unannotated promoter or a GC-rich
repeat. (b) Whether DDX11L presence at the tip predicts lower paralog/allelic divergence on those specific
arms - the test that would connect the gene placement to the exchange regime.
