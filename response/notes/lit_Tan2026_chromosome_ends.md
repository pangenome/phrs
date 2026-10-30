# Companion paper: Tan et al., "Widespread structural variations at human chromosome ends"

**Tan K-T, Ong RJX, Kee BBR, Yap RKH, Ng AJT, Rajapaksha GK, Garrity-Janger M, Lin Q, Chong JR,
Kyi CT, Human Pangenome Reference Consortium, Meyerson M, Li H.** bioRxiv 2026.
doi:10.64898/2026.09.03.748401 (posted 2026-09-06). Europe PMC PPR1313947. Zenodo data:
10.5281/zenodo.21807156 ("Human Subtelomere Sequences"). CC BY-NC.

Affiliation: NUS / NUS Medicine (Singapore) with Dana-Farber + Harvard + Broad. Corresponding: Kar-Tong
Tan and Heng Li — the "formerly at Dana-Farber" author is Heng Li; Tan is now independent at NUS and was
a member of the Li lab.

**Status of this note:** built from the Europe PMC abstract and the journal landing page. The full text
is behind an anti-bot block for our fetchers (bioRxiv returns 429/406 to this host); the numbers quoted
in §2 therefore come from the abstract only, and any methods-level rate assumption on their side is
**not yet verified**. Obtain the PDF (drop it in `~/Downloads/papers/`) before quoting their dating
approach in the manuscript.

## 1. What they did

Complete subtelomeric sequences for **all non-acrocentric chromosome arms** from **860 haploid
assemblies** across **six ancestry groups**, i.e. roughly twice our 465 assemblies and with ancestry
labels we do not analyse.

## 2. What they report (abstract-level)

1. Subtelomeres are **mosaics of blocks shared between chromosome arms** — the same architecture we
   describe (764 arm pairs, 98.6% same telomere-relative orientation, 9.6 Gbp).
2. Individual arms diverge extensively: **54% (21 of 39) of non-acrocentric autosomal arms carry
   multiple haplotypes differing by up to 100–200 kb**.
3. Blocks are **broadly conserved across the great apes**.
4. In humans their diversity is **associated with chromosome arm rather than ancestry**, which they
   read as evidence that "cross-arm paralogy block duplications **predate human population
   divergence**", with some ancestry-specific enrichment of haplotypes.
5. Divergent haplotypes **differ in gene content** — gene copy-number variation between individuals in
   olfactory receptors and other genes.
6. They report **"rare subtelomeric recombination"**.
7. Their sequence set improves **telomere length measurement** per chromosome end from long reads.

## 3. Where we agree

| their observation | our corresponding result |
|---|---|
| mosaics of blocks shared between arms | 764 arm pairs, 98.6% same telomere-relative orientation, 9.6 Gbp compared; block repertoire Jaccard 0.83 within community vs 0.57 between |
| 21/39 arms carry multiple haplotypes differing by 100–200 kb | 15,668 PHRs across 18,827 flanks; only 41 of 48 arms carry signal; haplotype-specific presence is the norm (9q:h2 carries the chr3-like module, 9q:h1 does not) |
| blocks conserved across great apes | our block-level d_sub ≈ 0.6–1.5% in humans (5–12 Myr upper bound, 2–5 Myr likely) — the same order as great-ape conservation age; and within-community divergence ≈ between-community divergence, i.e. the *blocks* are old |
| arm-associated rather than ancestry-associated | consistent with the standing-polymorphism reading of chr9q that we reached from the pedigree (module in both parents, in the grandmother, absent from the grandparents' other chr9 copies) |
| gene copy-number variation in ORs | our six-copy module described below; enrichment analysis: 39 families hold 187 copies, 8 families hold 142, 124/209 terms single-family |
| non-acrocentric coverage | we exclude acrocentric p arms plus X arms from the 3.51 Mb total — same scoping decision |
| "rare subtelomeric recombination" | our rate: 2.3% of arms carry a strong candidate, 34 candidates from 36 arms, 7 validated by transmission — rare by the same standard, but **recurrent across the population** rather than singular |

## 4. Where we differ, and where each of us may be wrong

1. **Their unit is the arm's own haplotype diversity; ours is the cross-arm community.** They count
   *how many distinct haplotypes one arm carries*; we count *which arms share which blocks*. These are
   two projections of one structure, and their 100–200 kb haplotype classes should correspond to
   combinations of our communities. Falsifiable and worth doing: map their haplotype classes onto our
   arm-level communities.
2. **"Rare" may be under-called.** A per-arm haplotype census cannot see a block that entered an arm
   recently and is still at low frequency; a population-wide multi-alignment can. Our 9q/3q case is a
   54 kb module present on one chromosome in one family — invisible as a "haplotype" of that arm in a
   per-arm clustering, but a clear exchange candidate in our framework. Their rare is probably our
   rare-plus-recent.
3. **Dating.** Their inference is qualitative ("predate human population divergence", "conserved
   across great apes"). Our molecular-clock estimate
   (`community_coalescence_times.md`) puts block coalescence at **5–12 Myr upper bound, 2–5 Myr if the
   subtelomeric/paternal mutation rate is 3×10⁻⁸**, which sharpens their statement into a number and
   also bounds their "ancestry-specific enrichment" as a much younger overlay. If their full text does
   contain a rate-based estimate we should compare the rate assumption directly — that is exactly the
   kind of method borrowing worth doing before submission.
4. **Copy number.** They report haplotype-level 100–200 kb differences; we find copy-number ambiguity
   at a *single* locus is itself unresolvable from long reads (flagger: `Dup` on the carrier, `Col`
   everywhere else; ONT collapses terminal arrays in every sample). Both papers are measuring the same
   homogenization; neither can assign copy number per end without a copy-aware method. Any critique of
   our module architecture on copy-number grounds applies to their haplotype counts as well.

## 5. What to do next

1. **Get the full text** (bioRxiv blocks our fetchers). Check: their mutation rate / dating method; their
   definition of "rare"; whether haplotype resolution is sequence- or graph-based; whether they cite
   and criticise us explicitly, and on what point.
2. **Map their haplotype classes to our communities** (§4.1) — a real cross-validation opportunity and a
   natural "we confirm and extend" result.
3. **Cite them in the manuscript** as the companion resource: independent support for cross-arm block
   mosaicism and for OR gene copy-number variation, with our population-wide exchange framework and
   molecular-clock bound as the extension.
4. Consider a short **response/companion paragraph** in Discussion rather than Methods, and only after
   the full text is read.
