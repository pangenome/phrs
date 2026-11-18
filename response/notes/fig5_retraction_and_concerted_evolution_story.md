# The Fig. 5 retraction, and the story that replaces it

## What the preprint claimed

"Ectopic subtelomeric recombination within a pedigree": the daughter's paternal chr9q terminus was
reported to be a chr3q-derived ectopic recombination product, with the drop in match to the father's own
chr9q as the evidence, plus chr1p/chr5q as a second candidate described as ancestral homology.

## What complete evaluation shows

**1. The chr9q call was an assembly artefact.** The terminal ~150 kb of the father's chr9q haplotype 2 is
misassembled in both releases; release QC marks a misjoin at the start of the block.
**2. The chr1p/chr5q call was an alignment artefact.** At the scan's seed-frequency cutoff the father's
chr5q candidate never enters the alignment, so chr1p wins by default.

**3. The competitive test, on the metric the paper uses.** All candidates offered simultaneously to IMPG
similarity (sweepga `--num-mappings many:many`, `--fastga-frequency 5000`; 2 kb bins; `impg similarity
--fast`), including the corrected and uncorrected paternal chr9q, his other chr9q haplotype, chr3q,
chr19q, chr1p and a no-signal control:

| candidate | tract dice | tract est. identity | outside dice |
|---|---|---|---|
| SELF.p (self) | 1.0000 | 1.0000 | 1.00 |
| **DAD chr9 h2, corrected** | **1.0000** | **1.0000** | **1.00** |
| DAD chr3 h2 | 0.9915 | 0.9994 | 0.98 |
| DAD chr9 h2, original (erroneous) | 0.9211 | 0.9945 | 1.00 |
| DAD chr9 h1 | 0.9191 | 0.9944 | 0.87 |
| DAD chr1 h1 | 0.9144 | 0.9941 | 0.94 |
| DAD chr19 h2 | no coverage | - | 0.05 |

**Corrected, the homolog wins outright at dice 1.0000**; chr3q is second at 0.9915. The original
assembly scored 0.9211 against the daughter's own copy, which is what let chr3q win.

## The story that replaces it

**A non-homologous end lies inside the spread of the homologous ones.** The father's chr3q terminus is a
*closer* match to his daughter's chr9q (dice 0.9915) than his own other chr9q haplotype is (0.9191).
That is the concerted-evolution signal itself, and it is now measured on the same metric as everything
else in the paper: chr3q sits within the margin of difference between the two chr9q's.

**And no exchange is detectable in the pedigree outside the obligate compartments.** The only
heterologous sharing that behaves like obligate exchange in a three-generation family is
acrocentric-acrocentric and Xp/Yp: those are the only regions where non-homologous copies are as similar
as alleles (ratio 1.00 and 0.99). Everywhere else the sharing is family-level similarity between blocks
that have been exchanging for far longer than a pedigree can observe.

## Metric correction made along the way

`estimated.identity` in the IMPG `--fast` output is **not** dice: it is the Mash transform
`(2J/(1+J))^(1/k)` of the Jaccard, while dice is `2*shared/(A+B)`; the two coincide only when sketch
sizes are equal. This is why the figure's identity axis compressed the chr3q-vs-homolog difference
(0.9994 vs 1.0000) while dice separates it clearly (0.9915 vs 1.0000). **Fig. 5B should plot dice.**

## Applied to the manuscript (build exit 0, 0 undefined, commit in `BUILD_LOG.md`)

- Abstract: the pedigree sentence now states that no exchange is found outside acrocentrics and PARs, and
  that the one candidate is an assembly artefact.
- Results: `Ectopic subtelomeric recombination within a pedigree` replaced by `Direct test for exchange
  between non-homologous ends`, which reports the artefact, the corrected competitive test, and the
  concerted-evolution reading.
- Fig. 5 caption retitled `An assembly error in one paternal haplotype mimics an ectopic exchange`,
  with the corrected-assembly numbers and PAR1 as the only exchange; the legend swatch now reads
  `9q/3q artefact`.

## What this strengthens

The paper's claims no longer depend on a single pedigree event. What remains is stronger and broader:
15,668 PHRs across 465 assemblies, communities with a measured three-tier structure, an exchange rate of
~1 per 740,000 generations per locus pair from ~50 kb tracts, a coalescence bound of 5-12 Myr, and the
TERRA-promoter junction. The pedigree now contributes a *negative* control - a demonstration that calls
of this kind fail - which is a more useful contribution than a single unreplicated event.
