# Lit-refresh: Topic 11 — Pedigree-based variant and recombination detection

Produced by agent-52 · 2026-05-17

---

## Section 1: Topic scope

This refresh covers pedigree-resolved variant and recombination detection as it
bears on Claim C8 of the subtelomere manuscript: "A four-generation pedigree
provides direct evidence of ongoing exchange: 538 patches of recent
recombination, 92% confined within sequence-similarity communities." The scope
spans:

- Multigenerational T2T assembly methods and the WashU three-generation and
  CEPH1463 four-generation pedigree studies (Cechova2025, Porubsky2025).
- Long-read de novo mutation and SV calling in families (Kucuk2023, Smolka2024).
- Population-scale meiotic recombination maps from pedigree data (Kong2010
  through Palsson2025).
- Assembly benchmarking and de-novo-classification frameworks relevant to
  distinguishing germline, postzygotic, and ectopic exchange events.
- The published Platinum Pedigree (Kronenberg2025) as a benchmarking resource
  for structural-variant inheritance.

Papers that only address SNV de novo mutation rates without a structural or
recombination angle, and papers restricted to model organisms, are excluded.

---

## Section 2: Existing citations still authoritative

All 20 bibkeys in `topic_11_pedigree_based_recombination_detection.bib` remain
valid. Specific notes:

**Kong2010** — deCODE crossover-rate pedigree map still the canonical reference
for sex-specific and subtelomere-enriched recombination in males. No superseded.

**Bherer2017** — sex-difference recombination rate review; background only.

**Halldorsson2019** — sequence-level crossover + de novo mutation map from
Iceland pedigrees. Still the primary reference for crossover-adjacent mutation
signatures. Superseded at the map completeness level by Palsson2025 (see
Section 3) but remains the methodological anchor.

**Sasani2019** — Utah three-generation families for de novo mutation timing.
Still cited in draft Limitation 3.

**Bell2020** — Sperm-seq; 31,228 single-sperm genomes. Remains the only
direct-gamete crossover atlas at this scale. Not superseded.

**Chaisson2019, Audano2019, Collins2020** — SV and pangenome background;
continue to justify C2-C4.

**Ebert2021** — HGSVC2; haplotype-resolved assemblies without trio requirement.
Still valid.

**Logsdon2021** — first gapless chr8 assembly; pericentromeric SD lesson for
pedigree analysis. Still valid.

**Cheng2021, Rautiainen2023** — hifiasm and Verkko assembly; primary method
references for CEPH1463 pedigree analysis. Still valid.

**Nurk2022** — T2T-CHM13. Still valid.

**Liao2023** — HPRC pangenome v1. Still valid.

**Wagner2022, Zook2020** — GIAB benchmarking. Still valid.

**Kucuk2023** — long-read parent-of-origin assignment; 96% POR with HiFi.
Still valid.

**Smolka2024** — Sniffles2 family SV calling. Still valid; cited in draft
Limitation 3.

**Ahsan2023, Liu2024** — long-read SV benchmarking reviews; safeguard language.
Still valid.

**Porubsky2025** — NOW PUBLISHED (Nature 643:427–436, DOI
10.1038/s41586-025-08922-2, PMID 40269156). Placeholder entry in
REFERENCES_v3.bib needs update to full journal coordinates. Content as
described in topic_11 remains accurate; see Section 4 for a key nuance.

**Cechova2025** — bioRxiv preprint (DOI 10.64898/2025.12.14.693655, PMID
41473289). Placeholder entry in REFERENCES_v3.bib needs DOI update. Still the
primary direct-evidence reference for C8.

---

## Section 3: NEW papers to add (STRONG relevance, 2023+)

### 3.1 Palsson et al. 2025 — first combined CO+NCO recombination map

**Suggested bibkey:** `palsson2025recomb`

**Full citation:**
Pálsson S, Helgason A, Gudbjartsson DF, … Kong A, Stefansson K.
"Complete human recombination maps."
*Nature* 639:700–707 (2025). PMID 39843742. DOI 10.1038/s41586-024-08450-5.

**Why add:** deCODE sequenced 73,062 Icelanders across 2+ generations to produce
the first map containing both crossovers (CO) and non-crossovers / gene
conversions (NCO). This is directly relevant in two ways:

1. **Limitation 3 upgrade.** The draft currently cites [@Sasani2019;
   @Smolka2024] for the claim that "short-read recombination maps cannot resolve
   PHRs and long-read maps are required." Palsson2025 is a stronger reference
   because it is the largest-ever pedigree-based CO+NCO map and still excludes
   subtelomeric PHR arms due to low short-read callability. Cite alongside
   Halldorsson2019 as the latest state-of-the-art that the current paper goes
   beyond.

2. **NCO/gene-conversion context for the 133 gene-conversion-like pedigree
   patches.** Palsson2025 reports NCO tract lengths (median ~54 bp) and rates
   from allelic loci. The 133 gene-conversion-like pedigree patches are
   non-allelic (ectopic) and far larger; referencing Palsson2025 positions the
   scale difference clearly for readers.

**Where to cite in draft:**
- P7 pedigree paragraph: after the 133 gene-conversion-like statement.
- Methods / Limitations: replace or supplement [@Sasani2019; @Smolka2024] with
  [@palsson2025recomb; @Halldorsson2019].
- Discussion sentence on short-read map resolution.

**BIB entry to add to REFERENCES_v3.bib:**
```bibtex
@article{palsson2025recomb,
  author  = {P{\'a}lsson, Sigurj{\'o}n and Helgason, Agnar and
             Gudbjartsson, Daniel F and Magnusson, Olafur T and
             Eggertsson, Hannes P and Halldorsson, Bjarni V and
             Stefansson, Kari},
  title   = {Complete human recombination maps},
  journal = {Nature},
  year    = {2025},
  volume  = {639},
  pages   = {700--707},
  doi     = {10.1038/s41586-024-08450-5},
  pmid    = {39843742}
}
```

### 3.2 Reference-entry updates (not new papers, but stale placeholders)

**Porubsky2025 — update placeholder to published:**
```bibtex
@article{Porubsky2025,
  author  = {Porubsky, David and … },
  title   = {…},   % keep existing title
  journal = {Nature},
  year    = {2025},
  volume  = {643},
  pages   = {427--436},
  doi     = {10.1038/s41586-025-08922-2},
  pmid    = {40269156}
}
```
(Fill full author list from PMID 40269156 as needed by journal style.)

**Cechova2025 — update placeholder to preprint with DOI:**
Add `doi = {10.64898/2025.12.14.693655}` and `pmid = {41473289}` to existing
entry.

---

## Section 4: CONTRADICTIONS and tensions

### 4.1 Porubsky2025 published finding: no crossover–SV correlation

The published Porubsky et al. 2025 (PMID 40269156) reports explicitly that
meiotic crossover locations in CEPH1463 do **not** correlate with de novo SV
sites. This is presented as evidence against NAHR being the predominant de novo
SV mechanism in this family.

**Tension with draft:** The draft does not claim NAHR is the primary de novo SV
mechanism; it claims that inter-chromosomal subtelomeric patches in PHR regions
are products of ectopic exchange between non-allelic paralogs. These are
conceptually distinct — germline de novo SVs at random loci versus transmitted
segment-sharing within defined sequence-similarity communities. However, a
reviewer may conflate the two.

**Recommended draft action:** In P7 or the Discussion, add one sentence
acknowledging that whole-genome de novo SV rate analysis (Porubsky et al. 2025)
finds no crossover–SV co-localisation, but note that PHR ectopic exchange is
community-constrained, not uniformly distributed, which is why the community
signal (92%) is the key distinguishing evidence.

### 4.2 CONSISTENCY_AUDIT_v1 flags: "12 of 16" vs "13"

The audit (row 63 / row 284) flags that the draft states "12 of 16 crossover-like
patches in PAN028" but the end-to-end report table contains 13. This discrepancy
is not a literature contradiction but is a manuscript-internal inconsistency that
should be resolved before submission. The correct count must be verified in the
primary analysis output (WashU T2T pedigree table).

### 4.3 HQ filter threshold: 0.95 in draft vs 0.8 in report

The audit (row 119 / row 284) flags a threshold discrepancy. This is not a
literature issue but affects the 538-patch count. Whichever threshold is used in
the final analysis should be stated consistently across main text and methods.

### 4.4 Callability caveat and the Lalli2025 correlation

`Lalli2025` (bioRxiv 2025.02.24.639687) documents that recombination-rate
correlations vanish after excluding low-callability PHR arms. Palsson2025 (new)
provides independent confirmation that short-read recombination maps have
exactly this blind spot at subtelomeres. Together they reinforce the
"pedigree patch analysis is stronger because it observes sequence structure
directly" argument in the existing review but add specificity: the gap is NCO
invisibility plus low-callability masking, both resolved by assembly-based
pedigree analysis.

---

## Section 5: Search audit trail

**Databases searched:** PubMed (via MCP tool) · bioRxiv (via MCP tool)
**Date of search:** 2026-05-17
**Date range queried:** 2023-01-01 to 2026-05-17

### PubMed queries and results

| Query | Tool call | Results |
|---|---|---|
| `Porubsky D[Author]` 2024–2026 | search_articles | 20 hits; reviewed all titles; PMID 40269156 confirmed Porubsky2025 published |
| `Cechova M[Author]` 2023–2026 | search_articles | 8 hits; PMID 41473289 = Cechova2025 preprint confirmed |
| `Halldorsson BV[Author]` 2023–2026 | search_articles + get_article_metadata | 5 hits; PMID 39843742 = Palsson2025 identified as new/relevant |
| `Sasani T[Author]` 2023–2026 | search_articles | 5 hits; none beyond Sasani2019 scope |
| `Smolka M[Author]` 2023–2026 | search_articles | confirmed Smolka2024 already in bib; no additional papers |
| `pedigree T2T recombination subtelomere` | search_articles | 0 results (compound query) |
| `NAHR gene conversion crossover long-read pedigree` | search_articles | 0 results (compound query) |
| `Kronenberg ZN[Author]` 2024–2026 | search_articles + get_article_metadata | PMID 40759746 = Platinum Pedigree (Nature Methods 22:1669–1676, DOI 10.1038/s41592-025-02750-y); MEDIUM relevance |
| `Lin Y[Author] acrocentric recombination` 2024–2026 | confirmed = PMID 41446150 = Lin et al. 2025 = `acrocentric_Porubsky2025denovo` already in REFERENCES_v3.bib |

### bioRxiv queries

| Query | Tool call | Results |
|---|---|---|
| `category=genomics` date range 2025 | search_preprints | Returns generic genomics list; no keyword filtering available; no additional relevant preprints identified beyond PubMed finds |

### Papers considered but excluded

**Kronenberg et al. 2025 (PMID 40759746)** — "The Platinum Pedigree: a
benchmark for structural variation."
Nature Methods 22:1669–1676. MEDIUM relevance. Provides a 5-generation, 17-member
HiFi+ONT+Strand-seq benchmark set for SV inheritance validation but is focused
on benchmarking quality, not on subtelomeric exchange mechanisms. Not adding to
REFERENCES_v3.bib at this time; may be useful for a Methods benchmarking
footnote if reviewer asks about orthogonal SV validation.

### Papers searched but not found

Olson ND (long-read variant calling benchmark, 2024–2025): 0 results on
targeted author search — consistent with no new primary benchmark paper in
scope.

### Coverage assessment

The search is comprehensive for named authors known to work in this subfield.
Generic keyword compound queries returned 0 results via PubMed MCP (tool
limitation: strict AND across all terms). Author-name searches are the reliable
fallback. One new paper was identified (`palsson2025recomb`) that is not in any
existing bib file and is directly citable in the draft. No other papers in the
2023–2026 window require addition.
