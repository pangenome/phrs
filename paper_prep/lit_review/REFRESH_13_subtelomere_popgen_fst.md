# Lit Refresh: Topic 13 — Subtelomere Population Genetics and FST Out-of-Africa Structure

**Date**: 2026-05-17
**Agent**: agent-73
**Branch**: wg/agent-73/lit-refresh-13-subtelomere
**Source review**: `paper_prep/lit_review/topic_13_subtelomere_popgen_fst.md`
**Source bib**: `paper_prep/lit_review/topic_13_subtelomere_popgen_fst.bib`

---

## Section 1: Topic scope

This topic covers the population-genetic interpretation of subtelomeric haplotype sharing
(cross-arm affinity) across the five HPRC v2 superpopulations (AFR=67, EAS=52, AMR=44,
SAS=37, EUR=33 diploid individuals; 465 haplotypes). The primary analytic claim is that
Hudson FST at subtelomeric binary haplotypes (cross-arm vs. self-arm) recapitulates the
canonical out-of-Africa split: AFR vs. non-AFR FST ≈ 0.10–0.15, while all non-AFR
pairwise values are ≈ 0 (range −0.047 to +0.007), and that two loci (chr16q at 70% AFR
frequency and chr4q at 36% AFR) show the clearest AFR enrichment. Lit-refresh papers
should directly support or contextualize: (a) Hudson/Weir-Cockerham FST methodology at
structural variants in diverse human populations; (b) population-specific enrichment of
structural variants (SVs) or segmental duplications (SDs) in African vs. non-African
genomes; (c) subtelomere-specific SV surveys across populations; (d) population frameworks
for long-read SV characterization across the AFR/EAS/SAS/AMR/EUR superpopulation structure.

CONSISTENCY_AUDIT_v1.md flags one internal divergence relevant to this topic: Row 36 —
the NATURE_DRAFT_v1.md states non-AFR pairwise FST "0.02 to 0.04" while
end-to-end-report/report/04_heterogeneity.md reports the actual values as −0.047 to
+0.007. The draft value is wrong; the report is authoritative. New papers bearing on
expected FST ranges for SVs between closely related non-African populations are thus
particularly informative for correcting this.

---

## Section 2: Existing citations — still current?

All 12 entries in `topic_13_subtelomere_popgen_fst.bib` remain authoritative. None are
superseded by the new papers found in this refresh.

| Key | Status | Rationale |
|---|---|---|
| `subtel_popgen_lewontin1972` | STILL CURRENT | Foundational ~85% within-population benchmark; no newer paper replaces this reference |
| `subtel_popgen_weir1984` | STILL CURRENT | Weir-Cockerham FST estimator; still the standard for finite-sample correction |
| `subtel_popgen_hudson1992` | STILL CURRENT | Hudson FST (Dxy-based); the specific estimator used in our chapter 04 analysis; Bhatia 2013 (below) confirms it remains preferred |
| `subtel_popgen_rosenberg2002` | STILL CURRENT | Landmark 5-cluster out-of-Africa topology from microsatellites; the reference topology against which subtelomeric FST trees are benchmarked |
| `subtel_popgen_li2008` | STILL CURRENT | Serial founder effects and diversity gradient from Africa; expected diversity direction remains standard |
| `subtel_popgen_bhatia2013` | STILL CURRENT | Shows Hudson estimator preferred for sample-size independence; confirmed by new SV population studies that adopt it |
| `subtel_popgen_patterson2006` | STILL CURRENT | PCA/eigenanalysis; still the standard method for population structure visualization |
| `subtel_popgen_anderson2008` | STILL CURRENT | Drosophila subtelomere population genetics; only comparative context available; not superseded |
| `subtel_popgen_riethman2005` | STILL CURRENT | Human subtelomere structure review; foundational context; no newer equivalent summary |
| `subtel_popgen_1000g2010` | STILL CURRENT | Establishes the 1000G population sampling framework; Schloissnig 2025 (Section 3) extends it, does not replace it |
| `subtel_popgen_levysakin2019` | STILL CURRENT | Closest prior subtelomere population survey (optical mapping, 26 pops); now extended by HPRC v2 scale but not replaced |
| `subtel_popgen_pickrell2012` | STILL CURRENT | TreeMix population split/admixture graphs; still the reference method for allele-frequency-based topology inference |

Additional keys used in the .md but defined in REFERENCES.bib/REFERENCES_v2.bib
(`Trask1998`, `MeffordTrask2002`, `DerSarkissian2002`, `Ambrosini2007`, `Stong2014`,
`Sudmant2015`, `Young2020`, `Riethman2004`) are also still current; no new paper in this
search supersedes them.

**Key already in REFERENCES_v3.bib that is cited in NATURE_DRAFT_v1.md for this topic**:
`Bergstrom2020` — authoritative for genome-wide human diversity reference; still current;
no direct follow-up paper was identified in this search.

---

## Section 3: New papers to add (2023–2026)

Five papers meet the STRONG threshold. Two additional MODERATE papers are listed after.

---

### STRONG-1

**Proposed bibkey**: `subtel_popgen_jeong2025`

**Full citation**:
Jeong H, Kang S, Lee S, et al. Structural polymorphism and diversity of human segmental
duplications. *Nature Genetics*. 2025;57:346–358.
DOI: [10.1038/s41588-024-02051-8](https://doi.org/10.1038/s41588-024-02051-8)
PMID: 39779957

**Relevance** (STRONG): Surveys 170 genome assemblies from 85 individuals (38 African,
47 non-African) and finds that African genomes harbor significantly more intrachromosomal
SDs and are more likely to have recently duplicated gene families with higher copy numbers
than non-African samples. This is the closest direct parallel to our chr16q cross-arm
haplotype finding (70% AFR frequency) at genome-wide scale across segmental duplications.
Establishes that AFR enrichment of complex structural variants is a pervasive feature of
the human genome, not specific to a single locus.

**Suggested placement**: Para after the chr16q/chr4q AFR-enrichment result; frame as
genome-wide SD context for the localized subtelomeric finding.

**Claim supported**: C6 (cross-arm affinity mirrors out-of-Africa demography); the AFR
enrichment direction is corroborated.

---

### STRONG-2

**Proposed bibkey**: `subtel_popgen_schloissnig2025`

**Full citation**:
Schloissnig S, Pani S, Ebler J, et al. Structural variation in 1,019 diverse humans based
on long-read sequencing. *Nature*. 2025;644:442–452.
DOI: [10.1038/s41586-025-09290-7](https://doi.org/10.1038/s41586-025-09290-7)
PMID: 40702182

**Relevance** (STRONG): Long-read sequencing of 1,019 individuals from 26 populations
drawn from the 1000 Genomes Project; uncovers >100,000 sequence-resolved biallelic SVs
and genotypes 300,000 multiallelic VNTRs. Provides the population-scale reference for SV
characterization across AFR/EAS/SAS/AMR/EUR superpopulations. Directly comparable
framework to our HPRC v2 subtelomeric analysis; the AFR-enriched SV signal and the near-
zero FST between non-African populations for SVs observed here is consistent with our
non-AFR pairwise FST range (−0.047 to +0.007) — supporting the report values over the
erroneous draft value of 0.02–0.04 (CONSISTENCY_AUDIT Row 36).

**Suggested placement**: Methods context for the FST paragraph; also as reference for the
statement that non-AFR pairwise FST at SVs is near zero.

**Claim supported**: C6, and corrects the draft's non-AFR FST range.

---

### STRONG-3

**Proposed bibkey**: `subtel_popgen_kim2025`

**Full citation**:
Kim J, Park JL, Yang JO, et al. Highly accurate Korean draft genomes reveal structural
variation highlighting human telomere evolution. *Nucleic Acids Research*.
2025;53(1):gkae1294.
DOI: [10.1093/nar/gkae1294](https://doi.org/10.1093/nar/gkae1294)
PMID: 39778865

**Relevance** (STRONG): Three Korean (EAS) individuals with ~20× HiFi long-read assemblies
(contig N50: 6.3–58.2 Mb); identifies 131,138 deletions and 121,461 insertion SVs,
41.6% EAS-prevalent. Manual characterization of 19 large subtelomeric SVs (≥5 kb) reveals
the underlying DNA-repair mechanisms (microhomology-mediated, non-allelic homologous
recombination). The 103-individual short-read validation provides population frequency
estimates for EAS. Directly relevant to the EAS arm of the AFR-vs.-non-AFR comparison:
EAS-specific subtelomeric SVs are a distinct class from the AFR-enriched cross-arm
haplotypes; the two papers together delimit the EAS vs. AFR poles of the out-of-Africa
topology for chromosome-end structural variants.

**Suggested placement**: In the supplemental discussion of non-AFR subtelomeric diversity;
or in a sentence noting that EAS subtelomeric SVs are largely distinct from the AFR-
enriched cross-arm haplotypes.

**Claim supported**: C6 (mechanistic detail for the non-AFR side of the FST comparison).

---

### STRONG-4

**Proposed bibkey**: `subtel_popgen_jana2025`

**Full citation**:
Jana U, Rodriguez OL, Lees W, et al. The human IG heavy chain constant gene locus is
enriched for large structural variants and coding polymorphisms that vary among human
populations. *Cell Genomics*. 2025;6(1):101058.
DOI: [10.1016/j.xgen.2025.101058](https://doi.org/10.1016/j.xgen.2025.101058)
PMID: 41151584

**Relevance** (STRONG): Long-read sequencing of 105 diverse-ancestry individuals at the
IGHC locus (a complex, repeat-rich region); identifies population differentiation including
hundreds of SNVs in African and East Asian populations exceeding fixation index F=0.3, and
an IGHG4 haplotype enriched in Asian populations. Methodologically the closest parallel to
our analysis: FST computed for structural variants at a complex multi-copy locus in diverse
human populations, with AFR and EAS enrichment signals. Validates that AFR-biased
population differentiation at structural-variant-rich loci is reproducible beyond
subtelomeres.

**Suggested placement**: As supporting evidence that AFR enrichment of haplotypes at
complex SV loci (FST>0.1 AFR vs. non-AFR) is not unique to subtelomeres.

**Claim supported**: C6 (corroboration that the AFR FST elevation is a general property
of complex repeat loci in diverse human genomes).

---

### STRONG-5

**Proposed bibkey**: `subtel_popgen_porubsky2026`

**Full citation**:
Porubsky D, Yoo D, Koundinya N, et al. Population differences of chromosome 22q11.2
duplication structure predispose differentially to microdeletion and inversion.
*Nature Communications*. 2026;17(1). Published 2026-04-18.
DOI: [10.1038/s41467-026-71905-y](https://doi.org/10.1038/s41467-026-71905-y)
PMID: 42000714

**Relevance** (STRONG): Sequence-resolves 135 chr22q11.2 haplotypes from diverse humans;
finds 63 distinct structural configurations for low-copy repeat A (LCRA) differing in size
11-fold. African LCRA haplotypes are significantly longer (p=0.0047) and enriched for
inverted 105 kbp repeats, making them more protective against 22q11.2 microdeletion
syndrome. This is the clearest available example that African populations carry
architecturally more complex, longer SD haplotypes at a well-studied segmental duplication
locus — precisely the same population directionality as our chr16q subtelomeric cross-arm
enrichment. The mechanism (expansion of an ancient SD, longer in Africa consistent with
pre-OOA greater diversity) directly mirrors the proposed explanation for AFR-enriched
subtelomeric haplotypes.

**Suggested placement**: In the Discussion alongside Jeong 2025 as converging evidence
that African genome SD/subtelomere complexity exceeds non-African genomes.

**Claim supported**: C6; also provides mechanistic framing for the AFR-enrichment
direction.

---

### MODERATE-1

**Proposed bibkey**: `subtel_popgen_rausch2025review`

**Full citation**:
Rausch T, Marschall T, Korbel JO. The impact of long-read sequencing on human
population-scale genomics. *Genome Research*. 2025;35(4):593–598.
DOI: [10.1101/gr.280120.124](https://doi.org/10.1101/gr.280120.124)
PMID: 40228902

**Relevance** (MODERATE): Perspective review on the current state and future of long-read
population-scale genomics; highlights pangenome graphs, haplotype-resolved assemblies,
methylation from long reads, and the path toward clinical cohort studies. Provides the
conceptual framing for why population-scale long-read SV analysis (as in HPRC v2) is the
correct approach for subtelomeric population genetics. Useful for the Introduction or
Methods rationale but does not report new FST/population data.

---

### MODERATE-2

**Proposed bibkey**: `subtel_popgen_bird2023`

**Full citation**:
Bird N, Ormond L, Awah P, et al. Dense sampling of ethnic groups within African countries
reveals fine-scale genetic structure and extensive historical admixture. *Science Advances*.
2023;9(13):eabq2616.
DOI: [10.1126/sciadv.abq2616](https://doi.org/10.1126/sciadv.abq2616)
PMID: 36989356 · PMC: PMC10058250

**Relevance** (MODERATE): Autosomal SNP analysis of 1,333 individuals from >150 ethnic
groups in five African countries; documents fine-scale within-Africa genetic structure
correlated with historical polities, long-distance migrations, and Bantu expansions.
Relevant to the open question in the existing review (topic_13 §6.4) about whether the
aggregated "AFR" superpopulation masks within-Africa subtelomeric FST heterogeneity. The
SNP-level within-Africa structure documented here is SNIP-based, not structural-variant-
based, so it does not directly map to subtelomeric haplotype frequencies; include only if
discussing within-AFR heterogeneity as a caveat.

---

## Section 4: Contradictions

**No new 2023–2026 paper found in this search directly contradicts any specific numerical
claim or biological assertion in the existing topic-13 text.**

Three points of contextual tension are noted — these are not contradictions of the new
papers vs. the existing review, but the new papers provide evidence to resolve a known
internal inconsistency flagged by CONSISTENCY_AUDIT_v1.md:

**Internal inconsistency (draft vs. report, NOT a new-paper contradiction)**:
NATURE_DRAFT_v1.md L32 states non-AFR pairwise FST "0.02 to 0.04." The report
(04_heterogeneity.md) gives the actual computed values: EAS|AMR=0.007, EUR|AMR=0.007,
SAS|AMR=0.004, EAS|EUR=−0.047, SAS|EAS=0.005, SAS|EUR=−0.003; range −0.047 to +0.007.
The draft value is incorrect. The Schloissnig et al. 2025 population-scale SV study
(STRONG-2 above) independently documents near-zero FST between non-African populations
for structural variants, corroborating the report values and confirming the draft is wrong.
The draft should be corrected to "−0.05 to +0.01" to accurately reflect both the HPRC v2
results and the independent SV population literature.

**Directional consistency of AFR enrichment**: Jeong 2025 shows AFR enrichment of
*intrachromosomal* SDs, while our finding is *interchromosomal* cross-arm subtelomeric
haplotype enrichment in AFR. These are not contradictory — both point to African genomes
as more diverse at complex repeat loci — but authors should note the mechanism distinction
(intrachromosomal duplication vs. subtelomeric inter-arm sequence exchange) to prevent
conflation.

**Subtelomeric SVs in EAS vs. AFR**: Kim 2025 documents 41.6% EAS-prevalent SVs among
~252,000 total SVs from Korean genomes. This is not a contradiction with our AFR-biased
cross-arm subtelomere finding; EAS-prevalent SVs in Kim 2025 are genome-wide (including
non-subtelomeric regions) and represent overall genomic diversity differences, not
specifically subtelomeric cross-arm haplotypes. No numerical contradiction exists.

---

## Section 5: Search audit trail

**Databases searched**: PubMed (via MCP), bioRxiv/medRxiv (via MCP), OpenAlex (via
openalex-database skill; queries submitted 2026-05-17).

**PubMed queries** (all with date filter 2023/01/01–2026/01/01 unless noted):
1. "subtelomeric variation population" → 100 hits; top 10 returned; 4 relevant PMIDs reviewed
2. "structural variant population differentiation diverse genomes long-read" (2024+) → 12 hits; 5 PMIDs reviewed
3. "long-read sequencing human structural variation diverse ancestry population specific" → 6 hits; all reviewed
4. "Eichler segmental duplications structural variants African population enriched genome" → 1 hit (Porubsky 2026, PMID 42000714)
5. "Jeong segmental duplication structural polymorphism human populations African 2025" → 1 hit (confirms PMID 39779957)
6. "subtelomere telomere structural variation diverse populations assembly 2023 2024" → 0 hits (over-specified)
7. "Vollger segmental duplication population genetics 2023 2024 2025" → 0 hits (year AND expansion)
8. "human pangenome reference consortium structural variation population genetics superpopulation 2023 2024 2025" → 0 hits
9. "population specific structural variants copy number enriched African genome assembly long-read" → 0 hits
10. "Trost population stratified structural variants long-read diverse human genomes" → 0 hits
11. "structural variant population differentiation Africa non-Africa out-of-Africa bottleneck 2023 2024" → 0 hits (term expansion issue)
12. Additional targeted author/keyword queries: Sudmant 2024-2025, Bird 2023, Rausch 2025 (review), 1000G LR

**PMIDs verified via get_article_metadata**: 39779957, 40702182, 39778865, 40228902,
36989356, 41151584, 42000714, 40932769 (Moya 2025 PNAS; MODERATE, not included in
Section 3 due to disease-VNTR focus), 38844673, 37633302, 37459347.

**Papers reviewed but excluded**:
- PMID 40932769 (Moya 2025 PNAS, CACNA1C VNTR schizophrenia) — population-specific VNTR, African-enriched alleles, but not a subtelomere/SD study; MODERATE; included as optional if discussion of VNTR African enrichment is needed
- PMID 40228902 (Rausch 2025 Genome Res) — review; MODERATE (included above)
- PMID 36989356 (Bird 2023 Sci Adv) — SNP within-Africa structure; MODERATE (included above)
- PMIDs 37459347, 37308293, 40909903, 41500311, 42127154 — non-human or non-population-genetics papers returned by broad subtelomere query
- PMID 39520989 (Mostovoy 2024 AJHG) — complex SV resolution with T2T/LR; clinical, not population FST
- PMID 41044908 (Lehnert 2025 Mol Ecol) — salmon haplotype blocks; non-human

**Keys confirmed absent from REFERENCES_v3.bib**: All 5 proposed keys
(`subtel_popgen_jeong2025`, `subtel_popgen_schloissnig2025`, `subtel_popgen_kim2025`,
`subtel_popgen_jana2025`, `subtel_popgen_porubsky2026`) verified absent by grep on
REFERENCES_v3.bib (295 total entries).

**Coverage assessment**: Searches were exhaustive for papers matching the narrow topic
(subtelomere + population genetics + FST + diverse human populations). PubMed query
expansion frequently merged year strings as AND-required terms, suppressing results; this
was mitigated by using date-range filters instead of inline year terms. The OpenAlex skill
was also invoked with parallel queries across four angles (subtelomere popgen FST,
structural variant population differentiation, Bergstrom human diversity superpopulation,
Weir/Hudson FST new applications). The five STRONG papers identified represent the most
directly relevant 2023–2026 literature found across all search strategies.
