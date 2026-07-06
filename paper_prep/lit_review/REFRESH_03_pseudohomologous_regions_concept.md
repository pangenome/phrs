# Literature Refresh: Topic 03 — Pseudohomologous Regions Concept

**Agent:** lit-refresh-03-pseudohomologous (agent-55)
**Date:** 2026-05-17
**Scope window:** 2023–2026

---

## Section 1: Topic scope

Topic 03 covers the intellectual lineage of the pseudohomologous region (PHR) concept: from early telomere-associated polymorphism (Brown & Wilkie 1990–1991) through the two-domain Flint/Mefford architectural model (1997–2002), the sequence-block catalogue (Linardopoulou 2005, Ambrosini 2007), and on to the T2T/HPRC enabling literature (Logsdon 2021, Nurk 2022, Liao 2023). It also anchors the best-characterized PHR pair — 4q/10q D4Z4 — via the clinical FSHD literature (Lemmers 2002, Ma 2024) and validates ongoing exchange conceptually (Rudd 2007, 2009). In the Nature draft this topic underpins paragraphs P2–P3 of the main text (L22–L28): the historical justification for why inter-chromosomal subtelomeric sharing exists and why a pangenome-scale, reference-free approach was needed. Draft claims directly traceable to topic 03 are **C4** (extended interchromosomal homology at nearly all subtelomeres), **C5** (named community clades, especially 4q/10q C1 and 10p/18p C2), and **C8** (ongoing exchange and concerted evolution). The end-to-end report sections most relevant are §01 (pipeline, PHR detection rationale) and §07 (integrated 3D interpretation, D4Z4-CTCF-lamin model for C1). In the slide deck (`zoom_review_deck.typ`) the concept appears in PHR definition slides (l.257–284), the community-method slide (l.370–376), and the PHR-length and community-labelling slides (l.537–595).

---

## Section 2: Existing citations still authoritative

All 14 entries from `topic_03_pseudohomologous_regions_concept.bib` assessed below.

| Bibkey | Status | Reason |
|---|---|---|
| Brown1990 | STILL CURRENT | Foundational; earliest description of polymorphic telomere-associated DNA in humans. No newer paper replaces its historical priority. |
| Wilkie1991 | STILL CURRENT | 260 kb 16p length polymorphism; still the canonical early-era evidence for subtelomeric exchange. |
| Martin2002 | STILL CURRENT | Evolutionary origin of subtelomeric homologies via primate chromosome history; no newer synthesis supersedes this framing. |
| Fan2002 | STILL CURRENT | Chr2 ancestral fusion site paralogy; still cited routinely in chr2-fusion papers (e.g., Poszewiecka 2023, Yang 2025). |
| Rudd2007 | STILL CURRENT | Sister-chromatid exchange rates at chromosome ends; mechanism remains cited as support for terminal recombinogenicity. |
| Rudd2009 | STILL CURRENT | Primate fission comparative analysis; best available non-homologous repair and retrotransposition evidence at subtelomeres. |
| Young2020 | STILL CURRENT | Optical mapping of 154 genomes; still the broadest pre-HPRC survey of large-scale subtelomeric structural variation. |
| Grigorev2021 | STILL CURRENT | Telomeric repeat-motif haplotype diversity; long-read evidence for subtelomeric heterogeneity predating T2T. |
| Logsdon2021 | STILL CURRENT | First complete chr8; enabling literature for T2T claims; still canonical. |
| Nurk2022 | STILL CURRENT | T2T-CHM13; foundational enabling literature; cannot be superseded in role. |
| Liao2023 | STILL CURRENT | HPRC draft pangenome; enabling literature for HPRC v1 context; companion HPRC v2 paper is the direct successor but does not replace this. |
| Lemmers2002 | STILL CURRENT | 4qA haplotype FSHD association; remains the primary citation for the polyadenylation signal distinction. Newer FSHD papers (Tardy 2026, Delourme 2023) extend rather than replace. |
| Vollger2022 | STILL CURRENT | T2T-CHM13 segmental duplication inventory; 68 Mbp previously unresolved SDs; still the definitive SD census from CHM13. |
| Ma2024 | STILL CURRENT | CRISPR-induced 4q/10q subtelomeric exchange reverting FSHD genotype; 2024 is the newest entry and remains the only experimental proof of therapeutic PHR exchange. |

No entry is SUPERSEDED. Three entries (Lemmers2002, Vollger2022, Ma2024) now have new companion papers (see Section 3) that deepen but do not replace the original claims.

---

## Section 3: NEW papers to add (STRONG relevance, 2023+)

---

**[Salsi2026fshd]** Salsi V, Losi F, Pini S, Chiara M, Tupler R. 2026. Rethinking genomics of facioscapulohumeral muscular dystrophy in the telomere-to-telomere era: pitfalls in the hidden landscape of D4Z4 repeats. *European Journal of Human Genetics* 34(3):357–367. PMID:41535478 / DOI:10.1038/s41431-025-02000-x

**Claim it supports/contradicts:** Draft P2, L22: "The D4Z4 macrosatellite of 4q35 is found in degenerate copies on 10q26 and underlies facioscapulohumeral dystrophy." Also report §07 D4Z4-CTCF-lamin model: "C1 sequences carry median 22 DUX4L genes while all 7 non-C1 outlier sequences have 0–2 on their own arm."

**Extends the claim.** Salsi et al. used T2T-CHM13 v2.0 to systematically reannotate D4Z4-like repeats genome-wide. They found D4Z4-like loci on **at least ten additional chromosomes** beyond 4q35 and 10q26, several carrying intact DUX4 open reading frames or polyadenylation signals. Standard primer sets (used for DUX4/DBE-T detection) amplify paralogous arrays on these other chromosomes. This directly supports the draft's C1 community framing while also warning that copy-number estimates relying on short-read or hg38 references may misattribute signal: the current paper's DUX4L census (median 22 per C1 sequence) depends on T2T-aware mapping that Salsi et al. show is necessary.

**Recommendation:** ADD. Cite alongside Lemmers2002 and Ma2024 in P2 (4q/10q D4Z4 background) to acknowledge that the "copies on 10q26" statement from legacy literature is now understood as the tip of a broader genome-wide D4Z4 paralogy landscape visible only in T2T assemblies.

```bibtex
@article{Salsi2026fshd,
  author  = {Salsi, Valentina and Losi, Francesca and Pini, Sara and
             Chiara, Matteo and Tupler, Rossella},
  title   = {Rethinking genomics of facioscapulohumeral muscular dystrophy
             in the telomere-to-telomere era: pitfalls in the hidden
             landscape of {D4Z4} repeats},
  journal = {European Journal of Human Genetics},
  year    = {2026},
  volume  = {34},
  number  = {3},
  pages   = {357--367},
  doi     = {10.1038/s41431-025-02000-x},
  pmid    = {41535478}
}
```

---

**[Tardy2026fshd]** Tardy C, Trani JP, Murcia Pienkowski V, et al. 2026. Benchmarking long-read sequencing approaches to resolve facioscapulohumeral dystrophy locus complexity. *Brain* 149(5):1750–1767. PMID:41642686 / DOI:10.1093/brain/awag029

**Claim it supports/contradicts:** Draft main text L38: "16 are crossover-like: a single reciprocal exchange in which the query haplotype switches its source chromosome at a discrete breakpoint." Also report §07 D4Z4-CTCF-lamin model, and Ma2024 framing in topic_03 narrative: "exchange between these arms is not merely historical residue but a documented ongoing process."

**Directly supports C8.** Tardy et al. resolved 7 complex structural variants (duplicated alleles) at 4q35/10q26 using Oxford Nanopore and PacBio long-read sequencing, characterizing breakpoints at nucleotide resolution and linking methylation patterns to pathogenicity. Key finding: duplicated alleles arise from **intrachromosomal recombination between LSau elements within D4Z4 and distal subtelomeric β-satellite elements**, with patient-specific breakpoints. This provides mechanistic detail for how 4q/10q subtelomeric exchange leaves structural rearrangements — exactly the class of exchange the draft calls "crossover-like" at population scale. It also confirms that standard optical genome mapping fails to detect these variants, reinforcing the draft's claim that long-read, assembly-level approaches are required for PHR analysis.

**Recommendation:** ADD. Cite in P2/methods discussion of 4q/10q exchange mechanics, alongside Ma2024 and Delourme2023fshd. Particularly useful for the "ongoing exchange with measurable medical consequences" sentence in the draft.

```bibtex
@article{Tardy2026fshd,
  author  = {Tardy, Charlotte and Trani, Jean Philippe and
             Murcia Pienkowski, Victor and Morin, Loeva and
             Castro, Christel and Souville, Louis and
             Humbert, Camille and G{\'e}rard, Laur{\`e}ne and
             Eudes, Nathalie and Assoumani, Amire and
             Bertaux, Karine and Verebi, Camille and
             Nectoux, Juliette and Salort Campana, Emmanuelle and
             Jacquemont, Marie Line and Toutain, Annick and
             Mallaret, Martial and Tard, C{\'e}line and
             Fradin, M{\'e}lanie and Attarian, Shahram and
             Nguyen, Karine and Bernard, Rafa{\"e}lle and
             Magdinier, Fr{\'e}d{\'e}rique},
  title   = {Benchmarking long-read sequencing approaches to resolve
             facioscapulohumeral dystrophy locus complexity},
  journal = {Brain},
  year    = {2026},
  volume  = {149},
  number  = {5},
  pages   = {1750--1767},
  doi     = {10.1093/brain/awag029},
  pmid    = {41642686}
}
```

---

**[Delourme2023fshd]** Delourme M, Charlene C, Gerard L, et al. 2023. Complex 4q35 and 10q26 Rearrangements: A Challenge for Molecular Diagnosis of Patients With Facioscapulohumeral Dystrophy. *Neurology Genetics* 9(3):e200076. PMID:37200893 / DOI:10.1212/NXG.0000000000200076

**Claim it supports/contradicts:** Topic_03 review (existing): "Van Deutekom et al. established in 1996 that D4Z4 tandem-repeat units are exchanged between 4q35 and 10q26 at approximately 20% frequency." Also draft C5 (4q/10q community) and C8 (ongoing exchange).

**Directly supports C5 and C8.** Delourme et al. investigated 2,363 FSHD diagnostic cases using molecular combing (MC). They found 147 individuals (6.2%) carrying atypical organization of 4q35 or 10q26 loci, including 4q-10q translocations, D4Z4 array duplications (1–2% of cases), and somatic mosaicism. In 54 of these cases the structural rearrangement was the only genetic defect, suggesting pathogenicity. This large clinical cohort provides population-scale evidence (beyond the original van Deutekom 1996 Dutch cohort) that 4q/10q interchromosomal exchange is ongoing and medically consequential. It also highlights that optical genome mapping (Bionano) fails to resolve these complex rearrangements — consistent with the draft's message that reference-free long-read approaches are required.

**Recommendation:** ADD. Cite alongside Ma2024 in P2 (4q/10q D4Z4 background) as contemporary clinical-scale evidence for ongoing interchromosomal exchange at the C1 PHR community. Also useful for footnoting the "medically consequential" framing.

```bibtex
@article{Delourme2023fshd,
  author  = {Delourme, Meg{\"a}ne and Chaix, Charlene and
             G{\'e}rard, Lauren{\`e} and Ganne, Benjamin and
             Perrin, Pierre and Vovan, Catherine and
             Bertaux, Karine and Nguyen, Karine and
             Bernard, Rafa{\"e}lle and Magdinier, Fr{\'e}d{\'e}rique},
  title   = {Complex 4q35 and 10q26 Rearrangements: {A} Challenge for
             Molecular Diagnosis of Patients With Facioscapulohumeral
             Dystrophy},
  journal = {Neurology Genetics},
  year    = {2023},
  volume  = {9},
  number  = {3},
  pages   = {e200076},
  doi     = {10.1212/NXG.0000000000200076},
  pmid    = {37200893}
}
```

---

**[Zhuang2026dux4]** Zhuang Z, Ueda MT, Yamaguchi K, Kochi Y. 2026. A new integrated genetic and transcriptomic approach for investigating DUX4 and DUX4C. *Journal of Human Genetics* [Epub ahead of print]. PMID:41540238 / DOI:10.1038/s10038-025-01450-x

**Claim it supports/contradicts:** Topic_03 narrative: "the FSHD literature … ectopic exchange between chr4_q and chr10_q D4Z4 arrays can modify FSHD alleles." Also draft reference to haplotype-level diversity in the 4q subtelomere.

**Supports haplotype complexity at 4q subtelomere.** Zhuang et al. identified two distinct DUX4C haplotypes (4qα and 4qβ) with eQTL effects, generated D4Ref-T2T (a T2T long-read reference for 4q subtelomere), and found strong linkage disequilibrium between DUX4C and DUX4 haplotypes (r = 0.86). This is directly relevant to the draft's use of C1 community, which has wide copy-number diversity (0–22 DUX4L copies across 465 near-complete assemblies). The paper underscores that even within the 4q subtelomere, haplotype-level complexity at T2T resolution is only now being resolved — consistent with the draft's message that population-scale pangenome analysis reveals structure missed by earlier references.

**Recommendation:** ADD. Cite in P2 alongside Lemmers2002 when discussing 4q subtelomere haplotype structure. Also cite in the methods discussion of DUX4L copy-count diversity.

```bibtex
@article{Zhuang2026dux4,
  author  = {Zhuang, Zhaohui and Ueda, Mahoko Takahashi and
             Yamaguchi, Kensuke and Kochi, Yuta},
  title   = {A new integrated genetic and transcriptomic approach for
             investigating {DUX4} and {DUX4C}},
  journal = {Journal of Human Genetics},
  year    = {2026},
  doi     = {10.1038/s10038-025-01450-x},
  pmid    = {41540238}
}
```

---

**[Kanoh2023subtel]** Kanoh J. 2023. Subtelomeres: hotspots of genome variation. *Genes & Genetic Systems* 98(3):155–160. PMID:37648502 / DOI:10.1266/ggs.23-00049

**Claim it supports/contradicts:** Topic_03 narrative: "Linardopoulou et al. explicitly calling human subtelomeres hot spots of interchromosomal recombination and segmental duplication." Also draft P2 framing of subtelomeres as dynamic exchange substrates.

**MEDIUM — provides current synthesis.** Kanoh 2023 is a peer-reviewed review explicitly linking subtelomere sequence variation to eukaryotic biology (including yeasts and mammals), with recent discoveries on copy-number variation and exchange. It is a concise (6-page) up-to-date synthesis that supports C4 and C8 as background, and could serve as a modern secondary citation alongside Linardopoulou2005 and Ambrosini2007 where the draft invokes the "hotspot" framing. It does not provide new human pangenome data.

**Recommendation:** ADD as background context citation. Useful when first invoking the subtelomere-hotspot concept in P2, alongside the historical Linardopoulou2005.

```bibtex
@article{Kanoh2023subtel,
  author  = {Kanoh, Junko},
  title   = {Subtelomeres: hotspots of genome variation},
  journal = {Genes \& Genetic Systems},
  year    = {2023},
  volume  = {98},
  number  = {3},
  pages   = {155--160},
  doi     = {10.1266/ggs.23-00049},
  pmid    = {37648502}
}
```

---

**[Kim2025korean]** Kim J, Park JL, Yang JO, et al. 2025. Highly accurate Korean draft genomes reveal structural variation highlighting human telomere evolution. *Nucleic Acids Research* 53(1). PMID:39778865 / DOI:10.1093/nar/gkae1294

**Claim it supports/contradicts:** Topic_03 narrative reference to Young2020 (optical mapping of 154 genomes showing widespread large-scale subtelomeric structural variation). Draft P2: subtelomeric variation "cannot be inferred reliably from short-read references alone."

**MEDIUM — population-scale subtelomeric SVs.** Kim et al. generated ~20× HiFi long-read assemblies for three Korean individuals and characterized 19 large (≥5 kb) subtelomeric SVs, resolving the underlying DNA damage-repair mechanisms (non-allelic homologous recombination, microhomology-mediated end joining, etc.) at nucleotide resolution. They found 41.6% of SVs prevalent in the East Asian population. This extends Young2020's optical-mapping view with repair-mechanism detail and provides a non-European population complement. Directly supports the statement that subtelomeric variation requires long-read, assembly-level analysis (C2 methodology justification).

**Recommendation:** ADD alongside Young2020 when justifying the limitation of reference-based approaches for subtelomeric regions (draft P2, report §10 limitations).

```bibtex
@article{Kim2025korean,
  author  = {Kim, Jun and Park, Jong Lyul and Yang, Jin Ok and
             Kim, Sangok and Joe, Soobok and Park, Gunwoo and
             Hwang, Taeyeon and Cho, Mun-Jeong and Lee, Seungjae and
             Lee, Jong-Eun and Park, Ji-Hwan and Yeo, Min-Kyung and
             Kim, Seon-Young},
  title   = {Highly accurate {Korean} draft genomes reveal structural
             variation highlighting human telomere evolution},
  journal = {Nucleic Acids Research},
  year    = {2025},
  volume  = {53},
  number  = {1},
  doi     = {10.1093/nar/gkae1294},
  pmid    = {39778865}
}
```

---

## Section 4: CONTRADICTIONS

No papers found in the 2023–2026 search window that contradict the core claims of topic 03 (C4: extended interchromosomal subtelomeric homology; C5: named community clades; C8: ongoing exchange). The strongest potential challenge is **Salsi et al. 2026**, which reveals that D4Z4-like loci on ≥10 chromosomes beyond 4q/10q may have been mis-attributed to 4q/10q in prior short-read studies — but this does not contradict the PHR community claim; it extends it (more chromosomes share D4Z4-like sequence than previously thought, consistent with C1 being embedded in a broader network). Similarly, **Tardy et al. 2026** reveals that complex structural variants at 4q35/10q26 are more common than previously diagnosed, but this is consistent with, not contradictory to, the draft's characterization of ongoing exchange.

**Potential nuance (not a contradiction):** Salsi et al. 2026 raises a diagnostic caveat that PCR-based D4Z4 copy counting may conflate signal from paralogous loci on other chromosomes. The draft should not claim that 4q/10q is the only location of D4Z4-like sequence; it should say "the canonical D4Z4-bearing PHR pair" and note that T2T reveals broader paralogy.

---

## Section 5: Search audit trail

### Tools used

1. **mcp__claude_ai_PubMed__search_articles** — primary search engine. Date filter applied: `date_from=2023`, `date_to=2026` (pdat mode) for all queries.
2. **mcp__claude_ai_PubMed__get_article_metadata** — used to fetch full metadata for shortlisted PMIDs.
3. **mcp__claude_ai_bioRxiv__search_preprints** — attempted; API does not support keyword search (date + category only). Genomics category search returned 100 results spanning 2023-2025 but none were keyword-matchable to PHR topic via category alone; no new hits.
4. **OpenAlex API** (openalex-database skill, Python client) — three keyword search rounds: "pseudohomologous region subtelomeric duplication interchromosomal", "Linardopoulou subtelomere chromosome homology", "subtelomere hotspot recombination duplication human chromosome", "Vollger subtelomere segmental duplication T2T pangenome". Identified candidate DOIs for complex 4q/10q rearrangements, Kanoh 2023 review, and Kim 2025 Korean SV paper.

### Query strings (PubMed)

| Query | Hits |
|---|---|
| `pseudohomologous region subtelomeric` 2023–2026 | 0 |
| `subtelomeric segmental duplication pangenome human` 2023–2026 | 5 |
| `subtelomere assembly evolution` 2023–2026 | 23 |
| `subtelomere FSHD DUX4` 2023–2026 | 4 |
| `Complex 4q35 10q26 Rearrangements subtelomeric` 2023–2026 | 2 |
| `subtelomere hotspot genome variation review` 2023–2026 | 1 |
| `FSHD telomere-to-telomere T2T D4Z4 genomics` 2023–2026 | 1 |
| `non-allelic homologous recombination subtelomere human chromosome rearrangement` 2023–2026 | 1 |
| `meiotic bouquet telomere clustering recombination nuclear envelope` 2023–2026 | 3 |
| `Zhuang DUX4 DUX4C haplotype T2T long-read subtelomere` 2024–2026 | 1 |
| `subtelomere telomere variation structure population human long read 2024` | 1 |

### Total hits before relevance filter

Combined unique PubMed PMIDs retrieved across all queries: **~45** (many overlapping across queries).

### Hits after relevance filter

**STRONG (added to Section 3):** 4 papers — Salsi2026fshd, Tardy2026fshd, Delourme2023fshd, Zhuang2026dux4.
**MEDIUM (added to Section 3):** 2 papers — Kanoh2023subtel, Kim2025korean.

### Papers dropped and why

| Paper | PMID / DOI | Reason dropped |
|---|---|---|
| Yang et al. 2025, Cell Genomics — chr2 fusion site SDs | 41338219 | About chr2 fusion-specific ILS, not general PHR concept; tangential. |
| Poszewiecka et al. 2023, Genome Biology — PhaseDancer assembler | 37697406 | Tool paper for chr2/great ape SD assembly; WEAK for PHR concept. |
| Lee 2026, Mobile DNA — SEPTIN14P pseudogene propagation | 41699652 | About secondary propagation of processed pseudogenes; topic-03-adjacent but not direct PHR concept claim. |
| Brann et al. 2024, BMC Genomics — Schistosoma mansoni subtelomere | 38413905 | Non-human parasite; irrelevant to human PHR concept. |
| Bakewell et al. 2025, Leukemia — DUX4-rearranged B-ALL | 40940582 | Cancer rearrangement biology; does not address subtelomeric PHR concept claims. |
| Attarian et al. 2024, J Neurol — FSHD French national protocol | 38955828 | Clinical management review; no new PHR concept insight. |
| Fernández-Álvarez 2023, Front Cell Dev Biol — bouquet non-canonical functions | 38020928 | Bouquet biology is relevant to topic 08 (meiotic bouquet) not topic 03 (PHR concept). |
| Zhang et al. 2025, bioRxiv — rhesus macaque T2T | 41019632 | SATR satellite subtelomeres in macaque; not human PHR. |
| You et al. 2024, Plant Physiol — rice chromosome end pairing | 38478471 | Plant species; irrelevant. |
| Various plant/fungal subtelomere papers | multiple | Non-human organisms; WEAK for topic 03. |
