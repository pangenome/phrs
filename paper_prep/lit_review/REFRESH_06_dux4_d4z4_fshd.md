# REFRESH_06: dux4_d4z4_fshd

Generated: 2026-05-17  
Agent: agent-51 (lit-refresh-06-dux4)  
Source searches: PubMed MCP (2023-2026), bioRxiv MCP (date-range genetics; no keyword filter available)  
Spot-verified PMIDs: 41535478, 37798116, 37703328 (all confirmed via mcp__claude_ai_PubMed__get_article_metadata)

---

## Section 1: Topic scope

Topic 06 covers the D4Z4 macrosatellite repeat array on chromosome 4q35, its near-identical paralog on chromosome 10q26, and the disease mechanism of facioscapulohumeral muscular dystrophy (FSHD). In our materials this topic anchors three draft claims: **C5** (the NJ tree recovers a 4q/10q DUX4 clade with wide copy-number diversity, corresponding to Leiden community C1), **C7** (D4Z4 arrays at 4q35 and 10q26 are tethered to the nuclear lamina via CTCF and lamin A/C, explaining peripheral Dip-C radial position 0.732), and **C8** (ongoing interchromosomal exchange is directly observed, e.g., the chr10q←chr4q gene-conversion event in PAN028). In the draft (NATURE_DRAFT_v1.md), topic-06 citations appear in paragraph 1 (NJ tree/historical framing, citing `dux4_d4z4_fshd_lemmers2010worldwide` and `dux4_d4z4_fshd_lemmers2007`), in the 3D-proximity paragraph (citing `Ottaviani2009`, `OttavianiGilson2008`, `Masny2004`, `Cabianca2012`), and in the causal-loop closing paragraph. In the end-to-end report, sections 03 (gene enrichment: C1 DUX4L pseudogenes, CTCF insulator, pLAM model) and 07 (integrated: D4Z4-CTCF-lamin tethering model, C1 co-peripheralisation, CTCF/cohesin testable prediction) carry the primary topic-06 content. The NJ-tree slide in `zoom_review_deck.typ` labels the 4q/10q branch and names D4Z4/DUX4 copy-number diversity as the clade marker.

---

## Section 2: Existing citations still authoritative

All bibkeys from `topic_06_dux4_d4z4_fshd.bib` evaluated. Four keys listed in the bib file header as "in REFERENCES.bib" (`vanDeutekom1996`, `Stout1999`, `Lemmers2010`, `Masny2004`, `Ottaviani2009`) are present in `REFERENCES_v3.bib` and are also evaluated.

| Bibkey | Status | Reason |
|--------|--------|--------|
| `dux4_d4z4_fshd_wijmenga1992` | STILL CURRENT | Foundational 4q35 linkage mapping; no paper supersedes the original chromosomal assignment |
| `dux4_d4z4_fshd_hewitt1994` | STILL CURRENT | First full D4Z4 repeat sequence; still the primary citation for unit architecture. Caveat: salsi2026t2t now shows additional D4Z4-like loci across the genome using T2T, but Hewitt 1994 remains the structural reference |
| `dux4_d4z4_fshd_vanGeel2002` | STILL CURRENT | Primary citation for 4q/10q comparative genomics (>42 kb homology, >99% D4Z4 unit identity) |
| `dux4_d4z4_fshd_lemmers2002` | STILL CURRENT | FSHD specificity of 4qA haplotype; no replication challenge found |
| `dux4_d4z4_fshd_lemmers2007` | STILL CURRENT | Specific sequence variants within 4q35 associated with FSHD; cited directly in draft P1 |
| `dux4_d4z4_fshd_snider2009` | STILL CURRENT | D4Z4 transcript complexity; recent nanopore data (xiao2026d4z4, butterfield2023nanopore) extends but does not contradict |
| `dux4_d4z4_fshd_degreef2010` | STILL CURRENT | FSHD1/FSHD2 symmetric vs asymmetric hypomethylation; foundational for epigenetic model |
| `dux4_d4z4_fshd_lemmers2010worldwide` | STILL CURRENT | Four interchromosomal transfer events in human evolution; cited in draft P1 and abstract. New population-scale data (HPRC v2) is consistent, not contradictory |
| `dux4_d4z4_fshd_cabianca2012` | STILL CURRENT | DBE-T lncRNA at 4q35; salsi2025genotoxic extends the RNA biology but does not supersede Cabianca |
| `dux4_d4z4_fshd_geng2012` | STILL CURRENT | DUX4 target gene landscape; germline/DEFB103/retroelement pathways still authoritative |
| `dux4_d4z4_fshd_lemmers2012` | STILL CURRENT | SMCHD1/FSHD2 digenic model; sikrova2023smchd1 adds mechanistic detail to LRIF1 crosstalk but does not supersede |
| `dux4_d4z4_fshd_hendrickson2017` | STILL CURRENT | DUX4 in ZGA; still the primary reference for embryonic context. Recent work (Gao2026phasesep, hamm2026klf18) adds downstream factors but does not challenge the core finding |
| `dux4_d4z4_fshd_ma2024` | STILL CURRENT | CRISPR 4q/10q subtelomere exchange reverts FSHD phenotype; unique proof-of-concept, not yet replicated or challenged. **Note:** this key duplicates `Ma2024` in REFERENCES_v3.bib (same PMID 38691604, DOI 10.1126/sciadv.adl1922); when assembling REFERENCES_v4.bib, retain the `Ma2024` key and drop the prefixed duplicate |
| `vanDeutekom1996` (in REFERENCES_v3.bib) | STILL CURRENT | First 4q35 gene identification (FRG1); historical, not a primary claim anchor |
| `Stout1999` (in REFERENCES_v3.bib) | STILL CURRENT | Dual-color FISH for 4q/10q hybrids in 20% of population; cited in topic_06 review |
| `Masny2004` (in REFERENCES_v3.bib) | STILL CURRENT | Lamin A/C tethering of 4q35 at nuclear periphery; cited in draft P6 for C7 |
| `Ottaviani2009` (in REFERENCES_v3.bib) | STILL CURRENT | CTCF/lamin insulator at D4Z4; cited in draft P6 and report 07. butterfield2023nanopore and xiao2026d4z4 provide single-CpG resolution validation of the CTCF site |

---

## Section 3: NEW papers to add (STRONG relevance, 2023+)

---

**[salsi2026t2t]** Salsi V et al., 2026. Rethinking genomics of facioscapulohumeral muscular dystrophy in the telomere-to-telomere era: pitfalls in the hidden landscape of D4Z4 repeats. *European Journal of Human Genetics* 34(3):357-367. PMID:41535478 / DOI:10.1038/s41431-025-02000-x

Claim it supports/qualifies: Draft abstract: "the D4Z4 macrosatellite of 4q35 is found in degenerate copies on 10q26" (P1); report 03_gene_enrichment.md: "DUX4L pseudogenes (22), AGGF1P1, CLUHP4, DBET" (C1 annotation).

Using T2T-CHM13 v2.0, this paper annotates the full D4Z4-like repertoire across the genome and finds clusters and monomers on at least ten chromosomes beyond 4q/10q, several harboring intact DUX4 open reading frames or polyadenylation signals. In silico PCR shows that widely used DUX4 and DBE-T primer sets amplify multiple off-target loci. The draft's formulation "found in degenerate copies on 10q26" understates the genome-wide D4Z4-like distribution; the 4q/10q community C1 reflects the dominant arms but the genome contains additional loci. This does not contradict the community analysis but is relevant to the CTCF/insulator prediction and to any DUX4 expression assays cited as support.

Recommendation: ADD to paragraph P1 (framing) as a FLAG/qualification: the binary 4q/10q model is incomplete in the T2T era; the community C1 still reflects the dominant 4q/10q arms, but additional DUX4-capable loci exist. Also relevant to report 07 "testable prediction" note about CTCF ChIP-seq realignment.

```bibtex
@article{salsi2026t2t,
  author    = {Salsi, Valentina and Losi, Francesca and Pini, Sara and Chiara, Matteo and Tupler, Rossella},
  title     = {Rethinking genomics of facioscapulohumeral muscular dystrophy in the telomere-to-telomere era: pitfalls in the hidden landscape of {D4Z4} repeats},
  journal   = {European Journal of Human Genetics},
  year      = {2026},
  volume    = {34},
  number    = {3},
  pages     = {357--367},
  doi       = {10.1038/s41431-025-02000-x},
  pmid      = {41535478}
}
```

---

**[xiao2026d4z4]** Xiao LC et al., 2026. Complete genetic and epigenetic architecture of D4Z4 macrosatellites in FSHD, BAMS, and reference cohorts with D4Z4End2End. *Genome Research* 36(4):827-848. PMID:41871882 / DOI:10.1101/gr.280907.125

Claim it supports: Report 07_integrated.md L64-67: "CTCF binds within D4Z4 repeat units (Ottaviani et al. 2009); D4Z4-proximal sequences are tethered to the nuclear periphery via lamin A/C interaction (Masny et al. 2004; Ottaviani et al. 2009)"; draft P6: "The D4Z4 array on 4q35 binds CTCF and is tethered to the nuclear lamina." Also supports C8 (structural variant diversity from recombination).

D4Z4End2End uses ultra-long whole-genome and Cas9-targeted nanopore sequencing to achieve full-span sequencing of arrays up to 40 repeat units (~132 kb). The study reveals: (1) length- and SMCHD1-dependent methylation gradients across the array, reaching a threshold at ~10 repeats consistent with the known FSHD1 boundary; (2) in-cis duplications visible as discrete structural variants within the array; (3) allele-level resolution of methylation patterns in FSHD1, FSHD2, and BAMS (Bosma arhinia microphthalmia syndrome, also SMCHD1-associated) cohorts. This provides the highest-resolution characterization of the D4Z4 epigenetic landscape to date, directly supporting the CTCF-methylation axis in claim C7.

Recommendation: ADD to support C7; also supports the "D4Z4 array architecture at nucleotide resolution" section of topic_06 review's open questions, which specifically anticipated this class of result.

```bibtex
@article{xiao2026d4z4,
  author    = {Xiao, Lucinda C and Semwal, Ayush and St John, Brianna and Zeglinski, Kathleen and Su, Shian and Lancaster, James and Xue, Shifeng and Reversade, Bruno and Ritchie, Matthew E and Magdinier, Fr{\'e}d{\'e}rique and Blewitt, Marnie E and Gouil, Quentin},
  title     = {Complete genetic and epigenetic architecture of {D4Z4} macrosatellites in {FSHD}, {BAMS}, and reference cohorts with {D4Z4End2End}},
  journal   = {Genome Research},
  year      = {2026},
  volume    = {36},
  number    = {4},
  pages     = {827--848},
  doi       = {10.1101/gr.280907.125},
  pmid      = {41871882}
}
```

---

**[butterfield2023nanopore]** Butterfield RJ et al., 2023. Deciphering D4Z4 CpG methylation gradients in facioscapulohumeral muscular dystrophy using nanopore sequencing. *Genome Research* 33(9):1439-1454. PMID:37798116 / DOI:10.1101/gr.277871.123

Claim it supports: Report 07_integrated.md L65: "CTCF binds within D4Z4 repeat units (Ottaviani et al. 2009)" and REFERENCES_v3.bib note at dux4_d4z4_fshd_lemmers2012 (line 1137): "input for the testable D4Z4-CTCF prediction in Discussion." The draft P6 Discussion mentions that "CTCF/cohesin density at PHR boundaries should correlate with Hi-C contact strength" and that the Gershman ENCODE CTCF ChIP-seq realignment to T2T-CHM13 should be used.

Using Cas9-targeted nanopore enrichment of 4q and 10q D4Z4 arrays, this study resolves methylation at every CpG with base-pair precision. Key findings: (1) asymmetric, length-dependent methylation gradients reaching hypermethylation at ~10 D4Z4 units (the clinically established threshold); (2) a discrete region of specifically LOW methylation co-localizes with the known CTCF/insulator binding site, providing single-molecule validation that the CTCF binding site is epigenetically accessible; (3) high methylation immediately before the DUX4 transcriptional start site; (4) 180-nt periodic methylation within DUX4 exons, consistent with phased nucleosomes. This is the empirical basis for the "testable CTCF prediction" in the draft discussion.

Recommendation: ADD to support C7 (CTCF/lamin mechanism); the CTCF insulator low-methylation finding directly validates the Ottaviani2009 model at nucleotide resolution. This paper is the closest direct empirical test of the CTCF density prediction mentioned in report 07.

```bibtex
@article{butterfield2023nanopore,
  author    = {Butterfield, Russell J and Dunn, Diane M and Duval, Brett and Moldt, Sarah and Weiss, Robert B},
  title     = {Deciphering {D4Z4} {CpG} methylation gradients in facioscapulohumeral muscular dystrophy using nanopore sequencing},
  journal   = {Genome Research},
  year      = {2023},
  volume    = {33},
  number    = {9},
  pages     = {1439--1454},
  doi       = {10.1101/gr.277871.123},
  pmid      = {37798116}
}
```

---

**[lemmers2024dup]** Lemmers RJLF et al., 2024. Autosomal dominant in cis D4Z4 repeat array duplication alleles in facioscapulohumeral dystrophy. *Brain* 147(2):414-426. PMID:37703328 / DOI:10.1093/brain/awad312

Claim it supports: Draft P8 (pedigree): "The new shared sequence enters the population, strengthens sequence-similarity edges and reinforces the original 3D proximity." Report 07_integrated.md L69: "ectopic exchange between chr4_q and chr10_q D4Z4 arrays can modify FSHD alleles." C8 (ongoing exchange at 4q35 detectable at population scale).

This paper characterizes in cis D4Z4 repeat array duplications (two adjacent arrays interrupted by a spacer) in a European FSHD cohort; prevalence 1.5%. Nanopore sequencing resolves breakpoints within LSau elements flanking D4Z4 units, identifying intrachromosomal recombination as the mechanism. Specific combinations of proximal and distal array sizes determine pathogenicity; an algorithm is provided for diagnostic interpretation. DUX4 expression is confirmed from both SMCHD1-variant and SMCHD1-wild-type duplication carriers. This provides direct evidence that intrachromosomal recombination restructures the D4Z4 array at a frequency detectable in clinical cohorts (1.5%), consistent with C8's claim of ongoing exchange.

Recommendation: ADD to support C8 (structural evidence of ongoing recombination within the D4Z4-containing subtelomere). Also updates the "D4Z4 array architecture at nucleotide resolution" open question in the existing review.

```bibtex
@article{lemmers2024dup,
  author    = {Lemmers, Richard J L F and Butterfield, Russell and van der Vliet, Patrick J and de Bleecker, Jan L and van der Pol, Ludo and Dunn, Diane M and Erasmus, Corrie E and D'Hooghe, Marc and Verhoeven, Kristof and Balog, Judit and Bigot, Anne and van Engelen, Baziel and Statland, Jeffrey and Bugiardini, Enrico and van der Stoep, Nienke and Evangelista, Teresinha and Marini-Bettolo, Chiara and van den Bergh, Peter and Tawil, Rabi and Voermans, Nicol C and Vissing, John and Weiss, Robert B and van der Maarel, Silv{\`e}re M},
  title     = {Autosomal dominant in cis {D4Z4} repeat array duplication alleles in facioscapulohumeral dystrophy},
  journal   = {Brain},
  year      = {2024},
  volume    = {147},
  number    = {2},
  pages     = {414--426},
  doi       = {10.1093/brain/awad312},
  pmid      = {37703328}
}
```

---

**[sikrova2023smchd1]** Sikrova D et al., 2023. SMCHD1 and LRIF1 converge at the FSHD-associated D4Z4 repeat and LRIF1 promoter yet display different modes of action. *Communications Biology* 6(1):677. PMID:37380887 / DOI:10.1038/s42003-023-05053-0

Claim it supports: "Lemmers et al. 2012 identified SMCHD1 mutations as the cause of FSHD2, establishing a digenic model" (existing topic_06 review). Draft P1: "the D4Z4 macrosatellite of 4q35 is found in degenerate copies on 10q26" (SMCHD1 methylates both loci symmetrically in FSHD2).

SMCHD1 and LRIF1 form an auxiliary layer of D4Z4 repression. SMCHD1 directly silences LRIF1 expression by binding the LRIF1 promoter; the SMCHD1-LRIF1 interdependency differs between the D4Z4 locus and the LRIF1 promoter. Somatic loss-of-function of either protein alone does not produce the D4Z4 chromatin changes seen with germline mutations, establishing that developmental timing of SMCHD1/LRIF1 activity matters. This mechanistic detail adds to the FSHD2 digenic model (dux4_d4z4_fshd_lemmers2012) without overturning it.

Recommendation: ADD to extend the FSHD2 mechanistic discussion; particularly relevant to the lamin/nuclear-envelope story (C7) since SMCHD1 operates within the same heterochromatin maintenance pathway as the lamin A/C-D4Z4 tethering complex.

```bibtex
@article{sikrova2023smchd1,
  author    = {Sikrov{\'a}, Darina and Testa, Alessandra M and Willemsen, Iris and van den Heuvel, Anita and Tapscott, Stephen J and Daxinger, Lucia and Balog, Judit and van der Maarel, Silv{\`e}re M},
  title     = {{SMCHD1} and {LRIF1} converge at the {FSHD}-associated {D4Z4} repeat and {LRIF1} promoter yet display different modes of action},
  journal   = {Communications Biology},
  year      = {2023},
  volume    = {6},
  number    = {1},
  pages     = {677},
  doi       = {10.1038/s42003-023-05053-0},
  pmid      = {37380887}
}
```

---

**[tardy2026lrseq]** Tardy C et al., 2026. Benchmarking long-read sequencing approaches to resolve facioscapulohumeral dystrophy locus complexity. *Brain* 149(5):1750-1767. PMID:41642686 / DOI:10.1093/brain/awag029

Claim it supports: C8 (ongoing interchromosomal exchange) and the existing topic_06 review's "open questions" about structural variants; report 07_integrated.md L69 (ectopic exchange modifying FSHD alleles).

Long-read sequencing (Oxford Nanopore and PacBio) resolves seven representative cases of clinically diagnosed FSHD carrying complex structural variants at 4q35 or 10q26 not detectable by Bionano optical genome mapping. Duplicated alleles arise from intrachromosomal recombination between LSau elements within D4Z4 and distal subtelomeric beta-satellite elements, producing variable proximal D4Z4 deletions with patient-specific breakpoints. Pathogenicity of structural variants requires integration of structural and epigenetic features. This study parallels and extends lemmers2024dup (different cohort, different methods) and together they establish that structural rearrangements within the 4q35 locus are mechanistically driven by sequence-mediated intrachromosomal recombination -- the same forces that drive the interchromosomal exchange modeled in C8.

Recommendation: ADD alongside lemmers2024dup as evidence for C8.

```bibtex
@article{tardy2026lrseq,
  author    = {Tardy, Charlotte and Trani, Jean Philippe and Murcia Pienkowski, Victor and Morin, Loeva and Castro, Christel and Souville, Louis and Humbert, Camille and G{\'e}rard, Laur{\`e}ne and Eudes, Nathalie and Assoumani, Amire and Bertaux, Karine and Verebi, Camille and Nectoux, Juliette and Salort Campana, Emmanuelle and Jacquemont, Marie Line and Toutain, Annick and Mallaret, Martial and Tard, C{\'e}line and Fradin, M{\'e}lanie and Attarian, Shahram and Nguyen, Karine and Bernard, Rafaelle and Magdinier, Fr{\'e}d{\'e}rique},
  title     = {Benchmarking long-read sequencing approaches to resolve facioscapulohumeral dystrophy locus complexity},
  journal   = {Brain},
  year      = {2026},
  volume    = {149},
  number    = {5},
  pages     = {1750--1767},
  doi       = {10.1093/brain/awag029},
  pmid      = {41642686}
}
```

---

**[salsi2025genotoxic]** Salsi V et al., 2025. Posttranscriptional RNA stabilization of telomeric RNAs FRG2, DBE-T, D4Z4 at human 4q35 in response to genotoxic stress and D4Z4 macrosatellite repeat length. *Clinical Epigenetics* 17(1):73. PMID:40320530 / DOI:10.1186/s13148-025-01881-5

Claim it supports: Cabianca et al. 2012 (dux4_d4z4_fshd_cabianca2012): "DBE-T recruits Trithorax-group protein Ash1L to the 4q35 locus." Draft P6: "The D4Z4 array on 4q35 binds CTCF and is tethered to the nuclear lamina." Report 03_gene_enrichment.md: C1 annotation includes "DBET."

The 4q35 subtelomere is subdivided into discrete chromatin domains: centromeric genes (SLC25A4, FAT1, FRG1) carry active histone marks, while telomeric loci (FRG2, DBE-T, D4Z4) carry poised or repressed marks. Under DNA damage, the telomeric RNAs FRG2, DBE-T, and D4Z4-derived transcripts are induced and stabilized posttranscriptionally in a D4Z4-copy-number-dependent manner (more induction with fewer repeats). This adds a genotoxic-stress-dependent RNA stabilization layer to the Cabianca2012 Polycomb/Trithorax switch and provides new evidence for the discrete chromatin domain organization of 4q35 that underlies the CTCF insulator model (C7).

Recommendation: ADD to extend Cabianca2012 and support the 4q35 domain architecture in C7.

```bibtex
@article{salsi2025genotoxic,
  author    = {Salsi, Valentina and Losi, Francesca and Salani, Monica and Kaufman, Paul D and Tupler, Rossella},
  title     = {Posttranscriptional {RNA} stabilization of telomeric {RNAs} {FRG2}, {DBE-T}, {D4Z4} at human 4q35 in response to genotoxic stress and {D4Z4} macrosatellite repeat length},
  journal   = {Clinical Epigenetics},
  year      = {2025},
  volume    = {17},
  number    = {1},
  pages     = {73},
  doi       = {10.1186/s13148-025-01881-5},
  pmid      = {40320530}
}
```

---

### MEDIUM relevance papers (do not add to REFERENCES_v4.bib; background context only)

- **Gao et al. 2026** (PMID 41832971, DOI 10.1093/procel/pwag014): DUX/DUX4 phase-separated condensates recruit CTCF to MERVL/MT2 super-enhancers; 3D genome reorganization. MEDIUM: adds mechanistic detail to DUX4-CTCF link (C7) but is an embryology/totipotency study not directly cited in the draft.
- **Hamm et al. 2026** (PMID 41986233, DOI 10.1101/gad.353253.125): KLF18 (VNTR-containing) is a DUX4 feed-forward component in ZGA. MEDIUM: extends dux4_d4z4_fshd_hendrickson2017 but not directly relevant to C5/C7/C8.
- **Fox et al. 2024** (PMID 39627769, DOI 10.1186/s13395-024-00361-3): SIX1/SIX2/SIX4 transcription factors required for DUX4 expression in differentiating FSHD myotubes. MEDIUM: mechanistic but not tied to the structural/genomic claims.
- **Belayew et al. 2025** (PMID 40855454, DOI 10.1186/s13395-025-00388-0): 25-year DUX4 review. MEDIUM: useful general reference but not claim-specific.
- **Smith et al. 2023** (PMID 37691147, DOI 10.1016/j.celrep.2023.113114): DUX4 in cancer induces metastable ZGA program + MHC suppression. MEDIUM: extends ZGA biology but not draft-claim-specific.

---

## Section 4: CONTRADICTIONS

**Partial contradiction (salsi2026t2t):** Salsi et al. 2026 (PMID 41535478) contradicts the implicit claim that D4Z4 repeats are "predominantly confined to 4q35 and 10q26 loci." Using T2T-CHM13, they find D4Z4-like loci with intact DUX4 ORFs or polyadenylation signals on at least ten additional chromosomes. The FSHD literature and our draft treat the 4q/10q pair as the only pathologically relevant D4Z4 loci; this assumption is operationally valid (no patient series with 10+ chromosomal D4Z4 disease origin), but the T2T genome reveals that "degenerate copies" are not confined to 10q26. The impact on C5 is limited: community C1 clustering reflects the dominant 4q/10q arms whose cross-arm sequence rate drives the NJ clade. However, the primer contamination finding (PCR assays for DUX4/DBE-T amplify off-target loci) is a methodological concern that affects interpretation of some expression data cited in the existing review (Snider2009, Geng2012, Cabianca2012 if those studies used standard short-read or PCR approaches). The draft should acknowledge that the T2T era requires locus-resolved, repeat-aware approaches for all DUX4 expression and methylation studies.

No contradictions to C7 (lamin/CTCF mechanism) or C8 (ongoing exchange) were found in the searched window.

---

## Section 5: Search audit trail

**Tools used:**
1. `mcp__claude_ai_PubMed__search_articles` — primary discovery
2. `mcp__claude_ai_PubMed__get_article_metadata` — metadata verification
3. `mcp__claude_ai_bioRxiv__search_preprints` (genetics category, 2023-01-01 to 2026-05-17) — no DUX4/FSHD preprints returned; bioRxiv MCP does not support keyword search, only date+category filter; category "genetics" returns general genetics preprints, not FSHD-specific hits

**Date range filter applied:** 2023/01/01 — 2026/12/31 (PubMed `pdat` filter)

**PubMed queries and hit counts:**

| Query | Hits before relevance filter | Hits after filter |
|-------|------------------------------|-------------------|
| `DUX4 FSHD facioscapulohumeral dystrophy` | 113 | 30 retrieved; 7 STRONG, 5 MEDIUM, 18 WEAK |
| `D4Z4 macrosatellite repeat FSHD` | 19 | 19 retrieved; overlaps with above; 3 new STRONG |
| `DUX4 zygotic genome activation totipotency embryo` | 7 | 7 retrieved; 2 MEDIUM, 5 WEAK |
| `SMCHD1 FSHD2 digenic D4Z4 hypomethylation` | 0 | — (query too restrictive; SMCHD1 papers captured via first query) |
| `D4Z4 4q35 subtelomere CTCF lamin nuclear lamina chromatin` | 0 | — (query too restrictive) |
| `Lemmers FSHD subtelomere 4q 10q haplotype population` | 0 | — (no new Lemmers population paper in window; Lemmers2024dup found via first query) |

**Total unique PMIDs retrieved:** ~43 across all queries  
**After relevance filter (STRONG + MEDIUM):** 12  
**Added to Section 3 (STRONG):** 7

**Spot-checks performed (3 required):**
- PMID 41535478: confirmed via `get_article_metadata` — Salsi et al. 2026, Eur J Hum Genet 34(3):357-367, DOI 10.1038/s41431-025-02000-x. VERIFIED.
- PMID 37798116: confirmed via `get_article_metadata` — Butterfield et al. 2023, Genome Research 33(9):1439-1454, DOI 10.1101/gr.277871.123. VERIFIED.
- PMID 37703328: confirmed via `get_article_metadata` — Lemmers et al. 2024, Brain 147(2):414-426, DOI 10.1093/brain/awad312. VERIFIED.

**Why specific hits were dropped:**

| PMID | Title (abbreviated) | Reason dropped |
|------|---------------------|----------------|
| 42066431 | Spatiotemporal DUX4 toxicity — biomarker-driven therapies | Clinical/therapeutic review; no structural genomic claim relevance |
| 41994867 | AOC 1020 antibody-oligonucleotide conjugate for FSHD | Therapeutic tool; irrelevant to C5/C7/C8 |
| 41944163 | Transgenic mouse models for DUX4 | Mouse model tool; WEAK for structural claims |
| 41781309 | Clinical overview FSHD features | Clinical review; no structural claim |
| 41889948 | AI drug screen BAZ1A bromodomain | Therapeutic; off-topic |
| 41536811 | ASO gapmer subcutaneous delivery | Therapeutic; off-topic |
| 41510809 | KHDC1L plasma biomarker | Biomarker; off-topic |
| 41652446 | DUCKS4 Nanopore workflow | Diagnostic tool; methods paper; not claim-specific |
| 41649965 | Losmapimod phase 3 REACH study | Clinical trial (failed); off-topic |
| 41329166 | DUX4-induced HSATII RNA aggregation | Interesting mechanism but no direct C5/C7/C8 link |
| 41540238 | DUX4/DUX4C integrated genetic approach T2T | MEDIUM; DUX4C haplotypes and LD with DUX4; subsumed by salsi2026t2t for T2T angle |
| 39603552 | OGM maternal mosaicism FSHD1 siblings | Diagnostic/case report; no structural claim |
| 37565369 | Methylation protocol optimization D4Z4 DR1 | Methods; no new claim |
| 38002249 | OGM for FSHD molecular diagnosis | Diagnostic validation; no new claim |
| 40855454 | DUX4 at 25 (review) | General review; MEDIUM, not claim-specific enough for Section 3 |
| 40226918 | ZBTB24-CDCA7-HELLS axis suppresses 2C-like reprogramming | Mouse mESC/Dux; not human 4q35 directly |
| 38915486 | PARP-DUX4 regulatory axis TIRN stem cells | Preprint; broad ZGA biology, not 4q35-specific |
| 37356343 | Human 8-cell embryo-like cells review | General ZGA review; WEAK |
| 41513943 | SUMOylation inhibition activates FSHD locus | Mechanism (SUMOylation independent of methylation/SMCHD1); interesting but not C5/C7/C8 |
| 41756954 | SMCHD1 loss rewires MYOD1 enhancer nexuses | MEDIUM; SMCHD1 3D chromatin; addressed by sikrova2023smchd1 which is more directly relevant |
| 37691147 | DUX4 in cancer, metastable ZGA program | MEDIUM; ZGA biology not structural claim |
