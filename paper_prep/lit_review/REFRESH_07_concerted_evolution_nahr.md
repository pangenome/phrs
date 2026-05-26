# REFRESH 07: Concerted Evolution, NAHR, Gene Conversion

**Date:** 2026-05-17
**Task:** lit-refresh-07-concerted
**Agent:** agent-79

---

## Section 1: Topic scope

Topic 07 supplies the mechanistic vocabulary — molecular drive (Dover), non-allelic homologous recombination (NAHR), gene conversion, break-induced replication (BIR), and GC-biased gene conversion (gBGC) — that underpins the abstract's central claim: human subtelomeres are maintained at high inter-chromosomal identity by ongoing recombination exchange, i.e., concerted evolution in the loose sense of Arnheim, Ohta, and Charlesworth. This topic bears directly on **draft paragraph P3** (the neighbour-joining tree of six clades, each of which is predicted to be maintained by NAHR/gene conversion), **draft paragraph P5** (pedigree resolution of 538 patches: 133 gene-conversion-like, 16 crossover-like, 229 acros_like — all classified according to the NAHR/gene-conversion mechanistic vocabulary), and **draft paragraph P6** (the causal-loop close, where "concerted evolution" is explicitly attributed to Arnheim1980, Ohta1984, and Charlesworth1994, and the limitations note that mechanism discrimination between NAHR, gene conversion, and BIR remains open). The open question of whether BIR tracts of 100+ kb are miscategorised as gene conversion (raised in Anand2013) directly concerns the sandwich-pattern calls in section 14 of the end-to-end report.

---

## Section 2: Existing citations still authoritative

All bibkeys from `topic_07_concerted_evolution_nahr.bib` are assessed below. The prior universe is the 22 entries in that file plus the 9 `concerted_evolution_nahr_`-prefixed keys in `REFERENCES_v3.bib` that the draft already cites.

| BibKey | Status | Notes |
|---|---|---|
| Smith1976 | STILL CURRENT | Foundational theory; no successor |
| Brown1972 | STILL CURRENT | Founding empirical case (Xenopus rDNA); no successor |
| Dover1982 | STILL CURRENT | Named molecular drive; no successor |
| Dover1986 | STILL CURRENT | Molecular drive synthesis; no successor |
| Liao1999 | STILL CURRENT | Field synthesis; no replacement review at this level |
| Lupski1998 | STILL CURRENT | Foundational genomic-disorders paper |
| StankiewiczLupski2002 | STILL CURRENT | Canonical NAHR architecture paper; still primary citation |
| LupskiStankiewicz2005 | STILL CURRENT | Mechanisms review; no new synthesis at this scope |
| Lupski2009 | STILL CURRENT | 10-year review; no newer equivalent general review |
| BaileyEichler2006 | STILL CURRENT | Primate SDs pericentromeric/subtelomeric bias; still primary citation |
| ChenCooper2007 | STILL CURRENT | Gene-conversion review of record; no 2023+ replacement review found |
| Hastings2009 | STILL CURRENT | CNV mechanisms review (NRG); still the standard reference |
| EickbushEickbush2007 | STILL CURRENT | rDNA homogenization case study; no successor |
| GanleyKobayashi2007 | STILL CURRENT | Whole-genome rDNA intra-array quantification; no successor |
| NeiRooney2005 | STILL CURRENT | Birth-and-death alternative; foundational |
| Galtier2001 | STILL CURRENT | gBGC isochore hypothesis; foundational |
| DuretGaltier2009 | STILL CURRENT | gBGC mechanisms review; supplemented but not superseded by Clessin2025 |
| JeffreysMay2004 | STILL CURRENT | Gene-conversion tract length at crossover hotspots; unique measurement |
| Williams2015 | STILL CURRENT | Non-crossover gene conversions with GC bias; no closer replacement |
| Halldorsson2019 | STILL CURRENT | Sequence-level genetic map; closest pedigree-resolution external comparator |
| Llorente2008 | STILL CURRENT | BIR mechanism review; no replacement |
| Anand2013 | STILL CURRENT | BIR review (CSHPB); still standard BIR reference |
| Guarracino2023 | STILL CURRENT | Acrocentric PHR recombination (Nature 2023); the direct precursor paper |
| Logsdon2024 | STILL CURRENT | T2T centromere/subtelomere completion; architectural complement |
| concerted_evolution_nahr_Arnheim1980 | STILL CURRENT | Foundational loose-sense concerted evolution; explicitly named in draft P6 |
| concerted_evolution_nahr_Charlesworth1994 | STILL CURRENT | Repetitive DNA dynamics; explicitly named in draft P6 |
| concerted_evolution_nahr_Eichler2001 | STILL CURRENT | SDs and genomic instability (NRG); still primary citation |
| concerted_evolution_nahr_Hastings2009 | STILL CURRENT | MMBIR model (PLoS Genet); still standard MMBIR reference |
| concerted_evolution_nahr_Hillis1991 | STILL CURRENT | gBGC in concerted evolution; foundational |
| concerted_evolution_nahr_Myers2010 | STILL CURRENT | PRDM9 and hotspot erosion; still primary PRDM9 citation |
| concerted_evolution_nahr_Ohta1984 | STILL CURRENT | Multigene family gene conversion models; explicitly named in draft P6 |
| concerted_evolution_nahr_SamonteEichler2002 | STILL CURRENT | SD evolution of primate genome (NRG); still primary SD evolution citation |
| concerted_evolution_nahr_Vollger2023 | STILL CURRENT | Increased mutation and gene conversion in SDs (Nature 2023); primary 2023 SD reference |

No existing citation is superseded. The field is classically oriented — foundational papers remain load-bearing, and new papers extend rather than replace them.

---

## Section 3: NEW papers to add (STRONG relevance, 2023+)

---

**[noyes2026sd]** Noyes MD, Sui Y, Kwon Y et al., 2026. Long-read sequencing of families reveals increased germline and postzygotic mutation rates in repetitive DNA. *Nature Communications* 17(1). PMID:41803180 / DOI:[10.1038/s41467-026-70342-1](https://doi.org/10.1038/s41467-026-70342-1)

Claim it supports: Draft paragraph P5: "133 are gene-conversion-like sandwiches … the central chrM:hY block aligning at 1.000 identity"; report section 14 (interlocus gene conversion as the dominant mechanism). Also supports the open question in the existing review (topic_07 §"Open questions"): "the relative contribution of NAHR vs gene conversion vs BIR … has not been resolved at HPRC scale."

Applied Illumina + ONT + PacBio long-read sequencing to 73 children from 42 autism families (157 individuals). Germline mutation rate 1.30 × 10⁻⁸ substitutions/bp/generation; postzygotic rate 0.23 × 10⁻⁸. Mutation enrichment in segmental duplications is strongly postzygotic (not germline), driven by faulty DNA repair and **interlocus gene conversion**. SD mutability is length- and percent-identity-dependent. Directly demonstrates that interlocus gene conversion between highly similar SDs operates at measurable rates in humans, extending Vollger2023 (already cited) to the postzygotic and germline decomposition level and providing the mechanism-specificity that Anand2013 flagged as unresolved.

Recommendation: **ADD** alongside `concerted_evolution_nahr_Vollger2023` and `Anand2013` as mechanistic support for the gene-conversion-like patch interpretation.

```bibtex
@article{noyes2026sd,
  author       = {Noyes, Michelle D. and Sui, Yang and Kwon, Youngjun and Koundinya, Nidhi and Wong, Isaac and Munson, Katherine M. and Hoekzema, Kendra and Kordosky, Jennifer and Garcia, Gage H. and Knuth, Jordan and Lewis, Alexandra P. and Eichler, Evan E.},
  title        = {Long-read sequencing of families reveals increased germline and postzygotic mutation rates in repetitive {DNA}},
  journal      = {Nature Communications},
  year         = {2026},
  volume       = {17},
  number       = {1},
  doi          = {10.1038/s41467-026-70342-1},
  pmid         = {41803180},
  note         = {Long-read + short-read trio sequencing of 42 autism families; interlocus gene conversion drives postzygotic SD mutability; directly supports the gene-conversion-like patch mechanism in PHR pedigree analysis}
}
```

---

**[chen2025paraphase]** Chen X, Baker D, Dolzhenko E et al., 2025. Genome-wide profiling of highly similar paralogous genes using HiFi sequencing. *Nature Communications* 16(1):2340. PMID:40057485 / DOI:[10.1038/s41467-025-57505-2](https://doi.org/10.1038/s41467-025-57505-2)

Claim it supports: Draft paragraph P3: "the named clades fit NAHR predictions in detail … extensive gene conversion and unequal crossing over contribute to highly similar gene copies" (existing topic_07 review §C5). Specifically, the 23 paralog groups with exceptionally low within-group diversity driven by gene conversion and unequal crossing-over match the prediction in the draft that sequence-level communities correspond to NAHR/gene-conversion-maintained paralog clusters.

Applied the Paraphase method (HiFi-based haplotype phasing) to 160 long (>10 kb) segmental duplication regions encoding 316 genes across five ancestral populations. Identified 23 paralog groups with exceptionally low within-group diversity attributable to extensive gene conversion and unequal crossing-over. In 36 trios, found 7 de novo SNVs and **4 de novo gene conversion events, 2 of which are non-allelic**. Non-allelic de novo gene conversions are direct in vivo observations of the mechanism the draft invokes for inter-chromosomal exchange. Population-variable copy numbers reinforce C5 (community structure from sequence identity).

Recommendation: **ADD** alongside `ChenCooper2007` and `concerted_evolution_nahr_Vollger2023` for de novo gene conversion evidence and population-scale paralogy maintenance.

```bibtex
@article{chen2025paraphase,
  author       = {Chen, Xiao and Baker, Daniel and Dolzhenko, Egor and Devaney, Joseph M. and Noya, Jessica and Berlyoung, April S. and Brandon, Rhonda and Hruska, Kathleen S. and Lochovsky, Lucas and Kruszka, Paul and Newman, Scott and Farrow, Emily and Thiffault, Isabelle and Pastinen, Tomi and Kasperaviciute, Dalia and Gilissen, Christian and Vissers, Lisenka and Hoischen, Alexander and Berger, Seth and Vilain, Eric and D{\'e}lot, Emmanu{\`e}le and Eberle, Michael A.},
  title        = {Genome-wide profiling of highly similar paralogous genes using {HiFi} sequencing},
  journal      = {Nature Communications},
  year         = {2025},
  volume       = {16},
  number       = {1},
  pages        = {2340},
  doi          = {10.1038/s41467-025-57505-2},
  pmid         = {40057485},
  note         = {Paraphase method; 23 paralog groups with low within-group diversity from gene conversion; 4 de novo gene conversion events in trios, 2 non-allelic; directly demonstrates ongoing NAHR/gene conversion between SD paralogs}
}
```

---

**[hinch2023meiotic]** Hinch R, Donnelly P, Hinch AG, 2023. Meiotic DNA breaks drive multifaceted mutagenesis in the human germ line. *Science* 382(6674):eadh2531. PMID:38033082 / DOI:[10.1126/science.adh2531](https://doi.org/10.1126/science.adh2531)

Claim it supports: Draft paragraph P5 / report section 14: the pedigree calls classify patches as "gene-conversion-like" and "crossover-like" based on the assumption that meiotic DSBs at these loci produce repair events. The Hinch2023 paper directly quantifies how mutagenic meiotic break repair is and what repair mechanisms are invoked. Meiotic break repair is 8-fold more mutagenic for single-base substitutions than previously understood (de novo mutation in one in four sperm, one in 12 eggs); impact on indels and structural variants is 100- to 1300-fold higher per break. Error-prone repair mechanisms include translesion synthesis and end joining. This extends the mechanistic repertoire beyond the canonical NAHR/SDSA/BIR triad (Anand2013, Hastings2009) and introduces translesion synthesis as a contributor to mutagenesis at meiotic DSB sites — a mechanism not previously linked to the subtelomere/PHR context but now worth considering for patch interpretation.

Recommendation: **ADD** as background support for the open-question discussion of BIR vs gene conversion miscategorisation (already flagged in topic_07 §"Open questions"), citing it alongside Anand2013 and Llorente2008.

```bibtex
@article{hinch2023meiotic,
  author       = {Hinch, Robert and Donnelly, Peter and Hinch, Anjali Gupta},
  title        = {Meiotic {DNA} breaks drive multifaceted mutagenesis in the human germ line},
  journal      = {Science},
  year         = {2023},
  volume       = {382},
  number       = {6674},
  pages        = {eadh2531},
  doi          = {10.1126/science.adh2531},
  pmid         = {38033082},
  note         = {Meiotic break repair 8x more mutagenic than previously known; error-prone repair (TLS, end joining) at DSB sites; extends BIR/SDSA mechanistic framework relevant to PHR patch interpretation}
}
```

---

**[clessin2025gbgc]** Clessin A, Joseph J, Lartillot N, 2025. Evolution of GC-biased gene conversion by natural selection. *Genetics* 230(4). PMID:40510024 / DOI:[10.1093/genetics/iyaf111](https://doi.org/10.1093/genetics/iyaf111)

Claim it supports: The topic_07 review §"Open questions": "whether gBGC leaves a measurable fingerprint on subtelomeric base composition is testable." This paper models gBGC as a quantitative trait under stabilizing selection: gBGC is maintained at a positive value (favoured by selection despite the genetic load it imposes in high-recombination regions). The implication for the draft is that GC-rich block structure at subtelomeric PHR communities may be both a consequence of gBGC and a stable equilibrium, not a transient state — strengthening the claim that sequence-level communities retain GC-rich block structure (topic_07 §C5, last sentence about DuretGaltier2009). Does not supersede DuretGaltier2009 but adds a fitness/population-genetic mechanism layer.

Recommendation: **ADD** alongside DuretGaltier2009 in the gBGC discussion thread.

```bibtex
@article{clessin2025gbgc,
  author       = {Clessin, Augustin and Joseph, Julien and Lartillot, Nicolas},
  title        = {Evolution of {GC}-biased gene conversion by natural selection},
  journal      = {Genetics},
  year         = {2025},
  volume       = {230},
  number       = {4},
  doi          = {10.1093/genetics/iyaf111},
  pmid         = {40510024},
  note         = {gBGC under stabilizing selection; maintained positive by selection despite load; strengthens inference that GC-rich subtelomeric block structure reflects stable gBGC equilibrium}
}
```

---

**[hebbar2026marmoset]** Hebbar P, Potapova T, Loucks H et al., 2026. A Complete Genome for the Common Marmoset. *bioRxiv* (preprint). PMID:41929024 / DOI:[10.64898/2026.03.25.713844](https://doi.org/10.64898/2026.03.25.713844)

Claim it supports: Draft paragraph P5 / report section 12: the broader claim that ongoing recombinational exchange between heterologous chromosomes at PHRs is a conserved feature of primate genomes, not a human-specific artifact. T2T marmoset assembly shows acrocentric autosomes with PHRs; chromosomes sharing PHRs also share closely related centromeric satellite DNA, supporting ongoing recombinational exchange between heterologous chromosomes. The rDNA and PHR recombination pattern is conserved across primates. Directly extends Guarracino2023 beyond human/great apes to New World monkeys.

**Note:** This is a preprint (bioRxiv). The DOI pattern (10.64898/) is non-standard; verify final published version before inclusion in REFERENCES_v4.

Recommendation: **ADD** as MEDIUM evidence for primate conservation of the PHR recombination mechanism, clearly flagged as preprint. Pair with Guarracino2023.

```bibtex
@misc{hebbar2026marmoset,
  author       = {Hebbar, Prajna and Potapova, Tamara and Loucks, Hailey and Ray, Karina and Rodrigues, Murillo F. and Ryabov, Fedor and Malukiewicz, Joanna and Yoo, DongAhn and de Lima, Leonardo and Haber, Annat and others and Gerton, Jennifer L. and Alexandrov, Ivan and Paten, Benedict},
  title        = {A Complete Genome for the Common Marmoset},
  howpublished = {bioRxiv preprint},
  year         = {2026},
  doi          = {10.64898/2026.03.25.713844},
  pmid         = {41929024},
  note         = {PREPRINT. T2T marmoset genome; acrocentric PHRs shared with centromeric satellite DNA; supports model of ongoing recombinational exchange between heterologous chromosomes as conserved primate feature}
}
```

---

## Section 4: CONTRADICTIONS

No papers found in the searched 2023–2026 window that directly contradict claims made in the draft or in the existing topic_07 review.

Relevant null findings:
- No paper challenges the existence of ongoing inter-chromosomal gene conversion at PHRs or disputes the acrocentric recombination evidence (Guarracino2023).
- No paper proposes that subtelomeric sequence similarity is maintained by mechanisms other than recombination-based exchange (the birth-and-death alternative of NeiRooney2005 remains a caveat for specific gene families, not a general contradiction).
- The Noyes2026 finding that repeat-region mutations are predominantly postzygotic (not germline) does not contradict the draft's pedigree claims, which report both meiotic (PAN028 12-generation) and non-meiotic patterns.
- Clessin2025 on gBGC under selection refines rather than contradicts the DuretGaltier2009 framework; no claim in the draft depends on whether gBGC is neutral or selected.

---

## Section 5: Search audit trail

### Tools used
1. **mcp__claude_ai_PubMed__search_articles** — primary discovery tool
2. **mcp__claude_ai_PubMed__get_article_metadata** — metadata verification for all candidate PMIDs
3. **OpenAlex API** (via openalex-database skill, Python client) — cross-validation and catch additional papers
4. **mcp__claude_ai_bioRxiv__search_preprints** — attempted; bioRxiv API does not support keyword search, only date+category; returned general genomics papers, not useful for this topic

### Date range filter
2023-01-01 to 2026-05-17

### Query strings (PubMed)

| Query | Hits before filter | Hits after filter |
|---|---|---|
| `concerted evolution NAHR non-allelic homologous recombination gene conversion` | 0 | 0 |
| `gene conversion segmental duplications human genome[Title/Abstract]` | 1 | 1 (PMID 40057485, STRONG) |
| `non-allelic homologous recombination NAHR genomic rearrangements[Title/Abstract]` | 1 | 1 (PMID 37775806, WEAK — computational algorithm) |
| `NAHR segmental duplication rearrangement genomic disorder` | 1 | 1 (PMID 38576798, WEAK — case report) |
| `gene conversion recombination human pedigree sequence resolution` | 0 | 0 |
| `break-induced replication copy number variation genome instability[Title/Abstract]` | 0 | 0 |
| `biased gene conversion GC recombination mammalian genome[Title/Abstract]` | 0 | 0 |
| `biased gene conversion gBGC GC-content human genome recombination` | 1 | 1 (PMID 40510024, MEDIUM) |
| `PRDM9 recombination hotspot meiosis evolution` | 3 | 0 STRONG (all WEAK for this topic) |
| `acrocentric chromosome recombination heterologous rDNA T2T` | 2 | 1 STRONG (PMID 37165241 = Guarracino2023, already in bib); 1 MEDIUM (PMID 41929024 marmoset preprint) |
| `segmental duplication mutation rate interlocus gene conversion pangenome` | 2 | 2 (PMID 41803180, 40791370 — same paper published + preprint; STRONG) |
| `break-induced replication mechanism template switching genome instability` | 3 | 0 STRONG (53BP1 cancer context WEAK; TANGO2 WEAK) |
| `rDNA ribosomal repeat homogenization copy number concerted evolution pangenome` | 2 | 0 STRONG (duckweed WEAK; Daphnia WEAK) |
| `Hinch meiotic DNA breaks mutagenesis human germ line recombination` | 1 | 1 (PMID 38033082, MEDIUM-STRONG) |

### OpenAlex query strings

| Query | Hits before filter | Hits after filter |
|---|---|---|
| `concerted evolution NAHR gene conversion subtelomere pangenome` | 0 | 0 |
| `interlocus gene conversion segmental duplication human long-read pedigree` | 6 | 2 novel (Noyes2026 confirmed; Logsdon2024 complex variants WEAK for this topic) |
| `break-induced replication genome instability copy number variation` | 10 | 0 STRONG (off-topic papers dominating) |
| `GC-biased gene conversion recombination mammalian genome evolution` | 10 | 0 novel (Vollger2023 already in bib) |
| `non-allelic homologous recombination NAHR segmental duplication genomic disorder` | 10 | 0 novel |

### Dropped papers (1-line reason each)

| PMID / DOI | Title snippet | Reason dropped |
|---|---|---|
| PMID 37775806 | Constructing founder sets under allelic and NAHR | Computational algorithm paper; no biological insight on mechanisms |
| PMID 38576798 | Detection of 4p16.3 deletion… case report | Clinical case report; adds no mechanistic or population-level novelty |
| PMID 38768268 | PRDM9 intra-genomic Red Queen model | Theoretical model of hotspot evolution; Myers2010 already cited; tangential to subtelomere claims |
| PMID 37830496 | Down the Penrose stairs (PRDM9) | Alternative PRDM9 model; Myers2010 already covers hotspot context; tangential |
| PMID 39761307 | PRDM9 drives recombination hotspots in salmonid fish | Non-human vertebrate; too far from the draft's focus |
| PMID 39368985 | 53BP1 deficiency leads to hyperrecombination using BIR | BIR in cancer/DSB repair regulation context; does not address subtelomeric gene conversion |
| PMID 42063611 | Long-read HiFi resolves retrotransposon-mediated TANGO2 deletions | FoSTeS/MMBIR at a specific disease locus; no population-scale relevance |
| PMID 38711607 | 5S rDNA repeats in Spirodela polyrhiza (duckweed) | Plant model; no relevance to human subtelomere/acrocentric topic |
| PMID 41261963 | Copy number variation in Daphnia rDNA | Invertebrate model; recombination rate estimation tangential |
| DOI 10.1101/gr.277334.122 | Gaps and complex structurally variant loci (Porubsky) | SV calling in phased assemblies; methods paper, no new mechanistic claims on gene conversion |
| DOI 10.1101/2024.09.24.614721 | Complex genetic variation in nearly complete human genomes (Logsdon) | SV detection; no new concerted evolution/NAHR mechanistic content |
| DOI 10.1101/2024.09.26.615256 | Human-specific gene expansions (Soto) | Gene family evolution via SDs; birth-and-death context; no new NAHR/gene conversion mechanism |
