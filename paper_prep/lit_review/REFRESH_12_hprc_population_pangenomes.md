# REFRESH 12: HPRC Population Pangenomes — Literature Update (2023–2026)

Generated: 2026-05-17
Task: lit-refresh-12-hprc

---

## Section 1: Topic Scope

Topic 12 covers the Human Pangenome Reference Consortium (HPRC) v1–v2 lineage as the enabling substrate for the subtelomere analysis. In OUR materials it supports claim **C1** (the subtelomere paper is a companion to the HPRC v2 main manuscript) and **C3** (465 near-complete assemblies as the analysis substrate), and provides motivating background for **C4** (genome-wide identity survey) and **C6** (population-genetic analyses). The topic traces the intellectual arc from the 1000 Genomes Project short-read cohorts through HGSVC long-read pilots, T2T-CHM13, HPRC v1 (Liao et al. 2023), and finally HPRC v2 (data release May 2025, manuscript in submission). In the Nature draft it appears explicitly in main-text P1 (introductory paragraph citing `[@hprc_hprcv2_2025; @Liao2023]` and invoking the 232-individual / 465-assembly substrate) and P2 (Methods paragraph describing sample selection and assembly provenance). In the end-to-end report it underpins `01_pipeline.md` (contig classification and flank extraction) and `12_literature.md` (confirmed-literature claims map). In the slides deck `zoom_review_deck.typ` the HPRC v2 data release is the first substantive scientific slide after the title.

---

## Section 2: Existing Citations — Current or Superseded

| Bibkey | Status | Reason |
|--------|--------|--------|
| `hprc_1000g2012` | STILL CURRENT | Foundational population cohort; no successor replaces its historical role. |
| `hprc_1000g2015` | STILL CURRENT | Phase 3 reference for short-read population SV landscape; still the benchmark for large-cohort short-read approaches. |
| `hprc_sudmant2015` | STILL CURRENT | First population-stratified SV map; the 68,818-deletion count is still cited as the short-read ceiling. |
| `hprc_audano2019` | STILL CURRENT | 9× subtelomeric SV enrichment finding remains the key empirical anchor for why subtelomeres require long-read assembly. |
| `hprc_collins2020` | STILL CURRENT | gnomAD-SV remains the largest short-read SV resource; its ceiling still motivates the long-read shift. |
| `hprc_ebert2021` | STILL CURRENT | HGSVC2 (64 haplotypes); now directly preceded by HGSVC3 (Logsdon et al. 2025, see Section 3 #1). Still authoritative as the step between HGSVC pilot and HPRC. |
| `hprc_cheng2021` | STILL CURRENT | hifiasm is the primary HPRC v2 assembler; the 2021 paper is still the canonical citation. |
| `hprc_li2020` | STILL CURRENT | Minigraph underpins Minigraph-Cactus; the 2020 paper is the foundational tool citation. |
| `hprc_nurk2022` | STILL CURRENT | T2T-CHM13 is the reference anchor in the present manuscript (`CHM13v2.0`); not superseded. |
| `hprc_wang2022` | STILL CURRENT | HPRC vision paper sets the 350-individual goal; still the institutional framing reference. |
| `hprc_rautiainen2023` | STILL CURRENT | Verkko 2023 is the co-assembler for HPRC v2 complex loci; still the primary Verkko citation. |
| `hprc_liao2023` | STILL CURRENT | HPRC v1 main paper, 1,113 citations as of 2026-05; the definitive HPRC v1 reference. |
| `hprc_hickey2024` | STILL CURRENT | Minigraph-Cactus is the HPRC v1 primary graph; still authoritative. |
| `hprc_garrison2024` | STILL CURRENT | pggb/odgi is the PHR graph construction method; directly used in this manuscript. |
| `hprc_ebler2022` | STILL CURRENT | PanGenie k-mer genotyping; still the state-of-the-art short-read → pangenome genotyping tool. |
| `hprc_siren2025` | STILL CURRENT | Published Nature 644:442–452 (2025); the 1,019-genome long-read SV catalogue directly cited in the "Open questions" section. |
| `hprc_hprcv2_2025` | **UPDATE NEEDED** | Currently `@misc` (pre-publication data release, "manuscript in submission" as of May 2025). As of May 2026 the manuscript has not yet indexed in PubMed under its expected DOI. Recommend checking humanpangenome.org and bioRxiv for a preprint or formal publication DOI; upgrade from `@misc` to `@article` when available. The existing URL entry is still citable as a data release, but its `note` field must be updated if a formal paper has appeared. |

---

## Section 3: NEW Papers to Add (STRONG relevance, 2023+)

---

**[logsdon2025hgsvc]** Logsdon GA, Ebert P, Audano PA, Loftus M, Porubsky D, Ebler J, et al. 2025. Complex genetic variation in nearly complete human genomes. *Nature*. PMID:40702183 / DOI:[10.1038/s41586-025-09140-6](https://doi.org/10.1038/s41586-025-09140-6)

Claim it supports/contradicts:
- **Supports** draft P1/P2: "The Human Pangenome Reference Consortium v2 release delivers near-T2T haplotype-resolved assemblies for 232 individuals…" (main text) — the HGSVC3 paper is the companion SV analysis at 65 diverse genomes that demonstrates that combining HPRC v1 pangenome with T2T-quality assemblies enables detection of 26,115 SVs/individual at median QV 45 from short reads, a 2–3× increase over linear-reference calling. It directly extends the HGSVC2 (hprc_ebert2021) entry in the existing review's chronology.
- **Supports** draft Methods §"Sample selection": median 130 Mb assembly continuity, 39% of chromosomes at T2T status across 65 genomes, sets the quality expectation for HPRC v2 assemblies.
- **Strengthens** the claim that pangenome-based genotyping (hprc_ebler2022) now reaches biobank scale: the paper demonstrates this on 1,852 resolved complex SVs using the HPRC pangenome.

Recommendation: **ADD** as a new entry for the HGSVC3 lineage. Does not replace hprc_ebert2021 (predecessor) but adds the current endpoint of the HGSVC series. Cite in the "Era overview" narrative of topic 12 between hprc_ebert2021 and hprc_liao2023.

Bibtex stanza:
```bibtex
@article{logsdon2025hgsvc,
  author  = {Logsdon, Glennis A and Ebert, Peter and Audano, Peter A and Loftus, Mark and Porubsky, David and Ebler, Jana and Yilmaz, Feyza and Hallast, Pille and Prodanov, Timofey and Yoo, DongAhn and others},
  title   = {Complex genetic variation in nearly complete human genomes},
  year    = {2025},
  journal = {Nature},
  volume  = {644},
  pages   = {430--441},
  doi     = {10.1038/s41586-025-09140-6},
  pmid    = {40702183}
}
```

---

**[jeong2025segdup]** Jeong H, Dishuck PC, Yoo D, Harvey WT, Munson KM, et al. 2025. Structural polymorphism and diversity of human segmental duplications. *Nature Genetics*. DOI:[10.1038/s41588-024-02051-8](https://doi.org/10.1038/s41588-024-02051-8)

Claim it supports/contradicts:
- **Supports** draft P2 (Methods): "we extracted 18,827 telomere-anchored 500 kb flanks… wfmash all-vs-all at 95% minimum identity" — the Jeong et al. paper characterises 173.2 Mb of duplicated sequence across 170 HPRC-derived haplotype-resolved assemblies (85 samples: 38 African, 47 non-African), demonstrating that the majority of autosomal SDs are now fully resolved by long-read assembly and showing extensive African-enriched haplotype diversity. This directly supports the claim that HPRC v2 assemblies capture population-scale SD structure previously inaccessible to short reads.
- **Supports** report §12 novel-contribution #1 ("Population-scale community structure"): the SD-level population genetics survey confirms that high-identity duplicated sequences show strong population stratification, coherent with the PHR Fst analysis.
- **Supports** draft §Limitations claim that "HPRC v2 achieves near-complete coverage for most chromosome arms": the paper shows that 170 assemblies from 85 samples fully resolve autosomal SDs at 173.2 Mb, providing independent validation of assembly completeness at segmental-duplication loci.

Recommendation: **ADD** to support C4 and C6. Natural citation in the population-genetic and assembly-quality context of topic 12.

Bibtex stanza:
```bibtex
@article{jeong2025segdup,
  author  = {Jeong, Hyeonsoo and Dishuck, Philip C and Yoo, DongAhn and Harvey, William T and Munson, Katherine M and others},
  title   = {Structural polymorphism and diversity of human segmental duplications},
  year    = {2025},
  journal = {Nature Genetics},
  volume  = {57},
  pages   = {390--401},
  doi     = {10.1038/s41588-024-02051-8}
}
```

---

**[gao2023chinesepangenome]** Gao Y, Yang X, Chen H, Tan X, Yang Z, et al. (Chinese Pangenome Consortium). 2023. A pangenome reference of 36 Chinese populations. *Nature*. DOI:[10.1038/s41586-023-06173-7](https://doi.org/10.1038/s41586-023-06173-7)

Claim it supports/contradicts:
- **Supports** draft P1 open-questions paragraph: "Broader population sampling — as envisioned by the Global Pangenome Alliance and related consortia — would allow asking whether specific PHR communities show signatures of local adaptation" (report §12 "Open questions"). The Chinese Pangenome Consortium Phase 1 delivers 116 haplotype-phased assemblies from 58 samples across 36 Chinese minority ethnic groups (average 30.65× HiFi coverage), demonstrating that Asian-ancestry groups underrepresented in HPRC can generate population-specific pangenomes with comparable quality and increased variant detection over GRCh38 and T2T-CHM13.
- **Supports** the framing in report §12 that HPRC v2's limited Oceanian/Middle Eastern representation is a known gap: CPC Phase 1 is the clearest example of a consortium acting to fill one such gap.
- **Supports** hprc_wang2022 claim about the 350-individual HPRC target being a first phase; CPC shows parallel efforts are live.

Recommendation: **ADD** as a representative population-coverage companion to the HPRC lineage, cited in the "Open questions" and "Population representation" subsections of topic 12.

Bibtex stanza:
```bibtex
@article{gao2023chinesepangenome,
  author  = {Gao, Yang and Yang, Xiaofei and Chen, Hao and Tan, Xinjiang and Yang, Zhaoqing and others and {Chinese Pangenome Consortium}},
  title   = {A pangenome reference of 36 {Chinese} populations},
  year    = {2023},
  journal = {Nature},
  volume  = {619},
  pages   = {112--121},
  doi     = {10.1038/s41586-023-06173-7}
}
```

---

**[rausch2025lrpop]** Rausch T, Marschall T, Korbel JO. 2025. The impact of long-read sequencing on human population-scale genomics. *Genome Research*. PMID:40228902 / DOI:[10.1101/gr.280120.124](https://doi.org/10.1101/gr.280120.124)

Claim it supports/contradicts:
- **Supports** report §12 "Open questions" paragraph on scaling long-read to large populations: "challenges remain in scaling long-read technologies to large populations due to cost, computational complexity, and the lack of tools to facilitate the efficient interpretation of SVs in graphs." This is a direct synthesis of the current state, with the review authors (Rausch, Marschall, Korbel) having led both HGSVC2 and EMBL's contributions to HPRC.
- **Supports** the framing in draft P1 that "per-chromosome alignment frames hid trans-chromosomal sequence sharing" — the review characterises reference-anchored SV calling as the dominant paradigm that HPRC challenges.

Recommendation: **ADD** as a background/context citation for the state-of-the-field commentary in topic 12 "Open questions". MEDIUM relevance: useful for the companion paper Methods context but not essential to specific claims.

Bibtex stanza:
```bibtex
@article{rausch2025lrpop,
  author  = {Rausch, Tobias and Marschall, Tobias and Korbel, Jan O},
  title   = {The impact of long-read sequencing on human population-scale genomics},
  year    = {2025},
  journal = {Genome Research},
  volume  = {35},
  pages   = {593--598},
  doi     = {10.1101/gr.280120.124},
  pmid    = {40228902}
}
```

---

**[kulmanov2025jasapage]** Kulmanov M, Ashouri S, Liu Y, Abdelhakim M, Alsolme E, Nagasaki M, et al. 2025. Phased genome assemblies and pangenome graphs of human populations of Japan and Saudi Arabia. *Scientific Data*. PMID:40796583 / DOI:[10.1038/s41597-025-05652-y](https://doi.org/10.1038/s41597-025-05652-y)

Claim it supports/contradicts:
- **Supports** report §12 "Population representation" subsection: "HPRC v2's 232 individuals prioritize AFR, EAS, EUR, SAS and AMR ancestry but have limited representation from Oceania, the Middle East, and isolated indigenous populations." JaSaPaGe constructs a 19-individual (9 Saudi + 10 Japanese) phased pangenome showing comparable or superior variant detection to HPRC for these underrepresented groups, and demonstrates population-specific variants missed by GRCh38 and T2T-CHM13.
- Motivates the framing that population-specific pangenomes (beyond HPRC) are needed to fully characterise subtelomeric diversity in non-European groups.

Recommendation: **ADD** as supporting evidence for the population-representativeness gap in HPRC v2, cited in the "Open questions" subsection. MEDIUM relevance.

Bibtex stanza:
```bibtex
@article{kulmanov2025jasapage,
  author  = {Kulmanov, Maxat and Ashouri, Saeideh and Liu, Yang and Abdelhakim, Marwa and Alsolme, Ebtehal and Nagasaki, Masao and Ohkawa, Yasuyuki and Suzuki, Yutaka and Tawfiq, Rund and Tokunaga, Katsushi and Katayama, Toshiaki and Abedalthagafi, Malak S and Hoehndorf, Robert and Kawai, Yosuke},
  title   = {Phased genome assemblies and pangenome graphs of human populations of {Japan} and {Saudi Arabia}},
  year    = {2025},
  journal = {Scientific Data},
  volume  = {12},
  pages   = {1316},
  doi     = {10.1038/s41597-025-05652-y},
  pmid    = {40796583}
}
```

---

## Section 4: CONTRADICTIONS

No papers in the 2023–2026 search window directly contradict any specific claim made in the Nature draft or in the topic-12 review. Three potential tension points were examined and resolved:

1. **466 vs. 465 haplotypes (C3):** The Logsdon 2025 HGSVC3 paper (Nature 644:430–441) makes no claim about HPRC v2 haplotype counts and cannot contradict C3. The existing review already documents that the 466 vs. 465 discrepancy is a naming convention (CHM13 counted differently), not a data error.

2. **HPRC pangenome coverage completeness:** Jeong et al. 2025 (Nature Genetics) reports that 170 assemblies from 85 samples fully resolve autosomal SDs. This is consistent with (not contradicting) the HPRC v2 near-T2T claim. The acrocentric short arms and sex chromosomes remain partially incomplete even in the Jeong et al. dataset, which also aligns with the draft's Limitations section.

3. **Short-read recombination vs. cross-arm affinity (Lalli 2025 anti-correlation):** No new paper in the searched window re-tests or contradicts the analysis in Lalli 2025 (already cited in the draft). The correlation collapses once the 7 low-callability arms are excluded, as reported in the current manuscript.

---

## Section 5: Search Audit Trail

### Tools Used

| Tool | Description |
|------|-------------|
| `mcp__claude_ai_PubMed__search_articles` | Primary PubMed keyword searches |
| `mcp__claude_ai_PubMed__get_article_metadata` | Article metadata retrieval by PMID |
| `mcp__claude_ai_bioRxiv__search_preprints` | bioRxiv date+category search |
| `mcp__claude_ai_bioRxiv__get_preprint` | Individual preprint DOI lookup |
| OpenAlex batch_lookup (skill) | DOI-based batch metadata + abstract retrieval |
| OpenAlex search_works (skill) | Keyword search with relevance ranking |

### Query Strings and Date Filters

All PubMed queries filtered `2023/01/01:2026/12/31` unless noted.

1. `Human Pangenome Reference Consortium pangenome assembly` → 58 total hits; 30 retrieved
2. `draft human pangenome reference haplotype population diverse` → 3 total hits
3. `Liao Wang Phillippy pangenome reference human genome diversity assembly 2023 2024` → 0 hits (overly specific)
4. `complex genetic variation nearly complete human genomes pangenome 2025` → 64 total hits; 10 retrieved (date 2024–2026)
5. `impact long-read sequencing human population scale genomics structural variation` → 7 hits; all retrieved
6. `haplotype-resolved genome assembly population pangenome structural variation long-read` → 13 hits; all retrieved
7. `Logsdon Eichler complex structural variation complete human genomes haplotype 2025` → 1 hit (= PMID 40702183, confirmed HGSVC3)
8. `Vollger segmental duplication population structural polymorphism pangenome` → 1 hit (= PMID 42000714, 22q11.2, tangential)
9. `phased genome assembly population Japan Saudi Arabia pangenome 2025` → 1 hit (= PMID 40796583, JaSaPaGe)
10. `pangenome graph GWAS genome-wide association disease variant discovery` → 0 hits
11. `structural polymorphism diversity segmental duplications HPRC pangenome 2024` → 0 hits (MeSH mismatch)
12. `HPRC pangenome reference 232 individuals haplotype-resolved diverse 2025 2026` → 0 hits
13. `pangenome reference high-quality haplotype assemblies diverse populations 2025 2026` → 0 hits (HPRC v2 paper not yet indexed)
14. bioRxiv category=genomics, date 2024-01-01 to 2026-05-17 → 50 results, no keyword filter available; scanned for HPRC-adjacent preprints
15. bioRxiv DOI lookup: `10.1101/2025.05.12.653384` (guessed HPRC v2 DOI) → not found
16. OpenAlex batch lookup: 6 target DOIs → confirmed publication metadata for 5 papers
17. OpenAlex keyword search: "Human Pangenome Reference Consortium population haplotype-resolved assembly" 2023–2026 → 20 results, 5 high-citations relevant
18. OpenAlex keyword search: "HPRC pangenome haplotype assembly diverse population 2025" 2024–2026 → 10 results

### Hit Counts After Relevance Filter

| Search | Total hits | After relevance filter | Notes |
|--------|-----------|----------------------|-------|
| PubMed HPRC consortium assembly | 58 | 5 | Kept HPRC-adjacent only |
| PubMed haplotype-resolved LRS | 13 | 4 | Discarded plant/animal genomes |
| PubMed complex genetic variation complete genomes | 64 | 2 | Only human + pangenome-related |
| PubMed long-read population-scale genomics | 7 | 2 | Kept human SV reviews |
| OpenAlex HPRC keyword | 20 | 3 | Discarded non-human, non-SV |
| OpenAlex HPRC 2025 keyword | 10 | 5 | Selected PMID >20 cites |
| **Total unique candidates evaluated** | **~175** | **10** | |
| **Strong recommendations** | | **5** | Sections 3–4 above |

### Dropped Papers (1-line rationale per paper)

| Paper | Reason dropped |
|-------|---------------|
| Ensembl 2025 (NAR 2024) | Database infrastructure update; not relevant to HPRC assembly claims |
| Complete Y chromosome (Nature 2023) | Already tracked in acrocentric/sex-chromosome topics; tangential to topic-12 |
| NextDenovo (Genome Biology 2024) | Non-HPRC assembler; no direct relevance to HPRC v1/v2 lineage |
| Maize T2T pangenome (Nature Genetics 2023) | Plant genome; no relevance |
| Pangenome graphs in biodiversity genomics (NatGen 2025) | Review of non-human pangenomes; tangential |
| Near-complete Middle Eastern genomes (NatGen 2025) | Primarily autozygosity/disease focus; population pangenome coverage overlap with JaSaPaGe already captured |
| vg Giraffe long-read mapper (bioRxiv 2025) | Mapping tool; supports pangenome utility but no new population assembly claims |
| DeepPolisher (Genome Research 2025) | Assembly polishing tool; supports HPRC v2 QV improvement but tool paper, not population pangenome |
| Beyond single references: pangenome in genomic medicine (Front Genet 2025) | Broad review; lower specificity than Rausch et al. 2025 |
| Scalable Nanopore sequencing (NatMeth 2023, PMID 37710018) | Protocol paper for ONT; relevant but superseded by Rausch review for topic-12 purposes |
| 22q11.2 duplication population differences (NatComm 2026, PMID 42000714) | Population-specific SD architecture; valuable but locus-specific, not HPRC-lineage |
| Porubsky 2026 JaSaPaGe alternative (was searching for Vollger) | Returned wrong paper; discarded |
| Recurrent amylase evolution (Science 2024, Nature 2024) | Amylase locus; relevant to concerted evolution topic but not to HPRC population pangenome lineage |
| Kolmogorov Napu ONT (NatMeth 2023) | Protocol paper; no HPRC claims |
| IGLoo immunoglobulin assembly (Cell Rep Methods 2025) | IG-locus-specific tool; not population pangenome |
| Rare disease pangenomics review (EJHG 2026) | Clinical review; background only, lower priority than Rausch 2025 |

### Note on HPRC v2 Main Paper

The `hprc_hprcv2_2025` entry in REFERENCES_v3.bib is currently `@misc` with note "Pre-publication data release, manuscript in submission. Announced May 12, 2025." Despite 13 PubMed searches and 3 OpenAlex queries, the HPRC v2 main manuscript was not found indexed in PubMed or OpenAlex as of this search (May 2026). The paper may be in press, published under a different title, or not yet indexed. **Recommendation for the integrator**: before finalising REFERENCES_v4.bib, check [humanpangenome.org](https://humanpangenome.org) and bioRxiv `10.1101/2025.*` for a deposited preprint; upgrade `@misc` to `@article` with full citation details.
