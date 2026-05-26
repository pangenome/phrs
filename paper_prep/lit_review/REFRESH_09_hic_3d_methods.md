# REFRESH_09: hic_3d_methods — Literature Refresh (2023–2026)

**Generated:** 2026-05-17
**Agent:** lit-refresh-09-hic (agent-63)
**Scope:** New papers 2023–2026 relevant to topic_09_hic_3d_methods

---

## Section 1: Topic Scope

Topic 09 covers the Hi-C, Pore-C, CiFi, Dip-C, and single-cell 3D methods that provide the 3D validation spine of the paper. In the Nature draft it supplies the evidence for claims **C7** (main text paragraphs P5–P6): "Bulk and single-cell Hi-C in GM12878 and sperm, together with mouse meiotic Hi-C peaking at zygotene (Mantel ρ=0.715), tie sequence homology to nuclear-envelope proximity through the meiotic bouquet." It underpins the meiotic bouquet mechanistic explanation for subtelomeric community-based nuclear co-localization. It also provides the methodological substrate for the multi-technology validation table (B/W ratios, Mantel correlations, O/E enrichments) summarized in P5–P6, the multi-mapping artifact control in P6 (the flanking-unique-sequence paradox), the single-cell distance enrichments in Dip-C and sperm (P6), and the mouse zygotene-peak result (P7). The literature sits at the intersection of 3D genome methodology, meiotic chromosome biology (LINC complex / bouquet), and the T2T reference improvement that unlocks the subtelomeric signal.

Draft paragraphs: **P5** (main text L34, 3D community enrichment results), **P6** (main text L36, flanking paradox, single-cell, bouquet mechanism), **P7** (main text L40, mouse meiotic result and zygotene peak).

Report sections: `end-to-end-report/report/05_hic_validation.md`, `06_dipc_validation.md`, `07_integrated.md`, `09_rpe1_self.md`.

---

## Section 2: Existing Citations — Current vs Superseded

All bibkeys from `topic_09_hic_3d_methods.bib` evaluated below. Cross-checked against the 2023–2026 new literature.

| Bibkey | Status | Reason |
|---|---|---|
| `hic3d_liebermanaiden2009` | STILL CURRENT | Foundational Hi-C paper; no supersession possible |
| `hic3d_dixon2012` | STILL CURRENT | TAD discovery paper; still the canonical primary reference |
| `hic3d_imakaev2012` | STILL CURRENT | ICE normalization; still the standard approach used in hicexplorer |
| `hic3d_nagano2013` | STILL CURRENT | First single-cell Hi-C; still the canonical reference for scHi-C rationale |
| `hic3d_rao2014` | STILL CURRENT | GM12878 kilobase Hi-C; still the cell-line baseline |
| `hic3d_stevens2017` | STILL CURRENT | Single-cell 3D structures; still relevant for interpreting Dip-C findings |
| `Tan2018` | STILL CURRENT | Dip-C; primary dataset used in the draft |
| `hic3d_alavattam2019` | STILL CURRENT | Mouse meiotic Hi-C with attenuated compartments; not superseded — but **new papers in 2023 now extend the stage-resolution and female sex coverage** (see Section 3) |
| `Patel2019` | STILL CURRENT | Mouse meiotic Hi-C loop sizes; not superseded; new papers confirm and extend findings |
| `hic3d_deshpande2022` | STILL CURRENT | Pore-C; primary methodology paper; still current |
| `Zuo2021` | STILL CURRENT | Stage-resolved mouse meiotic Hi-C; primary dataset for the zygotene peak result |
| `Xu2025` | STILL CURRENT | Sperm single-cell Hi-C; primary dataset for the sperm 3D result |
| `hic3d_cifi2025` | STILL CURRENT | CiFi; primary methodology paper |
| `hic3d_nurk2022` | STILL CURRENT | T2T-CHM13; methodological prerequisite |
| `hic3d_wolff2018` | STILL CURRENT | HiCExplorer; still the pipeline used |
| `hic3d_harper2004` | STILL CURRENT | Bouquet review; still the canonical description of bouquet formation across eukaryotes |
| `hic3d_koszul2008` | STILL CURRENT | Meiotic chromosome movements (yeast); still foundational |
| `hic3d_scnanoHiC2023` | STILL CURRENT | scNanoHi-C method paper; used in Methods |
| `hic3d_vara2019` | STILL CURRENT | Vara et al. 2019 spermatogenesis 3D structure; still relevant for background — **note: the 2025 RAD21L paper (same lab) extends this work** |
| `hic3d_scnanoHiC2_2025` | STILL CURRENT | scNanoHi-C2 applied to embryonic germ cells; used in Methods |

**Summary:** No existing citation is superseded. Three are extended by new 2023+ work (Alavattam 2019, Patel 2019, Vara 2019), and those new papers are added as companions, not replacements.

---

## Section 3: NEW Papers to Add (STRONG Relevance, 2023+)

### [1] STRONG relevance

**[hic3d_marinGual2025]** Marín-Gual L et al. (Vara C, Ruiz-Herrera A), 2025. Meiotic cohesin RAD21L shapes 3D genome structure and transcription in the male germline. *Science Advances* 11(40):eadv2283. PMID:41032613 / DOI:[10.1126/sciadv.adv2283](https://doi.org/10.1126/sciadv.adv2283)

**Claim it supports/contradicts:**
- Report §07 (meiotic bouquet section): "this alignment is not merely a manifestation of the transient bouquet conformation… it depends on the force-transmitting LINC complex (which modulates alignment range but does not alter loop sizes)."
- Draft P6 (testable prediction, implicit): "The highest-priority unresolved question is whether the bouquet mechanism specifically drives the subtelomeric community-based proximity signal."

**Why STRONG:** This paper performs a cohesin loss-of-function (RAD21L knockout) in mouse spermatocytes and directly measures the 3D genome. The key finding is that RAD21L deletion **disrupts bouquet formation** and **increases** telomeric inter-chromosomal interactions in primary spermatocytes. This is a direct mechanistic test of the claim that bouquet organization mediates inter-chromosomal telomeric contacts: when bouquet-associated cohesin is removed, inter-heterologous-chromosome telomeric contacts go up, not down — consistent with loss of directed proximity (community-structured co-localization) and gain of non-specific telomeric aggregation. This is the closest published analog to the SUN1-W151R LINC mutant experiment that the draft names as a testable prediction. The paper is from the same group as hic3d_vara2019.

**Recommendation:** ADD. Cite in draft P6 alongside the meiotic bouquet mechanistic sentence ("A mechanism is available. Meiotic prophase I telomeres are tethered to the nuclear envelope…"). Also cite in the testable-predictions section of report §07 as evidence that cohesin/bouquet perturbation changes inter-chromosomal telomeric contact patterns.

```bibtex
@article{hic3d_marinGual2025,
  author    = {Mar{\'i}n-Gual, Laia and Vara, Covadonga and Sainz-Urruela, Raquel and
               Cuartero, Yasmina and {\'A}lvarez-Gonz{\'a}lez, Luc{\'i}a and
               Felipe-Medina, Natalia and Garcia, Francisca and Llano, Elena and
               Marti-Renom, Marc A. and Pend{\'a}s, Alberto M. and
               Ruiz-Herrera, Aurora},
  title     = {Meiotic cohesin {RAD21L} shapes {3D} genome structure and
               transcription in the male germline},
  year      = {2025},
  journal   = {Science Advances},
  volume    = {11},
  number    = {40},
  pages     = {eadv2283},
  doi       = {10.1126/sciadv.adv2283},
  pmid      = {41032613},
}
```

---

### [2] STRONG relevance

**[hic3d_he2023]** He J, Yan A, Chen B, Huang J, Kee K, 2023. 3D genome remodeling and homologous pairing during meiotic prophase of mouse oogenesis and spermatogenesis. *Developmental Cell* 58(24):3009–3027.e6. PMID:37963468 / DOI:[10.1016/j.devcel.2023.10.009](https://doi.org/10.1016/j.devcel.2023.10.009)

**Claim it supports/contradicts:**
- Report §07 meiotic bouquet section: "Average meiotic chromatin loop sizes increase from ~500 kb at leptotene to ~700 kb at zygotene…" (Zuo 2021; Patel 2019 reported slightly larger values).
- Draft P7 (mouse meiotic result): the zygotene-peak result rests on the Zuo 2021 stage-sorted mouse Hi-C. He et al. provide an independent stage-resolved Hi-C (both sexes), with the key finding that homolog pairing changes strategy at the **leptotene-to-zygotene transition** — from LINE-enriched B-compartment contacts to SINE-enriched A-compartment contacts — which places the critical timing of subtelomeric contact reorganization at exactly the stage where the Mantel ρ peaks.
- Specifically: compartments and TADs "gradually disappeared and slowly recovered in both sexes," confirming the Alavattam 2019 observation and extending it to the female germline.

**Why STRONG:** An independent second stage-resolved meiotic Hi-C dataset in mouse (He et al. 2023, Tsinghua) that confirms the stage-specific chromatin reorganization backbone. This strengthens the claim that the Zuo 2021 zygotene-peak Mantel result is interpretable against a well-established compartment/TAD reorganization context. The pairing-strategy switch at leptotene–zygotene (from B- to A-compartment contacts) is directly relevant: PHRs may reside in or near B compartments at leptotene and shift to A contacts at zygotene, consistent with the bouquet clustering driving subtelomeric proximity.

**Recommendation:** ADD. Cite in report §08 (mouse section) and/or the meiotic bouquet section of §07 as a second independent confirmation of the leptotene-to-zygotene chromatin reorganization relevant to the zygotene Hi-C peak.

```bibtex
@article{hic3d_he2023,
  author    = {He, Jing and Yan, An and Chen, Bo and Huang, Jiahui and Kee, Kehkooi},
  title     = {{3D} genome remodeling and homologous pairing during meiotic prophase
               of mouse oogenesis and spermatogenesis},
  year      = {2023},
  journal   = {Developmental Cell},
  volume    = {58},
  number    = {24},
  pages     = {3009--3027.e6},
  doi       = {10.1016/j.devcel.2023.10.009},
  pmid      = {37963468},
}
```

---

### [3] STRONG relevance

**[hic3d_cheng2024]** Cheng G, Pratto F, Brick K, Li X, Alleva B, Huang M, Lam G, Camerini-Otero RD, 2024. High resolution maps of chromatin reorganization through mouse meiosis reveal novel features of the 3D meiotic structure. *bioRxiv* (NIDDK/NIH preprint). PMID:38903112 / DOI:[10.1101/2024.03.25.586627](https://doi.org/10.1101/2024.03.25.586627)

**Claim it supports/contradicts:**
- Report §07 meiotic bouquet section: "Average meiotic loop sizes increase… Most PHR sequences (105 kb median) would reside at the base of a single meiotic loop at leptotene."
- Draft P6 mechanistic sentence (the meiotic loop–PHR size argument).

**Why STRONG:** This NIH preprint provides the **highest temporal and spatial resolution Hi-C to date for mouse spermatogenesis**, using Micro-C–level resolution in stage-pure sorted cell populations from spermatogonia through all of meiotic prophase I. Key finding directly relevant to the draft's mechanistic argument: **CTCF sites anchor meiotic loops even when TADs are abolished**, and CTCF localizes to the meiotic axis (confirmed by ChIP), placing CTCF-anchored loop bases at the axis — and since the axis is at the nuclear envelope attachment points — at the telomere. This directly supports the claim that PHRs (105 kb median) fit at the base of a meiotic loop. Additionally, regulatory element contacts are preserved during meiosis despite TAD dissolution, confirming that the loss of bulk TADs is not a loss of all organized chromatin structure. The preprint is from the Camerini-Otero group at NIDDK, a major meiosis lab; its publication should be tracked.

**Recommendation:** ADD (preprint — flag for peer-reviewed version check before final submission). Cite in the meiotic bouquet section of §07 and in draft P6 alongside Patel 2019 as further mechanistic support for the loop-size argument. Note preprint status.

```bibtex
@article{hic3d_cheng2024,
  author    = {Cheng, Gang and Pratto, Florencia and Brick, Kevin and Li, Xin and
               Alleva, Benjamin and Huang, Mini and Lam, Gabriel and
               Camerini-Otero, R. Daniel},
  title     = {High resolution maps of chromatin reorganization through mouse meiosis
               reveal novel features of the {3D} meiotic structure},
  year      = {2024},
  journal   = {bioRxiv},
  doi       = {10.1101/2024.03.25.586627},
  pmid      = {38903112},
  note      = {Preprint; NIDDK/NIH},
}
```

---

### [4] MEDIUM relevance (include as background)

**[hic3d_kitamura2025]** Kitamura Y, Namekawa SH, 2025. The 3D genome during germline development and meiosis. *Trends in Genetics* 42(3):255–267. PMID:41407613 / DOI:[10.1016/j.tig.2025.11.004](https://doi.org/10.1016/j.tig.2025.11.004)

**Claim it supports/contradicts:**
- General background for topic 09 — synthesizes all recent progress in 3D genome organization in germlines, including meiotic chromosome dynamics.

**Why MEDIUM:** This is a current-state-of-the-field review (Namekawa lab, major meiosis group at UC Davis). It will be a standard reference for reviewers asking about the 3D germline landscape. It covers scNanoHi-C2, the new meiotic Hi-C datasets, and explicitly discusses "germline-specific 3D genome features" and "genome evolution" — directly aligning with the claims made in topic 09.

**Recommendation:** ADD as background citation. Cite in the era overview section of topic 09 and/or in report §07 as a synthesizing review.

```bibtex
@article{hic3d_kitamura2025,
  author    = {Kitamura, Yuka and Namekawa, Satoshi H.},
  title     = {The {3D} genome during germline development and meiosis},
  year      = {2025},
  journal   = {Trends in Genetics},
  volume    = {42},
  number    = {3},
  pages     = {255--267},
  doi       = {10.1016/j.tig.2025.11.004},
  pmid      = {41407613},
}
```

---

### [5] MEDIUM relevance (include as background)

**[hic3d_yin2024]** Yin L, Jiang N, Li T, Zhang Y, Yuan S, 2024. Telomeric function and regulation during male meiosis in mice and humans. *Andrology* 13(5):1170–1180. PMID:38511802 / DOI:[10.1111/andr.13631](https://doi.org/10.1111/andr.13631)

**Claim it supports/contradicts:**
- Draft P6 bouquet mechanism sentence: "Meiotic prophase I telomeres are tethered to the nuclear envelope by the MAJIN, TERB2 and TERB1 complex."
- Cites the LINC complex (SUN-KASH), SPDYA-CDK2, and TTM trimer (TERB1-TERB2-MAJIN) as critical regulators.

**Why MEDIUM:** Updated review covering the current state of LINC/bouquet/TERB complex biology in mice and humans, including clinical implications (azoospermia mutations). Relevant because it covers the human-specific aspects of bouquet regulation that topic 09 invokes. The draft's mechanistic paragraph (draft P6) lists MAJIN-SUN1-KASH5 but does not cite a recent human-focused review; this paper fills that gap.

**Recommendation:** ADD as background. Cite alongside `hic3d_harper2004` in the bouquet mechanism section.

```bibtex
@article{hic3d_yin2024,
  author    = {Yin, Lisha and Jiang, Nan and Li, Tao and Zhang, Youzhi and Yuan, Shuiqiao},
  title     = {Telomeric function and regulation during male meiosis in mice and humans},
  year      = {2024},
  journal   = {Andrology},
  volume    = {13},
  number    = {5},
  pages     = {1170--1180},
  doi       = {10.1111/andr.13631},
  pmid      = {38511802},
}
```

---

### [6] MEDIUM relevance (include as background)

**[hic3d_thadani2026]** Thadani R, Johnson N, Cooper JP, 2026. Chromosome Ends in Motion: Telomeres as Hazards and Hubs in Meiosis. *Cold Spring Harbor Perspectives in Biology* 18(2). PMID:41419316 / DOI:[10.1101/cshperspect.a041705](https://doi.org/10.1101/cshperspect.a041705)

**Claim it supports/contradicts:**
- Draft P6 bouquet mechanism: "the telomere bouquet orchestrates movements of meiotic chromosomes that facilitate pairing and recombination between homologous chromosomes."
- This review covers both canonical bouquet function and newly discovered meiotic telomere functions (hazards — non-homologous recombination risk) — directly relevant to the paper's thesis that telomere clustering creates risk of ectopic recombination.

**Why MEDIUM:** The Cooper lab (CU Anschutz) specifically covers the dual role of telomeres in meiosis as both "hubs" (bouquet function) and "hazards" (potential for non-homologous recombination). This is conceptually aligned with the draft's core argument: the same clustering that facilitates homologous pairing also pre-positions community arms for ectopic exchange. Using this framing would strengthen the mechanistic justification in draft P6.

**Recommendation:** ADD as background. Cite in the meiotic bouquet section of topic 09 alongside `hic3d_harper2004`.

```bibtex
@article{hic3d_thadani2026,
  author    = {Thadani, Rahul and Johnson, Noah and Cooper, Julia Promisel},
  title     = {Chromosome Ends in Motion: Telomeres as Hazards and Hubs in Meiosis},
  year      = {2026},
  journal   = {Cold Spring Harbor Perspectives in Biology},
  volume    = {18},
  number    = {2},
  doi       = {10.1101/cshperspect.a041705},
  pmid      = {41419316},
}
```

---

### [7] MEDIUM relevance (include as background)

**[hic3d_xie2025]** Xie W, Gowder M, Bazzano D, DeSantis M, Hammoud SS, 2025. Rewiring for movements in meiotic prophase: regulators, roles, and evolutionary pathways. *Current Opinion in Genetics & Development* 93:102366. PMID:40484002 / DOI:[10.1016/j.gde.2025.102366](https://doi.org/10.1016/j.gde.2025.102366)

**Claim it supports/contradicts:**
- Draft P6 (LINC complex): "Meiotic prophase I telomeres are tethered to the nuclear envelope by the MAJIN, TERB2 and TERB1 complex, which drives the bouquet stage of zygotene chromosome organisation."
- Also relevant to the "open questions" section of topic 09 on LINC complex mutant Hi-C.

**Why MEDIUM:** The Hammoud lab (Michigan) review covers meiotic prophase movements, RPMs, LINC complex hierarchy, and evolutionary variation in bouquet mechanisms. Relevant to the draft's mechanistic explanation and the testable predictions regarding LINC complex perturbation.

**Recommendation:** ADD as background. Cite in the LINC complex / bouquet mechanism text.

```bibtex
@article{hic3d_xie2025,
  author    = {Xie, Wenxin and Gowder, Manjunath and Bazzano, Dominic and
               DeSantis, Morgan and Hammoud, Saher S.},
  title     = {Rewiring for movements in meiotic prophase: regulators, roles,
               and evolutionary pathways},
  year      = {2025},
  journal   = {Current Opinion in Genetics \& Development},
  volume    = {93},
  pages     = {102366},
  doi       = {10.1016/j.gde.2025.102366},
  pmid      = {40484002},
}
```

---

## Section 4: Contradictions

**No papers found in the searched window (2023–2026) that directly contradict the draft's claims in topic 09.**

Specific checks performed:
- The claim that meiotic bouquet clustering positions subtelomeric regions for ectopic exchange: no paper contradicts this. The Marín-Gual 2025 RAD21L paper is consistent with bouquet-mediated co-organization (increased non-specific telomeric contacts when bouquet is disrupted by cohesin loss).
- The claim that TADs are lost in meiosis (Alavattam 2019, Zuo 2021): confirmed independently by He et al. 2023 and Cheng et al. 2024.
- The claim that CTCF/cohesin still anchors meiotic loops even without TADs: Cheng et al. 2024 confirms this (CTCF at loop bases on axes), does not contradict it.
- The claim that PHR loop-size argument (105 kb median fits within one meiotic loop): Cheng et al. 2024 Micro-C data confirms loop-base CTCF anchoring is present through meiotic prophase, and loop sizes reported are still consistent with the prior Zuo/Patel estimates (no contradictions on loop size).
- The Lalli et al. cM/Mb anti-correlation (mentioned in draft P9 as collapsing once low-callability arms removed): no new paper from 2023–2026 revives or re-establishes that anti-correlation.

---

## Section 5: Search Audit Trail

### Tools used
- **mcp__claude_ai_PubMed__search_articles** (PubMed MCP) — primary search tool with PMID and DOI verification
- **mcp__claude_ai_PubMed__get_article_metadata** — for full abstract and DOI verification of all candidates
- **openalex-database skill** — cross-validation queries for cited-by counts and additional discovery
- **mcp__claude_ai_bioRxiv__search_preprints** — searched by category (genomics, 2023-01-01 to 2026-05-01); returned general genomics preprints not keyword-searchable; supplementary only

### Date range filter
All PubMed searches: `date_from: 2023`, `date_to: 2026`.

### Query strings and hit counts

| Query | Hits (raw) | Relevant after filter |
|---|---|---|
| "Hi-C 3D genome chromatin conformation subtelomere OR telomere meiotic bouquet" (2023–2026) | 13 | 5 |
| "single cell Hi-C scHi-C sperm OR meiosis chromosome conformation 3D" (2023–2026) | 6 | 3 |
| "meiotic telomere bouquet LINC complex SUN1 KASH chromosome movement meiosis" (2023–2026) | 1 | 1 |
| "Hi-C 3D genome subtelomere territory" (2023–2026) | 0 | 0 |
| "RAD21L cohesin meiosis 3D genome telomere spermatocyte chromatin" (2023–2026) | 1 | 1 |
| "Pore-C CiFi long-read chromatin conformation repetitive 3D" (2023–2026) | 0 | 0 |
| "scNanoHi-C nanopore single cell chromatin conformation long read" (2023–2026) | 0 | 0 |
| "T2T telomere-to-telomere Hi-C chromatin 3D subtelomeric" (2023–2026) | 0 | 0 |
| "meiotic Hi-C spermatocyte prophase zygotene chromosome organization" (2023–2026) | 0 | 0 |
| "human meiosis Hi-C chromosome conformation prophase spermatocyte" (2023–2026) | 0 | 0 |
| "HiCExplorer pipeline normalization 2023 2024" (2023–2026) | 0 | 0 |
| "meiotic chromosome movement rapid prophase I SUN1 MAJIN telomere nuclear envelope" (2023–2026) | 0 | 0 |
| OpenAlex keyword searches (4 queries) | ~20 | 0 new (confirmed overlapping with PubMed set) |

**Total hits before relevance filter:** ~41 non-duplicate records.
**Total hits after relevance filter:** 11 records reviewed in depth; 7 included (3 STRONG, 4 MEDIUM); 4 dropped.

### Dropped hits (1-line reason each)

- **PMID 39833200** (Jiménez-Martín et al. 2025, Nat Comms) — Fission yeast (*S. pombe*) centromere/bouquet; too distant from mammalian/human claims in the draft.
- **PMID 36913831** (Mytlis et al. 2023, Curr Opin Cell Biol) — Bouquet MTOC review (zebrafish/mouse); focuses on cytoskeletal zygotene cilium anchoring, not Hi-C or sequence-similarity; tangential.
- **PMID 37181752** (Valero-Regalón et al. 2023, Front Cell Dev Biol) — Marsupial meiosis DSB patterns; too phylogenetically distant from human/mouse focus of the draft.
- **PMID 38520405** (Puerto et al. 2024, Nucleic Acids Res) — Drosophila somatic chromosome pairing (condensin II / Z4); different organism, somatic context, no relevance to meiotic bouquet or Hi-C methods used in draft.
- **PMID 38478471** (You et al. 2024, Plant Physiol) — Rice (*Oryza sativa*) chromosome end pairing in meiosis; plant model, not relevant to human/mouse claims.
- **PMID 39283979** (Zhou et al. 2024, Plant J) — Rice kinesin PSS1, homologous pairing; plant model, too distant.
- **PMID 39013853** (Cromer et al. 2024, Nat Comms) — Arabidopsis meiotic RPMs; plant model.
- **PMID 38809870** (Kameyama et al. 2024, Zoolog Sci) — Medaka TERB1 mutant; fish model, no Hi-C data.
- **PMID 36173570** (Sakashita et al. 2023, Methods Mol Biol) — Bioinformatics pipeline chapter for super-enhancers/3D chromatin in spermatogenesis; methods review, not original new findings.
- **PMID 40715639** (Cai et al. 2025, Nat Plants) — Arabidopsis LINC (SINE3/PSS1); plant model.

### Key gap identified
Human meiotic Hi-C remains absent from the 2023–2026 literature at the stage-resolution and cell-type purity needed to directly test the draft's cross-species inference from mouse zygotene to human prophase-I. The Vara lab has ongoing work in human spermatogenesis (hic3d_vara2019), but no published 2023–2026 human stage-resolved meiotic Hi-C was found in the search window. This gap is stated in the existing review and remains open.
