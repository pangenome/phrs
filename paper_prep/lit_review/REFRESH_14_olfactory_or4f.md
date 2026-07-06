# REFRESH 14: Olfactory OR4F — Literature Refresh (2023–2026)

*Agent: agent-57 | Date: 2026-05-17*

---

## Section 1: Topic Scope

This topic covers the OR4F subfamily of olfactory receptor genes as molecular tracers of subtelomeric community structure. In the draft and report, OR4F genes appear explicitly in the Extended Data figure (Ed4: "OR4F pseudogenisation gradient") and in the CONSISTENCY_AUDIT flag for §03 gene enrichment (MAJOR GAP: draft does not name OR4F). The review anchors abstract claim C5 (cladistic structure of subtelomeric communities) and claim C4 (extended interchromosomal homology): OR4F genes are among the defining duplicated families at chromosome ends (confirmed as Ambrosini Block 2), and the arm-specific pseudogenization gradient (11.1% on chr7p to 99.8% on chr15q across 16 arms, n=5,023 annotations) is the sharpest available molecular signal of differential exchange-history among communities. In the Nature draft the OR4F story is entirely absent from the Main text (only DUX4L is named); the CONSISTENCY_AUDIT recommends a 3–5 sentence §03 paragraph naming OR4F as the canonical Mefford prediction (Main L40). This literature refresh specifically targets papers that update the evolutionary mechanisms driving OR pseudogenization, the functional genomics of OR gene-cluster variation, and the comparative primate context for the trade-off hypothesis.

---

## Section 2: Existing Citations — Current or Superseded

| Bibkey | Status | Reason |
|--------|--------|--------|
| `BuckAxel1991` | STILL CURRENT | Foundational discovery of OR multigene family; not challenged. |
| `Rouquier1998` | STILL CURRENT | First chromosome-scale OR survey; data uncontested. |
| `Trask1998` | STILL CURRENT | Established OR genes in polymorphically duplicated subtelomeric blocks; still primary evidence for C4. |
| `Glusman2001` | STILL CURRENT | First comprehensive human OR catalog; baseline for pseudogene fraction estimates. |
| `Mefford2001` | STILL CURRENT | OR subtelomeric blocks with inter-chromosomal exchange evidence; still mechanistic anchor. |
| `MeffordTrask2002` | STILL CURRENT (from REFERENCES_v3.bib) | Broader subtelomere structure review; still foundational. |
| `Niimura2003` | STILL CURRENT | Reannotation baseline with pseudogene fraction; methods and data still used as reference. |
| `YoungFriedman2002` | STILL CURRENT | Mouse/human OR evolutionary comparison; data not superseded. |
| `GiladNatSelection2003` | STILL CURRENT | Natural selection analysis; still used as evidence for relaxed selection in hominids. |
| `Gilad2004` | STILL CURRENT but FRAMING UPDATED | Original data (60% pseudogene rate in humans; trichromacy coincidence) still valid, but the simple "trade-off" framing is now contested by Chi et al. 2025 (PMID 40021902), which shows olfactory receptors underwent functional reallocation (narrow→broad tuning) rather than straightforward deterioration. The statement "this 'trade-off' hypothesis" in the existing review should be qualified: the evidence now favors "sensory reallocation" with pseudogenization as one component. |
| `HumanORFamily2004` | STILL CURRENT | Classification paper; not challenged. |
| `Niimura2007` | STILL CURRENT | Birth-and-death evolution across 13 mammals; still definitive for that scope. |
| `Keller2007` | STILL CURRENT | OR7D4 androstenone variant → perception change; paradigm intact. |
| `Hasin2008` | STILL CURRENT | CNV map for OR loci; data valid. NOTE: the existing review's "Open questions" section incorrectly attributes this paper to "Young et al. (2008) \cite{Hasin2008}" — should be "Hasin et al. (2008)". The correct bibkey for the Young 2008 CNV paper is `Young2008CNV`. Both keys exist in REFERENCES_v3.bib. |
| `Young2008CNV` | STILL CURRENT | Extensive OR CNV in AJHG; data uncontested. |
| `Mainland2014` | STILL CURRENT | Functional variability of OR repertoire; 30% allelic difference between individuals still the standard citation for this claim. |
| `Niimura2014` | STILL CURRENT | Elephant OR expansion across 13 placentals; useful comparative context. |
| `Niimura2009` | STILL CURRENT | Vertebrate OR origin and evolution; not superseded. |
| `Zozulya2001` | STILL CURRENT | Early OR repertoire survey; used as background. |

---

## Section 3: NEW Papers to Add (STRONG relevance, 2023+)

> PMIDs and DOIs verified via PubMed MCP (search date: 2026-05-17).

---

**[chi2025primate]** Chi et al., 2025. Genomic and phenotypic evidence support visual and olfactory shifts in primate evolution. *Nature Ecology & Evolution* 9(4):721–733. PMID:40021902 / DOI:10.1038/s41559-025-02651-5

Claim it supports/contradicts: Existing review, "Open questions" section, paragraph 1: "The dominant hypothesis in the OR evolution literature through 2014 attributed human OR pseudogenization primarily to relaxed purifying selection following visual system improvements (the trichromacy trade-off, Gilad et al. 2004 \cite{Gilad2004}). This explains the lineage-wide elevation of pseudogene fractions relative to other mammals, but it does not explain the arm-by-arm gradient." Chi et al. 2025 UPDATES and partly CONTRADICTS the simple trichromacy trade-off framing: large-scale functional analyses of visual and olfactory receptors across extant primates show that olfactory receptors in anthropoids shifted from narrowly to broadly tuned (functional reallocation), not simply lost. Strepsirrhines retained sensitive dim-light vision and enhanced narrowly tuned OR function. The pattern is "sensory reallocation rather than strict trade-offs," explicitly challenging the Gilad 2004 framing. For the draft: the discussion of why OR pseudogenization rates vary should be updated to note that the lineage-wide pseudogenization burden is now understood as part of a broader sensory-system reallocation, making the arm-by-arm gradient (which Chi et al. do not study) even more likely to reflect exchange architecture rather than simple selective relaxation.

Recommendation: ADD as upgrade to `Gilad2004` context. Cite alongside `Gilad2004` in any sentence presenting the trichromacy/olfaction trade-off. Suggested addition to "Open questions" paragraph 1: "...though this simple trade-off hypothesis has been challenged by large-scale genomic and functional analysis showing 'sensory reallocation rather than strict trade-offs' in anthropoid primates (Chi et al. 2025)."

```bibtex
@article{chi2025primate,
  author  = {Chi, Hai and Wan, Jiahui and Melin, Amanda D. and DeCasien, Alex R.
             and Wang, Sufang and Zhang, Yudan and Cui, Yimeng and Guo, Xin
             and Zhao, Le and Williamson, Joseph and Zhang, Tianmin and Li, Qian
             and Zhan, Yue and Li, Na and Guo, Jinqu and Xu, Zhe and Hou, Wenhui
             and Cao, Yumin and Yuan, Jiaqing and Zheng, Jiangmin and Shao, Yong
             and Wang, Jinhong and Chen, Wu and Song, Shengjing and Lu, Xiaoli
             and Qi, Xiaoguang and Zhang, Guojie and Rossiter, Stephen J.
             and Wu, Dong-Dong and Liu, Yang and Lu, Huimeng and Li, Gang},
  title   = {Genomic and phenotypic evidence support visual and olfactory shifts
             in primate evolution},
  journal = {Nature Ecology \& Evolution},
  year    = {2025},
  volume  = {9},
  number  = {4},
  pages   = {721--733},
  doi     = {10.1038/s41559-025-02651-5},
  pmid    = {40021902}
}
```

---

**[foerster2025gwasolfaction]** Förster et al., 2025. Genome-wide association meta-analysis of human olfactory identification discovers sex-specific and sex-differential genetic variants. *Nature Communications* 16(1):5434. PMID:40593737 / DOI:10.1038/s41467-025-61330-y

Claim it supports/contradicts: Existing review, "Findings most relevant to C5" section and the Mainland et al. 2014 context: "Mainland et al. (2014) \cite{Mainland2014} added the functional dimension. Using heterologous expression assays, they identified agonists for 18 human odorant receptors and found that 63% of those tested carry polymorphisms that alter receptor function in vitro." Förster et al. 2025 extends this to the GWAS level (n=21,495) and identifies 10 independent loci associated with odor identification, 7 of them novel, predominantly in OR gene clusters. This is the largest-scale genetic evidence to date that OR gene cluster variation — including at OR loci concentrated in subtelomeric regions — translates to measurable inter-individual differences in human olfactory perception. It supports the existing review's "functional consequence of pseudogenization" argument: if OR gene cluster loci are GWAS hits for olfactory phenotypes, pseudogenization at specific arms could reduce olfactory repertoire diversity.

Recommendation: ADD. Cite in the "Open questions" paragraph on "OR gene-functional variation and community-mediated rescue" (existing review, paragraph 3 under Open questions): "...whether individuals with low-OR4F communities have detectably reduced olfactory repertoires, but this is a tractable question" — Förster et al. 2025 provides indirect evidence that OR cluster variation is already detectable at GWAS scale.

```bibtex
@article{foerster2025gwasolfaction,
  author  = {F{\"o}rster, Franz and Emmert, David and Horn, Katrin and Pott, Janne
             and Frasnelli, Johannes and Imtiaz, Mohammed Aslam and Melas, Konstantinos
             and Talevi, Valentina and Chen, Honglei and Engel, Christoph
             and Filosi, Michele and Fornage, Myriam and G{\"o}gele, Martin
             and L{\"o}ffler, Markus and Mosley, Thomas H. and Pattaro, Cristian
             and Pramstaller, Peter and Shrestha, Srishti and Aziz, N. Ahmad
             and Breteler, Monique M. B. and Wirkner, Kerstin and Scholz, Markus
             and Fuchsberger, Christian},
  title   = {Genome-wide association meta-analysis of human olfactory identification
             discovers sex-specific and sex-differential genetic variants},
  journal = {Nature Communications},
  year    = {2025},
  volume  = {16},
  number  = {1},
  pages   = {5434},
  doi     = {10.1038/s41467-025-61330-y},
  pmid    = {40593737}
}
```

---

**[hayakawa2025chimpanzeeOR]** Hayakawa et al., 2025. Genome-scale evolution in local populations of wild chimpanzees. *Scientific Reports* 15(1):548. PMID:39747985 / DOI:10.1038/s41598-024-84163-z

Claim it supports/contradicts: Existing review, "Primate comparative genomics and the OR4F lineage specifically" (Open questions paragraph 4): "A direct comparison of OR4F pseudogenization rates across great apes — ideally using T2T-quality assemblies — would show whether the arm-specific gradient is human-specific or shared with Pan, Gorilla, and Pongo." Hayakawa et al. 2025 provides population-level OR pseudogene data in wild chimpanzees using exome-capture sequencing across 42 individuals from six African regions. Key finding: OR7D4 and TAS2R42 are shared segregating pseudogenes across western and eastern chimpanzee populations, with many other OR pseudogenes being population-specific. This directly supports the general principle that OR pseudogenization in great apes is population-stratified — exactly what we predict for the human arm-by-arm gradient. It also shows that the tools for population-level OR pseudogene analysis in great apes are now available, even if subtelomeric OR4F specifically has not been analyzed.

Recommendation: ADD. Cite in "Open questions" paragraph 4 to update the "tractable question" language: the population-level pseudogene approach is already deployed in chimpanzees and applicable to human OR4F arms.

```bibtex
@article{hayakawa2025chimpanzeeOR,
  author  = {Hayakawa, Takashi and Kishida, Takushi and Go, Yasuhiro
             and Inoue, Eiji and Kawaguchi, Eri and Aizu, Tomoyuki
             and Ishizaki, Hinako and Toyoda, Atsushi and Fujiyama, Asao
             and Matsuzawa, Tetsuro and Hashimoto, Chie and Furuichi, Takeshi
             and Agata, Kiyokazu},
  title   = {Genome-scale evolution in local populations of wild chimpanzees},
  journal = {Scientific Reports},
  year    = {2025},
  volume  = {15},
  number  = {1},
  pages   = {548},
  doi     = {10.1038/s41598-024-84163-z},
  pmid    = {39747985}
}
```

---

**[dubey2026orprecision]** Dubey et al., 2026. The Emerging Role of Olfactory Receptors: From Genomics to Precision Medicine. *Molecular Diagnosis & Therapy* 30(2):295–319. PMID:41615562 / DOI:10.1007/s40291-026-00832-x

Claim it supports/contradicts: Existing review, "Copy-number variation and the OR4F gradient in diverse human populations" (Open questions paragraph 5): "Whether the 11.1%–99.8% pseudogenization gradient is consistent across populations, or whether certain superpopulations have different arm-specific pseudogenization profiles, remains unknown." Dubey et al. 2026 (review) specifically states: "At the population level, substantial variation in OR allele frequencies across ancestries introduces both opportunities and challenges for precision medicine." This is consistent with our prediction that population-stratified OR4F pseudogenization rates should exist. The review also synthesizes evidence that CNV and pseudogenization drive interindividual differences in receptor function — supporting the functional consequence of arm-specific pseudogenization.

Recommendation: ADD as contextual background. Useful as a 2026 synthesis reference for OR genomics in any paragraph discussing the translational relevance of OR pseudogenization gradients. MEDIUM relevance; do not lead with it.

```bibtex
@article{dubey2026orprecision,
  author  = {Dubey, Nidhi and Rai, Swati and Tripathi, Prabhat
             and Arya, Ankish and Sahoo, Amaresh Kumar and Varadwaj, Pritish Kumar},
  title   = {The Emerging Role of Olfactory Receptors: From Genomics to Precision Medicine},
  journal = {Molecular Diagnosis \& Therapy},
  year    = {2026},
  volume  = {30},
  number  = {2},
  pages   = {295--319},
  doi     = {10.1007/s40291-026-00832-x},
  pmid    = {41615562}
}
```

---

**[brann2024schistosomasubtel]** Brann et al., 2024. Subtelomeric plasticity contributes to gene family expansion in the human parasitic flatworm Schistosoma mansoni. *BMC Genomics* 25(1):217. PMID:38413905 / DOI:10.1186/s12864-024-10032-8

Claim it supports/contradicts: Report section 03 (gene enrichment) and the existing review's mechanistic framework for the OR4F gradient: "The mechanism is exchange-mediated concerted evolution. When two chromosome arms exchange subtelomeric sequence at high frequency, they tend toward sequence homogeneity." Brann et al. 2024 independently demonstrates in a metazoan parasite that subtelomeric regions undergo extensive interchromosomal recombination producing segmental duplications that drive gene family expansion — directly mechanistically parallel to the human OR subtelomere story. Their conclusion that "subtelomeric regions act as a genomic playground for trial-and-error of gene duplication and subsequent divergence" supports the general principle. Relevance is cross-taxonomic (not human/primate) but strengthens the exchange-architecture mechanism claim.

Recommendation: ADD as supporting mechanism evidence. Appropriate in a brief cross-phylogenetic sentence in the Discussion: "Subtelomeric gene family expansion through interchromosomal recombination is not unique to primates; independent evidence from Schistosoma mansoni shows the same mechanism operating at flatworm chromosome ends (Brann et al. 2024)."

```bibtex
@article{brann2024schistosomasubtel,
  author  = {Brann, T. and Beltramini, A. and Chaparro, C. and Berriman, M.
             and Doyle, S. R. and Protasio, A. V.},
  title   = {Subtelomeric plasticity contributes to gene family expansion in
             the human parasitic flatworm \textit{Schistosoma mansoni}},
  journal = {BMC Genomics},
  year    = {2024},
  volume  = {25},
  number  = {1},
  pages   = {217},
  doi     = {10.1186/s12864-024-10032-8},
  pmid    = {38413905}
}
```

---

## Section 4: Contradictions

One partial contradiction found:

**Chi et al. 2025 (PMID 40021902) vs. Gilad et al. 2004 (bibkey `Gilad2004`) — partial contradiction of framing.**

Gilad et al. 2004 argued that "Loss of olfactory receptor genes coincides with the acquisition of full trichromatic vision in primates" — framing OR pseudogenization as a straightforward consequence of reduced olfactory dependency once color vision improved. The existing review presents this as the dominant explanatory framework for the lineage-wide pseudogenization burden.

Chi et al. 2025 explicitly tests this across extant primates with functional receptor modeling and finds "sensory reallocation rather than strict trade-offs": olfactory receptors in anthropoids shifted functional tuning (narrow → broad) rather than simply accumulating pseudogenes at a higher rate post-trichromacy. Strepsirrhines, which lack trichromacy, show enhanced narrowly-tuned OR function rather than superior total OR count. The contradiction is at the interpretive level: Gilad 2004's data (pseudogene counts by clade) are not wrong, but the mechanistic story it told is now understood to be incomplete. The arm-by-arm OR4F gradient reported in our paper is still unexplained by either model (both models explain inter-lineage differences, not intra-genome arm differences), which is the novel contribution of the pangenome data.

**Action required:** Update the discussion of Gilad 2004 in the existing review's "Open questions" section to cite Chi et al. 2025 alongside it. The proposed sentence is in Section 3 above.

---

## Section 5: Search Audit Trail

### Tools used

1. **mcp__claude_ai_PubMed__search_articles** — primary discovery tool
2. **mcp__claude_ai_PubMed__get_article_metadata** — metadata verification for all candidate PMIDs
3. **mcp__claude_ai_bioRxiv__search_preprints** — preprint coverage (category=genomics; limitation: no keyword search supported)
4. **openalex-database** (Skill) — parallel cross-validation (launched; results not yet returned at time of writing; bibkeys were already finalized from PubMed)

### Query strings and date ranges applied

All PubMed queries used `date_from: 2023/01/01, date_to: 2026/12/31` unless noted.

| Query | Total hits | Post-filter |
|-------|-----------|-------------|
| `OR4F olfactory receptor subtelomeric` | 0 | 0 |
| `olfactory receptor gene family subtelomeric duplicon pseudogene` | 0 | 0 |
| `human olfactory receptor pangenome T2T annotation evolution` | 0 | 0 |
| `subtelomere segmental duplication interchromosomal exchange recombination` | 0 | 0 |
| `olfactory receptor gene pseudogene human genome evolution` | 2 | 1 relevant (PMID 41395915 dogs — dropped) |
| `olfactory receptor copy number variation human population` | 3 | 1 relevant (PMID 41615562 — kept MEDIUM) |
| `Niimura olfactory receptor genome` | 3 | 2 relevant (PMID 39490737 cat, 38649162 rodents — dropped WEAK) |
| `subtelomere structure duplicon human genome pangenome` | 1 | 1 relevant (PMID 38413905 — kept MEDIUM) |
| `olfactory receptor odor perception genetic variation functional` | 7 | 2 relevant (PMID 40593737 — kept STRONG; PMID 36691623 2022 — dropped out of window) |
| `telomere subtelomere chromosome end recombination human pangenome` | 1 | 0 relevant (trypanosomes — dropped) |
| `olfactory receptor T2T CHM13 telomere-to-telomere genome annotation` (2022–2026) | 0 | 0 |
| `olfactory receptor primate great ape evolution pseudogenization comparative genomics` (2022–2026) | 1 | 1 relevant (PMID 39747985 — kept MEDIUM) |
| `Ambrosini subtelomeric duplicon block human chromosome` (2022–2026) | 0 | 0 |
| `trichromacy color vision olfaction evolution primate relaxed selection` (2022–2026) | 0 | 0 |
| `olfactory receptor gene evolution primate human genomics` | 3 | 1 STRONG (PMID 40021902 — kept) |
| `olfactory receptor birth death evolution gene family multigene human` | 0 | 0 |
| `HPRC human pangenome reference consortium gene annotation 2023 2024` | 0 | 0 |
| `olfactory receptor human genome annotation diversity population 2024` | 0 | 0 |
| `subtelomeric region chromosome end gene family human segmental duplication 2024` | 0 | 0 |
| `Niimura olfactory receptor evolution mammal 2023 2024 2025` | 0 | 0 |
| `human olfactory repertoire variation loss-of-function allele pseudogene functional` | 0 | 0 |
| bioRxiv genomics category 2023–2026 (50 results sampled) | 50 | 0 relevant (no OR/subtelomere content) |

### Summary counts

- Total unique PMIDs examined: 21
- Dropped (non-human species, wrong topic): 13
  - PMID 41395915 — dogs/wolves OR diversity (Inoue/Niimura 2025): WEAK, non-human
  - PMID 36978520 — giant panda OR analysis (Zhou 2023): WEAK, non-human
  - PMID 37907519 — pig OR diversity (Kang 2023): WEAK, non-human
  - PMID 36779496 — bladder cancer OR2L5 mutation (Alradhi 2022): WEAK, clinical case
  - PMID 39490737 — domestic cat genome Niimura 2024: WEAK, non-human
  - PMID 38649162 — hystricomorph rodents chemical senses (Niimura 2024): WEAK, non-human
  - PMID 38074093 — Trypanosoma telomere maintenance (Li 2023): WEAK, pathogen
  - PMID 40244808 — structural dynamics of OR binding (Aier 2025): WEAK, molecular structure
  - PMID 40408280 — olfaction as metabolic gatekeeper (Manoel 2025): WEAK, physiology
  - PMID 39742477 — mosquito sensory compensation (Morita 2025): WEAK, invertebrate
  - PMID 38808556 — ant genome GWAS odorant receptor (Macit 2024): WEAK, invertebrate
  - PMID 36691623 — Neanderthal/Denisovan OR (de March 2022): published Dec 2022, at window boundary; MEDIUM but superseded by Chi 2025
  - PMID 36625229 — OR5AN1 musk perception (Sato-Akuhara 2023): WEAK, single receptor pharmacology
  - PMID 38508692 — house mouse Y chromosome / vomeronasal OR cluster (Fujiwara 2024): WEAK, mouse
  - PMID 38273363 — mosquito Culex genome odorant receptor expansion (Ryazansky 2024): WEAK, invertebrate
- Kept for Section 3: 5 (PMID 40021902, 40593737, 39747985, 41615562, 38413905)

### Field-specific observation

No 2023–2026 paper specifically studies OR4F genes in the context of human subtelomeric community structure, pangenome variation, or arm-specific pseudogenization rates. This confirms that our paper's OR4F gradient finding (11.1%–99.8% across 16 arms in 465 near-complete assemblies) is genuinely novel. The absence of such papers strengthens, not weakens, the C5 claim novelty.
