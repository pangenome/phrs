# Literature Refresh: Topic 10 — Pangenome Graphs and IMPG

**Task:** lit-refresh-10-pangenome  
**Date:** 2026-05-17  
**Reviewer:** agent-56 (Documenter role)

---

## Section 1: Topic scope

Topic 10 covers the methodological lineage of the pangenome graph tools — wfmash, seqwish, pggb, odgi, and IMPG — that together constitute the computational pipeline for this subtelomere study. In the Nature draft this topic appears in three places: (1) the abstract phrase "implicit pangenome graph — a reference-free, all-to-all alignment sampling roughly 12% of pairwise haplotype combinations without chromosomal partitioning"; (2) paragraph P3 (Methods: "We treat every haplotype as its own reference…") which cites `@pangenome_graphs_impg_GarrisonGuarracino2023`, `@pangenome_graphs_impg_GuarracinoHeumos2022`, `@pangenome_graphs_impg_IMPG2023`, `@pangenome_graphs_impg_Hickey2024`, `@Garrison2018`, and `@Garrison2024pggb`; and (3) the Methods section (L50–L80) which specifies all tool versions and parameters. The key claims anchored here are C2 (implicit pangenome graph design), the 12% pairwise sampling rate with Erdős-Rényi connectivity justification, the IMPG transitive-closure operation as the PHR extraction mechanism, and the pggb+odgi graph for Jaccard similarity. The CONSISTENCY_AUDIT_v1.md flags two NOT-IN-REPORT items for this topic: (i) the 11.6%/12% sampling rate claim has no external citation and is a documented methods gap (rows 13, 97); (ii) IMPG is cited only as a GitHub @misc, not a journal article — a weak citation for the core method.

---

## Section 2: Existing citations still authoritative

| Bibkey | Status | Reason |
|---|---|---|
| `pangenome_graphs_impg_Lee2002` | STILL CURRENT | POA (2002) remains the conceptual ancestor of smoothxg; no superseding treatment. |
| `pangenome_graphs_impg_Garrison2018` | STILL CURRENT | vg is the foundational variation graph toolkit; extended by Chang 2025 (preprint) but not superseded. |
| `pangenome_graphs_impg_Jain2018` | STILL CURRENT | MashMap provides the sketching basis for wfmash; no published wfmash-specific paper exists; this remains the only citable source for the prefiltering design. |
| `pangenome_graphs_impg_Eizenga2020` | STILL CURRENT | Best comprehensive review of pangenome graph concepts; no post-2023 review supersedes it. |
| `pangenome_graphs_impg_Li2020` | STILL CURRENT | Minigraph remains in active use; rGFA format unchanged; contrast with pggb intact. |
| `pangenome_graphs_impg_GuarracinoHeumos2022` | STILL CURRENT | odgi is in active use; Heumos2024 (SGD layout) and Heumos2024nfcore (nf-core pipeline) extend it without replacing it. The `odgi similarity --all -P` command used in this study comes from this paper. |
| `pangenome_graphs_impg_GarrisonGuarracino2023` | STILL CURRENT | seqwish graph induction algorithm is unchanged; still the graph-induction step in pggb. |
| `pangenome_graphs_impg_Nurk2022` | STILL CURRENT | T2T-CHM13 is the reference used for subtelomere extraction; not superseded. |
| `pangenome_graphs_impg_Hickey2024` | STILL CURRENT | Minigraph-Cactus is the reference-anchored foil for pggb; the contrast is the methodological argument for why no chromosomal partitioning is used. |
| `pangenome_graphs_impg_Liao2023` | STILL CURRENT | HPRC v1 is the direct predecessor dataset; HPRC v2 (hprc_hprcv2_2025) extends rather than replaces it. |
| `pangenome_graphs_impg_Garrison2024pggb` | STILL CURRENT | Formal pggb description in Nature Methods (2024); primary citation for the pipeline. |
| `pangenome_graphs_impg_Guarracino2023` | STILL CURRENT | Acrocentric study is the direct methodological precedent; its pipeline is the same one used here. |
| `pangenome_graphs_impg_IMPG2023` | STILL CURRENT but FLAG | IMPG is cited as a GitHub @misc. It remains the only citable entry for the tool. No journal paper for IMPG has been published as of this review. **Recommendation:** flag this in methods as a tool reference; consider a data-in-brief or methods note to convert it to a citable publication before submission. |
| `pangenome_graphs_impg_Siren2021` | STILL CURRENT | Giraffe pangenome aligner at 5,202 genomes; the Chang 2025 long-read extension is a preprint, not yet peer-reviewed. |
| `pangenome_graphs_impg_Siren2022` | STILL CURRENT | GBZ format for coordinate-preserving pangenome graphs; the coordinate-free implicit approach in this study is explicitly contrasted with GBZ. |

---

## Section 3: NEW papers to add (STRONG relevance, 2023+)

---

**[andreace2023pangenome]** Andreace et al., 2023. Comparing methods for constructing and representing human pangenome graphs. *Genome Biology* 24:274. PMID:38037131 / DOI:10.1186/s13059-023-03098-2

Claim it supports: Draft Methods P3, and topic-10 review section "Reference-free, without chromosomal partitioning": "Reference-seeded graph builders such as Minigraph [@Li2020] and Minigraph-Cactus [@Hickey2024] inherit the chromosomal coordinates of the reference genome… This is not a bug but a design choice: reference-anchored graphs prioritize coordinate stability over completeness." Andreace et al. directly compare pggb, Minigraph-Cactus, Minigraph, Bifrost, and mdbg on 52 human haplotypes plus CHM13 and GRCh38, showing that pggb captures variation at all scales (SNPs to megabase SVs) that reference-anchored methods miss in highly variable regions. This is the first independent controlled benchmark establishing pggb's completeness advantage and providing empirical support for the "reference-free" design choice.

Recommendation: ADD alongside `@pangenome_graphs_impg_Hickey2024` in Methods P3 wherever pggb is contrasted with reference-anchored methods. Also add to the topic-10 review section "Reference-free, without chromosomal partitioning."

```bibtex
@article{andreace2023pangenome,
  author    = {Andreace, Francesco and Lechat, Pierre and Dufresne, Yoann and Chikhi, Rayan},
  title     = {Comparing methods for constructing and representing human pangenome graphs},
  journal   = {Genome Biology},
  year      = {2023},
  volume    = {24},
  number    = {1},
  pages     = {274},
  doi       = {10.1186/s13059-023-03098-2},
  pmid      = {38037131},
  note      = {Independent controlled benchmark of pggb, Minigraph-Cactus, Minigraph, Bifrost, and mdbg on 52 human haplotypes; shows pggb captures SNP-to-SV variation missed by reference-anchored methods; empirical support for reference-free design choice.}
}
```

---

**[heumos2024nfcore]** Heumos et al., 2024. Cluster-efficient pangenome graph construction with nf-core/pangenome. *Bioinformatics* 40(11):btae609. PMID:39400346 / DOI:10.1093/bioinformatics/btae609

Claim it supports: End-to-end-report §01 Pipeline ("pggb -p 95 -D /scratch") and NATURE_DRAFT_v1.md Methods L56 ("pggb -p 95 -D /scratch"). The draft implicitly claims that "the current subtelomere study is among the first publications to demonstrate the pangenome graph methods stack at 466-haplotype scale" (topic-10 review, Open questions section). Heumos et al. describe nf-core/pangenome, a production-grade Nextflow wrapper for the exact pipeline (wfmash → seqwish → smoothxg → gfaffix → odgi) deployed on HPC clusters. The paper demonstrates construction of a 1000-haplotype chromosome-19 pangenome (2-3x the subtelomere dataset size) with a 2-3× speedup over standalone pggb. This directly validates the scalability of the pipeline at larger-than-current-study scale.

Recommendation: ADD to Methods section after the `@Garrison2024pggb` citation where pggb is described, to provide the production-deployment and scalability reference. Also cite in the "HPRC v2 and population-scale implicit graphs" open-questions paragraph.

```bibtex
@article{heumos2024nfcore,
  author    = {Heumos, Simon and Heuer, Michael L. and Hanssen, Friederike and Heumos, Lukas and
               Guarracino, Andrea and Heringer, Peter and Ehmele, Philipp and Prins, Pjotr and
               Garrison, Erik and Nahnsen, Sven},
  title     = {Cluster-efficient pangenome graph construction with nf-core/pangenome},
  journal   = {Bioinformatics},
  year      = {2024},
  volume    = {40},
  number    = {11},
  pages     = {btae609},
  doi       = {10.1093/bioinformatics/btae609},
  pmid      = {39400346},
  note      = {nf-core/pangenome: Nextflow HPC wrapper for the wfmash→seqwish→smoothxg→gfaffix→odgi pipeline; constructs 1000-haplotype chromosome-19 pangenome; 2-3x speedup over standalone pggb; validates scalability beyond 466-haplotype current study scale.}
}
```

---

**[leonard2023bovine]** Leonard et al., 2023. Graph construction method impacts variation representation and analyses in a bovine super-pangenome. *Genome Biology* 24:124. PMID:37217946 / DOI:10.1186/s13059-023-02969-y

Claim it supports: Topic-10 review section "What is contested" and "Reference-free, without chromosomal partitioning." Leonard et al. build super-pangenomes using pggb, Cactus, and Minigraph with taurine/indicine cattle and related bovids, finding that pggb and Cactus achieve ~95% exact matches with assembly-derived small variant calls while Minigraph (reference-seeded, no base-level variation) misses these. This independently demonstrates, in a non-human context, that reference-seeded methods underrepresent variation relative to all-vs-all approaches. It also validates the use of pggb specifically for downstream Jaccard similarity computation (base-level variation is required for meaningful graph node sharing).

Recommendation: ADD to Methods section alongside `@Garrison2024pggb` as independent validation that pggb's all-vs-all approach captures base-level variation critical for Jaccard similarity. MEDIUM-STRONG relevance.

```bibtex
@article{leonard2023bovine,
  author    = {Leonard, Alexander S. and Crysnanto, Danang and Mapel, Xena M. and Bhati, Meenu and Pausch, Hubert},
  title     = {Graph construction method impacts variation representation and analyses in a bovine super-pangenome},
  journal   = {Genome Biology},
  year      = {2023},
  volume    = {24},
  number    = {1},
  pages     = {124},
  doi       = {10.1186/s13059-023-02969-y},
  pmid      = {37217946},
  note      = {Multi-species bovine super-pangenome comparing pggb, Cactus, Minigraph; pggb and Cactus achieve ~95% match with assembly-derived small variants, Minigraph misses base-level variation; independent non-human validation that reference-seeded methods underrepresent variation relative to all-vs-all pggb.}
}
```

---

**[kaushan2026tracepoints]** Kaushan et al., 2026. Adaptive Tracepoints for Pangenome Alignment Compression. *bioRxiv*. PMID:41757015 / DOI:10.64898/2026.02.16.706236

Claim it supports: Topic-10 review Open questions section: "IMPG's transitive closure is efficient for the current dataset (18,827 sequences, ~88 GB PAF) but the depth of the closure…"; also "compact CIGAR-delta representations" in the IMPG description. By Kaushan, Marco-Sola, Garrison, Prins, and Guarracino (the IMPG development team). The paper proposes adaptive tracepoints — complexity-aware alignment encoding segmenting PAF/CIGAR data by edit distance or diagonal distance rather than fixed intervals. On real pangenomes with 390M alignments, it achieves 23-139× compression relative to uncompressed representations with no score degradation and linear-time reconstruction. This directly addresses the storage challenge for large all-vs-all PAF sets (the current study's ~88 GB PAF), and represents the next algorithmic step in the implicit graph infrastructure lineage. STRONG relevance to the "future work" discussion of IMPG scalability and efficient representation.

Recommendation: ADD to Open questions section "Scalability of transitive closure" and to the IMPG description footnote. Also cite in Methods L80 (software versions) if this encoding is planned for future pipeline versions. Note: preprint as of May 2026; verify publication status before submission.

```bibtex
@article{kaushan2026tracepoints,
  author    = {Kaushan, Hasitha and Marco-Sola, Santiago and Garrison, Erik and Prins, Pjotr and Guarracino, Andrea},
  title     = {Adaptive Tracepoints for Pangenome Alignment Compression},
  journal   = {bioRxiv},
  year      = {2026},
  doi       = {10.64898/2026.02.16.706236},
  pmid      = {41757015},
  note      = {Preprint. Complexity-aware PAF/CIGAR alignment encoding achieving 23-139x compression on 390M real pangenome alignments; directly relevant to IMPG's compact CIGAR-delta storage for ~88 GB all-vs-all PAF. From the IMPG development team (Garrison, Guarracino).}
}
```

---

**[edwards2025multispecies]** Edwards et al., 2025. Multispecies pangenomes reveal a pervasive influence of population size on structural variation. *Science* 390(6778):eadw1931. PMID:41379974 / DOI:10.1126/science.adw1931

Claim it supports: Topic-10 review Open questions "IMPG in broader population genetics contexts" and the general claim that pggb-based pangenome analysis is a mature, broadly applicable approach. Edwards et al. use 45 long-read assemblies of North American jays (55-fold range in N_e) and apply pggb-based pangenome methods (with Garrison and Guarracino as co-authors) to detect SVs shaped by population size. Published in Science, this demonstrates the pggb/pangenome tool stack at vertebrate pangenome scale outside of human genomics, providing the strongest citation for broad cross-species validation of the pipeline.

Recommendation: ADD to the Open questions section "IMPG in broader population genetics contexts" as the most recent and highest-impact demonstration of pggb-based pangenome analysis in a non-human system. MEDIUM relevance to core claims.

```bibtex
@article{edwards2025multispecies,
  author    = {Edwards, Scott V. and Fang, Bohao and Khost, Danielle and Kolyfetis, George E.
               and Cheek, Rebecca G. and DeRaad, Devon A. and Chen, Nancy and
               Fitzpatrick, John W. and McCormack, John E. and Funk, W. Chris and
               Ghalambor, Cameron K. and Garrison, Erik and Guarracino, Andrea and
               Li, Heng and Sackton, Timothy B.},
  title     = {Multispecies pangenomes reveal a pervasive influence of population size on structural variation},
  journal   = {Science},
  year      = {2025},
  volume    = {390},
  number    = {6778},
  pages     = {eadw1931},
  doi       = {10.1126/science.adw1931},
  pmid      = {41379974},
  note      = {45 long-read jay assemblies; pggb-based pangenome analysis (Garrison, Guarracino co-authors); SVs modulated by population size; Science-level validation of pggb across vertebrate diversity outside human genomics.}
}
```

---

## Section 4: CONTRADICTIONS

None found in the searched window (2023-01-01 to 2026-05-17). No paper in the retrieved set argues against the implicit graph design, the all-vs-all alignment approach, the IMPG transitive-closure operation, or the 95% identity threshold as a valid PHR cutoff. The Andreace 2023 benchmark confirms pggb's superiority for complete variation capture; Leonard 2023 independently corroborates this in bovines. The field's trajectory (nf-core/pangenome scaling to 1000 haplotypes, vg Giraffe extending to long reads) is consistent with increasing adoption of the methods stack used here, not contradiction.

**One flag (not a contradiction, but a methods gap):** The 12% wfmash sampling rate claim and the Erdős-Rényi connectivity argument remain uncited in any external publication (CONSISTENCY_AUDIT_v1.md rows 13 and 97 both flag this as NOT-IN-REPORT). No paper published 2023-2026 fills this gap. The argument is sound but needs to be formalized either in the companion paper methods section or in a brief note. No new literature weakens the argument.

---

## Section 5: Search audit trail

### Tools used

| Tool | Queries run | Notes |
|---|---|---|
| PubMed (mcp__claude_ai_PubMed__search_articles) | 15 queries | See query strings below |
| PubMed (mcp__claude_ai_PubMed__get_article_metadata) | 4 batch calls | Metadata verification of all retrieved PMIDs |
| bioRxiv (mcp__claude_ai_bioRxiv__search_preprints) | 1 query (genomics category) | No topic-specific hits; bioRxiv API does not support keyword search |
| bioRxiv (mcp__claude_ai_bioRxiv__get_preprint) | 3 DOI lookups | Three specific DOIs checked; none resolved |
| OpenAlex (openalex-database skill) | 3 queries | Cross-validation; caught Andreace 2023 missed by PubMed queries |

### Date range applied

2023-01-01 to 2026-05-17 (PubMed date_from/date_to); OpenAlex filter publication_year 2023-2026.

### PubMed query strings used

1. `pangenome graph variation Garrison Guarracino[Author]` → 9 hits
2. `Garrison E[Author] pangenome` → 51 total (20 retrieved, paged)
3. `Guarracino A[Author] pangenome genomics` → 40 total (20 retrieved, paged)
4. `Human Pangenome Reference Consortium 2025 assemblies haplotypes` → 9 hits
5. `Jain C[Author] whole genome sequence alignment pangenome` → 1 hit (mm2-plus, not relevant)
6. `Andreace F[Author] comparing methods pangenome graph construction representing human 2023` → 1 hit (PMID 38037131) ✓
7. `Leonard AS[Author] graph construction method variation representation pangenome bovine` → 1 hit (PMID 37217946) ✓
8. `Rautiainen M[Author] pangenome genome assembly graph` → 2 hits (Verkko2, not relevant)
9. `pangenome graph vg giraffe haplotype read alignment 2024` → 1 hit (PMID 39261641, not relevant)
10. Various additional queries (wfmash, IMPG implicit, MashMap3, nf-core pangenome) → 0 hits each (PubMed query translation too restrictive)

### OpenAlex query strings used

1. `pangenome graph wfmash PGGB seqwish odgi alignment` (2023-2026, sorted by citations) → ~20 hits
2. `IMPG implicit pangenome graph transitive closure interval` (2023-2026) → 0 hits
3. `Garrison pangenome graph variation unbiased reference-free` (2023-2026, sorted by citations) → ~15 hits

### Total hits before relevance filter

~95 unique papers across all searches and tools.

### Hits after relevance filter

5 papers (STRONG or MEDIUM relevance to topic-10 draft claims).

### Dropped papers (1 line each)

| PMID | Title (truncated) | Drop reason |
|---|---|---|
| 39990470 | Comparative population pangenomes (jays, preprint) | Published version is Edwards2025 above; dropped duplicate |
| 40506254 | Accurate short-read alignment through-index-based pangenome indexing | Read mapping optimization; no bearing on all-vs-all implicit graph claims |
| 40152239 | wgatools: ultrafast toolkit for manipulating whole-genome alignments | PAF/WGA processing tool; tangential to IMPG core operation |
| 39626271 | Panacus: fast pangenome growth and core size estimation | Analyses GFA from pggb; MEDIUM at best; does not directly support any draft claim |
| 40269156 | Human de novo mutation rates from a four-generation pedigree reference | Porubsky2025 pedigree paper; already cited in draft; not topic-10 |
| 40759746 | The Platinum Pedigree: a long-read benchmark | Kronenberg 2025 = Porubsky2025 in REFERENCES_v3; already cited |
| 41282249 / 41256518 | vg Giraffe long/short read mapping to large pangenome graphs | Preprint (Chang et al.); MEDIUM relevance; shows ecosystem maturation but no specific claim gap |
| 40389285 | Verkko2 integrates proximity-ligation data | Assembly tool; not relevant to graph analysis or implicit graph framing |
| 37217946 already included | | |
| 41757015 already included | | |
| 39400346 already included | | |
| 38037131 already included | | |
| 41379974 already included | | |
| All remaining PubMed/OpenAlex hits | Various | WEAK: non-human pangenome applications, reviews of unrelated biology, forensics, clinical genomics — none bearing on the draft's implicit graph methodology claims |
