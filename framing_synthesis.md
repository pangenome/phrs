# Framing synthesis: Genome-wide pseudo-homologous regions

*Prepared for Erik Garrison and Andrea Guarracino*
*For submission to Nature alongside the HPRCv2 data paper*

---

## 1. The three-paper arc

This paper is the third in a series establishing pseudo-homologous regions (PHRs) as a fundamental feature of human genome organization.

**Paper 1** (Guarracino et al. 2023, Nature 617:335-343): "Recombination between heterologous human acrocentric chromosomes." Used HPRCy1 (47 individuals) to discover PHRs on the short arms of acrocentric chromosomes (chr13, 14, 15, 21, 22). Defined PHRs via positional homology entropy — diversity in reference-relative phylogenies across pangenome contigs. Found megabase-scale PHRs (chr13: 4.5 Mb, chr14: 6.5 Mb, chr15: 719 kb, chr21: 3.8 Mb, chr22: 2.8 Mb). Showed faster LD decay in PHRs than flanking arms, indicating elevated recombination. Predicted that Robertsonian translocations break within PHRs at the SST1 macrosatellite, where chr14's inverted orientation relative to chr13/chr21 enables crossover-type fusion. Published alongside the HPRC draft pangenome paper (Liao et al. 2023).

**Paper 2** (de Lima, Guarracino et al. 2025, Nature 647:952-961): "The formation and propagation of human Robertsonian chromosomes." Generated T2T assemblies of three Robertsonian chromosomes (rob(14;21), two rob(13;14)). Confirmed all three break within SST1, exactly in the predicted PHR. The chr14 SST1 inversion — unique to humans among great apes — is the structural basis. PRDM9 motif enrichment makes SST1 permissive to meiotic recombination. Accompanied by a News & Views: "A common chromosome fusion in humans explained."

**Paper 3** (this work): Generalizes PHRs genome-wide using HPRCv2 (232 individuals, 465 near-complete assemblies). Shows that 41 of 48 chromosome arms share subtelomeric sequence at >=95% identity with non-homologous chromosomes, organized into 15 communities of preferential exchange. PARs (pseudoautosomal regions) are two of those 15 communities. Validates community structure with 3D nuclear organization data across 4 technologies (Hi-C, Pore-C, CiFi, Dip-C), 3 cell types (lymphoblastoid, retinal epithelium, sperm), and in mouse meiosis.

The trajectory — discover, prove, generalize — is clean and Nature editors value this kind of arc.

---

## 2. The core thesis

**PARs are the special case. PHRs are the general case.**

Every geneticist knows pseudoautosomal regions: two non-homologous chromosomes (X and Y) share sequence at their tips, recombine there, and this must be accounted for in every genomics pipeline. The acrocentric papers showed the same thing happens between five non-homologous autosomes. This paper shows it happens across 41 of 48 chromosome arms, organized into 15 communities — and PARs are literally just two of those 15 (C14 = PAR2, C15 = PAR1).

The conceptual move: PHRs subsume PARs. The PAR concept, treated as a sex-chromosome oddity since 1986, turns out to describe a genome-wide organizational principle.

---

## 3. The PAR analogy in detail

PARs are the perfect frame because they are universally understood, have practical consequences, and are structurally analogous to what this paper describes.

### What makes PARs a "named entity"

- **Obligate meiotic crossover**: PAR1 is where X and Y must cross over; without it, X-Y segregation fails and the result is male sterility.
- **Extreme recombination rate**: ~20-fold higher than the genome average in PAR1.
- **Escape from X-inactivation**: All characterized PAR1 genes escape inactivation.
- **Clinical significance**: SHOX haploinsufficiency causes Leri-Weill dyschondrosteosis.
- **Pipeline-breaking behavior**: Without Y-masking, variant calling in PARs completely fails — Webster et al. (2019) showed 0 variants called in PAR1 before masking vs. 7,563 after.

### PARs vs PHRs: the comparison table

| Feature | PARs (X/Y) | PHRs (genome-wide) |
|---|---|---|
| Chromosomes involved | 2 (X, Y) | Up to 22 per community; 41 arms total |
| Exchange groups | 2 (PAR1, PAR2) | 15 communities |
| Size range | PAR1 ~2.4 Mb; PAR2 ~334 kb | 5-500 kb, median 105 kb |
| Sequence identity | ~100% (X/Y identical) | >=95%, peak ~98% |
| Recombination | Obligate crossover | Ectopic, not obligate |
| Named entity since | 1986 (Burgoyne) | Guarracino et al. 2023 (acrocentrics) |
| Pipeline handling | Y-masking, ploidy switching | No standard approach exists |
| Clinical exemplar | SHOX/Turner/Leri-Weill | FSHD (D4Z4), Robertsonian translocations |
| 3D proximity | X/Y pair during meiosis | Same-community arms cluster in Hi-C, Dip-C, sperm, mouse |

The punchline: PARs have had 40 years of pipeline engineering. PHRs — involving more chromosomes, more communities, and comparable sequence sizes — have had none.

### Has anyone drawn this analogy before?

No. The connection exists implicitly — Linardopoulou et al. (2005) noted PAR homology in their subtelomeric paralogy map, and PARs are technically subtelomeric — but no paper has framed subtelomeric inter-chromosomal sharing as a generalization of pseudoautosomal behavior. The fact that PARs fall out as two of fifteen communities in an unbiased clustering analysis makes the analogy structurally precise, not just metaphorical.

---

## 4. The Gerton model: why this is more than annotation

Jennifer Gerton published a working model (J Cell Sci 2024) laying out three conditions for Robertsonian chromosome formation:

1. **Sequence homology** on non-homologous chromosomes
2. **Recombination initiation** during meiosis
3. **Physical proximity** of the homologous sequences in 3D space

This paper validates all three conditions genome-wide:

1. PHR communities = sequence homology (15 communities, >=95% identity, 15,668 PHR sequences)
2. Cross-arm affinity + discordance = evidence of recombination (11% of sequences resemble a foreign arm; up to 47% of individuals heterozygous for subtelomeric type)
3. Hi-C/Dip-C/Pore-C/CiFi = physical proximity (B/W ratios 0.027-0.074; Mantel rho up to -0.66; sperm W/B = 0.401; mouse meiotic validation across 4 stages)

The implication: the conditions that cause the most common chromosomal rearrangement in humans (Robertsonian translocations, ~1 in 1,000 births) exist at nearly every chromosome end. The Gerton model, proposed for one class of rearrangement at one class of chromosome, describes a general principle.

---

## 5. What makes this Nature-level

### The four criteria for Nature genome structure papers

Based on analysis of how TADs, A/B compartments, segmental duplications, the T2T genome, and the HPRC pangenome were introduced:

1. **A new organizational principle** — not just "subtelomeres are repetitive" but "subtelomeres form a structured community system of inter-chromosomal pseudo-homology"
2. **Universality** — 41 of 48 arms, 4 technologies, 3 cell types, 2 species, 232 individuals
3. **Functional consequences** — FSHD, Robertsonian translocations, and the uncharted clinical territory of 13 other communities
4. **A simple framework readers carry away** — "PARs are the special case; PHRs are the general case"

### What distinguishes this from adjacent work

| Prior work | What they had | What this paper adds |
|---|---|---|
| Mefford & Trask 2002 (NRG) | Qualitative "patchwork" model, single genome | Population-scale (232 individuals), community structure |
| Linardopoulou et al. 2005 (Nature) | "Hot spots" label, single genome | 15 discrete communities, not just "hot" |
| Ambrosini et al. 2007 | 11 duplicon blocks, bimodal identity | All blocks mapped to communities, population validation |
| Young/Stong 2020 (PLoS Genet) | Optical mapping, 154 genomes | Sequence resolution, pangenome graph, community detection |
| Guarracino et al. 2023 (Nature) | PHRs at 5 acrocentric arms | PHRs at 41 arms genome-wide |
| Vollger et al. 2024 (Nat Genet) | SD population genetics, 170 samples | Community structure + 3D validation |
| Yardimci et al. 2024 (Nat Comms) | 40K inter-chromosomal contacts, 62 Hi-C datasets | Connected to sequence-sharing communities |

Nobody has combined population-scale pangenome sequence analysis with 3D nuclear organization validation to show that inter-chromosomal sharing forms a structured community system. Vollger had the population scale but not the communities or 3D. Yardimci had the 3D but not the sequence communities.

---

## 6. Practical consequences for genomics

The "why should I care" for the genomics audience. If you are doing anything with short reads near chromosome ends, you are already affected by PHRs and do not know it.

### Read mapping

Reads from PHR regions map ambiguously to multiple non-homologous chromosomes. A read from chr21_p could map to chr13_p, chr14_p, chr15_p, or chr22_p (community C7). Unlike PARs, where only 2 chromosomes compete, PHR communities involve up to 22 chromosomes. MapQ scores will be near-zero; reads will scatter across chromosomes; coverage will appear artificially low on each individual chromosome.

### Variant calling

Paralogous sequence variants (PSVs) will be miscalled as SNPs. True variants will be filtered for low mapQ. Without PHR annotation, variant callers produce a mixture of false positives (PSVs called as SNPs), false negatives (true variants in mapQ=0 regions), and chromosome misassignment.

### Structural variant calling

Inter-chromosomal SV callers will detect split reads and discordant pairs spanning PHR boundaries between non-homologous chromosomes. These will be called as translocations or fusions. In reality, many are normal subtelomeric sharing. SV callsets near chromosome ends will be contaminated with false-positive translocations. This is especially acute for cancer genomics, where inter-chromosomal rearrangements are clinically actionable.

### Phasing

PHRs create phase deserts at chromosome tips, analogous to centromeric phase breaks. Reads that cannot be assigned to a single chromosome cannot contribute to haplotype phase. Hi-C-based phasing will be contaminated by cross-chromosome contacts.

### FSHD diagnostics

The D4Z4 locus (community C1, chr4q/chr10q) is the paradigmatic example. The subtelomeric regions share ~99% identity over ~150 kb. Translocations between these regions are common in the general population (~20% carry translocated arrays). Standard assays cannot reliably distinguish chromosome of origin. Ricci et al. (2024, Science Advances) showed that 4q-10q translocations can revert the FSHD phenotype, but detecting them requires subtelomere-aware analysis. A 2025 EJHG paper ("Rethinking genomics of FSHD in the telomere-to-telomere era") showed D4Z4-like repeats on at least ten additional chromosomes beyond 4q and 10q — the inter-chromosomal sharing network is far more extensive than previously appreciated.

### The pipeline prescription

For PARs, the community converged on: mask PARs on Y, map to sex-informed reference, call variants separately with appropriate ploidy, annotate PAR boundaries explicitly. The PHR equivalent: annotate all PHR boundaries genome-wide; use pangenome-graph-based mapping; treat PHR communities as explicit units in variant calling; flag subtelomeric variants with their community membership. The PHR community catalog this paper provides is a necessary first step.

---

## 7. How landmark papers framed similar discoveries

### The "missing 8%" pattern (T2T)

The T2T papers (Nurk et al. 2022, Science) framed findings as revealing what was hidden: 200 million new base pairs in the "remaining 8%" of the genome. The message: completion changes everything quantitatively. This paper has a version: these regions were unmappable before pangenome-quality assemblies, and now the full extent of inter-chromosomal sharing is visible for the first time.

### The "crucibles" framing (Eichler)

Evan Eichler's 2006 NRG review titled segmental duplications as "crucibles of evolution, diversity and disease" — a single word that reframed assembly artifacts as engines of innovation. PHR communities are subtelomeric crucibles: they generate both evolutionary novelty (gene family diversification, olfactory receptor variation) and pathology (FSHD, Robertsonian translocations, subtelomeric deletion syndromes).

### What makes a genomic term stick

Based on analysis of TADs, A/B compartments, segmental duplications, PARs, and the fractal globule:

1. **Descriptive transparency** — the name tells you what the thing is. "Pseudo-homologous region" = a region that is pseudo-homologous. This works.
2. **Pronounceable acronym** — TADs, PARs, SDs, PHRs. All monosyllabic or disyllabic.
3. **Fills a conceptual gap** — before "PHRs," people said "subtelomeric regions that share sequence between chromosomes." A mouthful.
4. **Introduced in a high-impact venue** — every one of these terms was introduced in Nature or Science.
5. **Has a visual/structural referent** — TADs have the Hi-C triangle pattern. PHR communities have network/clustering visualizations.

The term "pseudo-homologous regions" is well-established from the 2023 paper. The new concept to name is the community structure itself. "PHR communities" or "exchange communities" are natural descriptors.

### The ENCODE cautionary tale

ENCODE's 2012 claim that 80% of the genome was "functional" was the most famous overreach of the dark-matter reframing. The problem: conflating biochemical activity with biological function. This paper should frame PHRs as an organizational principle with structural and clinical consequences — not claim the shared sequences themselves are functional. The exchange system matters because it generates rearrangements, organizes nuclear architecture, and confounds genomic analysis. That is enough.

---

## 8. Recent competition and context

### Adjacent work (2024-2026)

- **Vollger et al. 2024 (Nature Genetics)**: 170 HPRC assemblies, 173.2 Mb duplicated sequence. SD population genetics. Does not identify community structure or validate in 3D.
- **Yardimci et al. 2024 (Nature Comms)**: 40,282 non-homologous chromosomal contacts across 62 Hi-C datasets. Shows inter-chromosomal contacts are real and enriched at chromosome ends. Does not connect to sequence-sharing communities.
- **Karimian et al. 2024 (Science)**: Telomere Profiling showing telomere length is chromosome end-specific. Had to develop specialized algorithms to handle subtelomeric variation — acknowledging the diversity this paper maps.
- **4D Nucleome Project 2025 (Nature)**: >140,000 looping interactions, single-cell 3D models. General framework for nuclear organization.

### Follow-up work on the acrocentric PHR story

- **Lin et al. (Dec 2025, bioRxiv)**: 156 phased acrocentric short arms, 10x elevated de novo SNV rate, one ectopic chr13-chr21 recombination breakpoint in a segmental duplication 1.6 Mb distal to SST1. Complementary — they see individual recombination events; this paper sees the population-scale community structure.
- **Solar, Guarracino et al. (Dec 2025, bioRxiv)**: PHRs in great ape acrocentrics across 25 million years of evolution. Centromere repositioning by whole-arm inversion.
- **Biobank-scale ROB genotyping (Mar 2026, bioRxiv)**: Types ROB carriers from short reads in UK Biobank (n=490K). 0.11-0.12% ROB frequency. Validates clinical relevance at population scale.
- **EJHG 2025**: "Rethinking genomics of FSHD in the telomere-to-telomere era." D4Z4-like repeats on at least ten additional chromosomes beyond 4q/10q.

None of these do what this paper does — genome-wide PHR mapping with 3D validation. But they create context that makes the timing ideal.

---

## 9. Suggested paper structure

For Nature format (~3,000 words main text + methods + extended data):

### Opening paragraph

> Pseudoautosomal regions, where non-homologous sex chromosomes share sequence and undergo obligate meiotic crossover, are a foundational concept in human genetics with direct consequences for variant calling, inheritance modeling, and clinical interpretation. We recently showed that analogous pseudo-homologous regions (PHRs) on the short arms of acrocentric chromosomes mediate recombination between non-homologous autosomes, predicting the breakpoints of Robertsonian translocations — a prediction subsequently confirmed by telomere-to-telomere assembly of three Robertsonian chromosomes. Here, using 465 near-complete assemblies from the Human Pangenome Reference Consortium and CHM13, we show that PHRs are a genome-wide feature of human subtelomeres. Forty-one of 48 chromosome arms share sequence at >=95% identity with non-homologous chromosomes, organized into 15 communities of preferential exchange — of which PARs are two. These communities are reflected in 3D nuclear co-localization across four proximity technologies, three cell types, and in mouse meiosis, establishing subtelomeric pseudo-homology as a pervasive system of inter-chromosomal organization with implications for genome stability, structural variant interpretation, and disease.

### Figure plan

**Figure 1: PHRs are genome-wide.** Community map showing all 15 communities. Chord diagram or network graph of the 41 arms. PARs (C14, C15) highlighted as the familiar reference point. The PAR/PHR comparison table as an inset or panel.

**Figure 2: Population-scale polymorphism.** Cross-arm affinity rates across arms. f7501 population distribution (reproducing Mefford & Trask's Fig 3 at 233x scale). Discordance rates. The message: this is a segregating polymorphism, not a frozen relic.

**Figure 3: 3D nuclear organization mirrors sequence communities.** Convergent evidence panel: Hi-C (multiple samples), Dip-C (single-cell), sperm, mouse meiosis. The flanking paradox as an inset (unique sequence shows stronger signal — rules out artifact). Per-arm-pair correlation (Jaccard vs contact).

**Figure 4: Cross-species conservation and meiotic context.** Mouse meiotic Hi-C across 4 stages (leptotene to diplotene). Zygotene peak consistent with bouquet-stage telomere clustering. Community structure in mouse (different architecture — TLC/L1-LINE driven — but same principle).

**Extended Data**: The exhaustive per-community tables, multi-resolution sensitivity, RPE-1 self-discovery (including t(X;10) translocation detection as positive control), sequence-level 50-community partition, gene enrichment, TAR1 analysis, two-domain model.

### Discussion points

1. **The Gerton model generalized**: The three conditions for Robertsonian translocation — sequence homology, recombination competence, 3D proximity — are met genome-wide.
2. **The feedback loop**: Sequence similarity leads to 3D proximity leads to ectopic exchange leads to increased similarity. Causal direction cannot be established from present data, but convergent evidence across technologies and species supports the model.
3. **The meiotic bouquet**: Zuo et al. (2021) showed chromosome-end alignment extends ~20% of chromosome length during meiotic prophase. Median PHR (105 kb) fits within a single meiotic loop (~500 kb at leptotene). Mouse meiotic Hi-C peaks at zygotene, consistent with bouquet-stage proximity enabling exchange.
4. **Practical consequences**: The PHR community catalog as a resource. The analogy to PAR annotations.
5. **What this paper does not show**: Exchange timing, causal direction, whether specific communities will produce clinically relevant rearrangements beyond the known cases (FSHD, Robertsonian translocations).

---

## 10. Title options

1. **"Pseudo-homologous regions across the human genome"** — clean, direct, echoes prior titles
2. **"Genome-wide pseudo-homology at human chromosome ends"**
3. **"Communities of inter-chromosomal exchange at human subtelomeres"**
4. **"Subtelomeric pseudo-homologous regions organize inter-chromosomal exchange across the human genome"**

Recommendation: option 1 for brevity and continuity with Guarracino et al. 2023. It says what it is and signals the genome-wide extension.

---

## 11. The one sentence for the cover letter

> We show that pseudo-homologous regions — previously described only at acrocentric chromosomes — are a genome-wide system of inter-chromosomal sequence sharing at human subtelomeres, organized into 15 exchange communities that predict 3D nuclear co-localization across cell types and species, with implications for variant calling, structural variant interpretation, and disease genetics at chromosome ends.

---

## 12. Key references

### The three-paper arc
- Guarracino et al. 2023. "Recombination between heterologous human acrocentric chromosomes." Nature 617:335-343. doi:10.1038/s41586-023-05976-y
- de Lima, Guarracino et al. 2025. "The formation and propagation of human Robertsonian chromosomes." Nature 647:952-961. doi:10.1038/s41586-025-09540-8
- Gerton 2024. "A working model for the formation of Robertsonian chromosomes." J Cell Sci 137:jcs261912.

### Foundational subtelomere literature
- Mefford & Trask 2002. "The complex structure and dynamic evolution of human subtelomeres." Nat Rev Genet 3:91-102.
- Linardopoulou et al. 2005. "Human subtelomeres are hot spots of interchromosomal recombination and segmental duplication." Nature 437:94-100.
- Ambrosini et al. 2007. "Human subtelomeric duplicon structure and organization." Genome Biol 8:R151.
- Flint et al. 1997. "The relationship between chromosome structure and function at a human telomeric region." Nat Genet 15:252-257.

### PAR references
- Mangs & Morris 2007. "The human pseudoautosomal region (PAR): origin, function and future." Curr Genomics 8:129-136.
- Webster et al. 2019. "Identifying, understanding, and correcting technical artifacts on the sex chromosomes in next-generation sequencing data." GigaScience 8:giz074.

### Adjacent pangenome/3D work
- Liao et al. 2023. "A draft human pangenome reference." Nature 617:312-324.
- Vollger et al. 2024. "Structural polymorphism and diversity of human segmental duplications." Nat Genet.
- Yardimci et al. 2024. "Inter-chromosomal contacts demarcate genome topology along a spatial gradient." Nat Comms.
- Zuo et al. 2021. "Meiotic chromosome organization and its role in recombination and cancer." (Mouse meiotic Hi-C)

### Clinical relevance
- Ricci et al. 2024. Science Advances. (4q/10q translocation reverts FSHD)
- EJHG 2025. "Rethinking genomics of FSHD in the telomere-to-telomere era."
- Biobank-scale ROB genotyping 2026. bioRxiv. (0.11% ROB frequency in UK Biobank)

### Concept framing precedents
- Eichler 2006. "Primate segmental duplications: crucibles of evolution, diversity and disease." Nat Rev Genet.
- Dixon et al. 2012. "Topological domains in mammalian genomes identified by analysis of chromatin interactions." Nature 485:376-380.
- Lieberman-Aiden et al. 2009. "Comprehensive mapping of long-range interactions reveals folding principles of the human genome." Science 326:289-293.
