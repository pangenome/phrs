## Community gene enrichment analysis

**What it does.** Tests which genes are shared across community member arms, which are community-specific, and which span multiple communities.

**How.** Gene annotations (the annotation section) grouped by arm-level community. Shared gene instances counted across arms within each community. Fisher's exact tests for per-community enrichment (116 tests, BH-corrected).

**Key metrics.** 374 unique genes across 39 arms (chr7_q and chr12_q have zero gene annotations). 576 shared gene instances across arms within communities. 93 genes community-specific (found in only one community). 216 genes in 2 or more communities.

### Biotype composition

**Result.** Subtelomeric genes are predominantly pseudogenes and ncRNA, with low protein-coding content:

| Community | Arms | Total genes | Protein-coding (%) | Pseudogene (%) | ncRNA (%) |
|-----------|------|-------------|---------------------|-----------------|-----------|
| C1 | chr4_q, chr10_q | 59 | 5.1 | 86.4 | 8.5 |
| C2 | chr10_p, chr18_p | 14 | 21.4 | 28.6 | 50.0 |
| C3 | chr3_q, chr7_p, chr9_q, chr11_p, chr16_q, chr19_p | 195 | 6.7 | 55.4 | 37.9 |
| C4 | chr7_q, chr12_q | 0 | — | — | — |
| C5 | chr6_p, chr9_p, chr12_p, chr20_q | 101 | 5.9 | 61.4 | 32.7 |
| C6 | chr1_q, chr13_q, chr17_q, chr19_q, chr21_q, chr22_q | 34 | 8.8 | 55.9 | 35.3 |
| C7 | chr13_p, chr14_p, chr15_p, chr21_p, chr22_p | 118 | 5.1 | 65.3 | 29.7 |
| C8 | chr15_q | 72 | 4.2 | 50.0 | 45.8 |
| C9 | chr16_p | 136 | 5.1 | 56.6 | 38.2 |
| C10 | chr17_p | 39 | 5.1 | 64.1 | 30.8 |
| C11 | chr1_p, chr5_q, chr6_q, chr8_p | 162 | 4.3 | 57.4 | 38.3 |
| C12 | chr2_q, chr20_p | 129 | 3.9 | 57.4 | 38.8 |
| C13 | chr4_p | 12 | 8.3 | 66.7 | 25.0 |
| C14 | chrX_q, chrY_q | 60 | 8.3 | 58.3 | 33.3 |
| C15 | chrX_p, chrY_p | 28 | 32.1 | 39.3 | 28.6 |

C15 (PAR1) has 32.1% protein-coding genes, reflecting the functional gene content of the pseudoautosomal region. C4 (chr7_q/chr12_q) has zero gene annotations — its similarity regions are confined to 5–25 kb at the telomeric tip. C2 (chr10_p/chr18_p) has 21.4% protein-coding and 50.0% ncRNA.

### Olfactory receptor genes

**What it does.** Tests whether olfactory receptor gene families are present across communities, as predicted by Mefford & Trask (2002) and Ambrosini et al. (2007).

**Key metrics.** 10 OR4F family genes detected across 7 communities (C3, C5, C8, C9, C11, C12, C14). OR4F5 and OR4F8P most widespread (14 arms each). IQSEC3 detected in C5 (chr12_p, 453 samples).

**Result.** Confirms at population scale that "human subtelomeres can contain genes, such as members of the olfactory receptor gene family" (Mefford & Trask 2002) and extends Ambrosini et al.'s (2007) OR duplicon architecture (Table 1, Block 2) to 465 haplotypes.

### Ambrosini subtelomere-specific blocks → Leiden communities

**What it does.** Maps each of Ambrosini et al.'s (2007) 11 subtelomere-specific duplicon block entries (Table 1, numbered 1–3, 5–8, 10–12; block 6' is a variant of block 6) to the present Leiden communities.

**Result.** Their anchor telomeres map systematically to the present Leiden communities:

| Ambrosini block | Anchor | Size | Copies | Community | Diagnostic gene(s) | Confirmed |
|-----------------|--------|------|--------|-----------|---------------------|-----------|
| 1 | chr1_p | 25 kb | 4 | C11 (chr1_p,chr5_q,chr6_q,chr8_p) | — | by arm membership |
| 2 | chr15_q | 88 kb | 7 | C8 (chr15_q singleton) | OR4F17 (1 arm, 2 samples), OR4F4 (1 arm, 409 samples) | gene-level |
| 3 | chr1_p | 35 kb | 1 | C11 | — | by arm membership |
| 5 | chr2_p | 17 kb | 5 | No PHR signal | RPL23AP7 (chr2_p excluded; no inter-chromosomal signal at ≥95%) | — |
| 6 | chr3_q | 38 kb | 5 | C3 | RPL23AP7-related | by arm membership |
| 6' | chr11_p | 11 kb | 1 | C3 | RYD5 (not in current annotations) | by arm membership |
| 7 (D4Z4) | chr4_q | 28 kb | 1* | C1 (chr4_q,chr10_q) | DUX4L1–DUX4L44 (28 genes, 2 arms) | gene-level |
| 8 (TUBB4q) | chr4_q | 14 kb | 6 | C1 | TUBB4q (not in current annotations) | by arm membership |
| 10 | chr2_q | 49 kb | 1 | C12 (chr2_q,chr20_p) | FBXO25 (not in current annotations) | by arm membership |
| 11 (IL9R) | chr9_q | 36 kb | 6 | C3 (chr3_q,chr7_p,chr9_q,chr11_p,chr16_q,chr19_p) | IL9RP1 (1 arm) | gene-level |
| 12 | chr12_p | 15 kb | 1 | C5 (chr6_p,chr9_p,chr12_p,chr20_q) | IQSEC3 (1 arm, 453 samples) | gene-level |

IL9R pseudogenes also appear in C2 (chr10_p/chr18_p: IL9RP2, IL9RP4), C6 (IL9RP4), C9 (chr16_p: IL9RP3, IL9R), and C14/C15 (PAR: IL9R, SPRY3), consistent with Mefford & Trask (2002) who reported IL9R pseudogenes at chr9_q, chr10_p, chr16_p, and chr18_p.

**Conclusion.** The concordance between single-genome duplicon blocks (Ambrosini 2007) and population-scale Leiden communities confirms that the community structure captures the known subtelomeric duplicon architecture.

### Ambrosini subterminal families → Leiden communities

**What it does.** Maps each of Ambrosini et al.'s (2007) 6 subterminal duplication families (Table 2, A–F) to the present Leiden communities.

**Result.** Their anchor telomeres map to the present Leiden communities:

| Family | Anchor | Size | Subterminal copies | Community | Key genes (Ambrosini) | Gene status |
|--------|--------|------|-------------------|-----------|----------------------|-------------|
| A | chr2_p | 7 kb | 6 | No PHR signal | RPL23AP7, FAM41C | chr2_p excluded (no inter-chromosomal signal at ≥95% identity); FAM41C detected in C11 (1 arm) |
| B | chr4_p | 17 kb | 10 | C13 (chr4_p) | RPL23AP7, FAM41C | C13 has DDX11L16 (1 arm); FAM41C not detected in C13 |
| C | chr9_p | 10 kb | 6 | C5 (chr6_p,chr9_p,chr12_p,chr20_q) | DDX11L-like, CXYorf1 | DDX11L family confirmed in C5 (7 members, 3–4 arms) |
| D | chr10_q | 22 kb | 10 | C1 (chr4_q,chr10_q) | RPL23AP7, FAM41C | RPL23A pseudogenes present in C1 (11 members, 1–2 arms); FAM41C not detected in C1 |
| E | chr17_p | 21 kb | 5 | C10 (chr17_p) | — | RPL23AP21/45/47/88 present (1 arm each) |
| F | chr18_p | 15 kb | 1 | C2 (chr10_p,chr18_p) | — | IL9RP2, IL9RP4 present (1 arm each) |

Family A (chr2_p) has no PHR signal — chr2_p is one of 6 arms excluded in the inter-chromosomal detection section for lacking inter-chromosomal signal at ≥95% identity. In Ambrosini et al.'s Table 2, Family A has 6 subterminal + 12 subtelomeric + 1 non-subtelomeric = 19 total copies, so it is predominantly subtelomeric. The lack of PHR signal at ≥95% identity likely reflects sequence divergence below the 95% threshold rather than absence of duplicated content.

FAM41C (lncRNA) is detected only in C11 (chr1_p/chr5_q/chr6_q/chr8_p, 1 arm, 215 samples), not at the anchor communities of families A, B, or D — indicating that the FAM41C copies associated with these subterminal families either lie below the 95% identity threshold or reside outside the PHR regions.

The RPL23A pseudogene family is the most widespread subtelomeric duplicon marker (RPL23AP45 spans 10 communities, 21 arms), but the specific member RPL23AP7 cited by Ambrosini is not annotated in the current Liftoff gene models; the RPL23A pseudogene nomenclature has been revised since Ambrosini's analysis.

**Family C's DDX11L-like gene signature is confirmed:** DDX11L family members (7 members) are present across 3–4 arms of C5, consistent with this family's inter-chromosomal distribution. Family F (chr18_p) maps to C2, which has the highest TAR1 density (mean 2.51 copies/sequence; the annotation section), consistent with chr18_p's known repeat-rich subtelomeric architecture.

**Conclusion.** Family A (chr2_p) has no PHR signal, consistent with sequence divergence below the 95% threshold. Family C's DDX11L-like gene signature is confirmed across 3–4 arms of C5.

### Community-specific genes

**What it does.** Identifies genes found in exactly one community, marking distinct subtelomeric identities.

| Community | N specific | Key genes | Biological significance |
|-----------|-----------|-----------|------------------------|
| C7 (acrocentric p-arms) | 48 | MTCO1P34, MTCO3P26, MTCO3P33, MTCO3P34, SNX18P15, SOWAHCP1, ASNSP2 | Mitochondrial pseudogenes (MTCO) and rDNA-associated loci; MTCO enrichment in acrocentric p-arm subtelomeres is a novel observation enabled by T2T-quality assemblies of these previously unresolved regions |
| C1 (chr4_q/chr10_q) | 26 | DUX4L pseudogenes (22), AGGF1P1, CLUHP4, DBET, LOC100996375 | D4Z4 macrosatellite array: DUX4L pseudogenes are copies of the DUX4 gene embedded within each D4Z4 repeat unit; the terminal copy on a permissive 4qA haplotype can produce pathogenic DUX4 protein causing FSHD (Lemmers et al. 2010); CTCF binds within D4Z4 as an insulator (Ottaviani et al. 2009). Chr10_q copies lack the stabilizing polyadenylation signal |
| C15 (PAR1) | 10 | SHOX, PPP2R3B, PLCXD1, GTPBP6, P2RY8, LOC124905300, LOC102724521, LINC00685, FABP5P13, KRT18P53 | 5 protein-coding genes including SHOX (short stature homeobox), a key growth regulator whose haploinsufficiency causes Leri-Weill dyschondrosteosis |
| C11 | 4 | FAM87B, LINC00115, LINC01409, LOC124903817 | lncRNAs |
| C3 (f7501) | 3 | FAM41AY1, FAM41AY2, LOC105375112 | lncRNAs specific to the 6-arm f7501 cluster |
| C14 (PAR2) | 1 | LOC124905309 | Single pseudogene |
| C6 | 1 | LOC124907874 | Single protein-coding gene |

### Hub genes

**What it does.** Identifies genes present in 3 or more communities — the common duplicon backbone spanning multiple community boundaries.

| Gene | N communities | N arms | Biotype |
|------|--------------|--------|---------|
| RPL23AP45 | 10 | 21 | pseudogene |
| SEPTIN14P22 | 9 | 22 | pseudogene |
| DDX11L16 | 9 | 20 | transcribed pseudogene |
| FAM138D | 9 | 17 | lncRNA |
| LOC101929828 | 9 | 21 | lncRNA |
| LOC102723681 | 9 | 23 | pseudogene |
| RPL23AP60 | 9 | 18 | pseudogene |
| RPL23AP87 | 9 | 19 | transcribed pseudogene |
| RPL23AP88 | 9 | 19 | pseudogene |

RPL23A pseudogenes (ribosomal protein L23a) and SEPTIN14 pseudogenes are the most widespread subtelomeric duplicon markers, consistent with their identification as core subtelomeric duplicon components by Ambrosini et al. (2007). The DDX11L (DEAD/H-box helicase), WASH (Wiskott-Aldrich syndrome protein homolog), and MIR6859 families also span 7–9 communities. These gene families are among the most widespread subtelomeric duplicon markers. The subtelomeric WASH copies are pseudogenes of WASHC1 (the catalytic subunit of the WASH complex).

**Conclusion.** The predominance of pseudogenes and ncRNA across all communities (28.6–86.4% pseudogene) is consistent with telomere position effect (TPE). Mefford & Trask (2002, citing Baur et al. 2001) noted that TPE operates in human cells, with "reporter genes near telomeres expressed at ten times lower levels," and suggested that subtelomeric regions might "buffer genes in chromosome-specific regions and in proximal subtelomeric domains from telomere-mediated repression." The exception — PAR1 (C15, 32.1% protein-coding) — involves a region with obligate crossover recombination. Fisher's exact tests for gene family enrichment per community (116 tests, BH-corrected) yield no significant results — the qualitative enrichments described above (MTCO in C7, DUX4L in C1, OR4F in C3) reflect presence patterns but do not survive multiple testing correction.

### Sequence-level vs arm-level enrichment comparison

**What it does.** Tests whether gene enrichment patterns are preserved when using the finer 50-community partition.

**Result.** The sequence-level enrichment uses 50 communities that are subpartitions of the 15 arm-level communities. Gene annotations are assigned at the arm level (the biotype section–the hub genes section), where each community groups 1–6 chromosome arms with shared gene content. At the 50-community level, most communities are pure single-arm (18/50) or near-pure (23/50), making gene enrichment largely redundant with the arm-level analysis. The arm-level enrichment (the biotype section–the hub genes section) provides the biologically meaningful gene characterization; the sequence-level partition's value is in resolving within-arm polymorphism (the polymorphic arm section, the heterogeneity section) rather than gene content.

**Conclusion.** The arm-level analysis captures the major duplicon block structure; the sequence-level analysis captures within-arm polymorphism.

### Singleton/doubleton QC

**What it does.** Identifies sequences assigned to a community where their arm is not the dominant arm (1–2 sequences only).

**Key metrics.** 29 singletons and 7 doubletons across 20 communities involving 25 arms.

**Conclusion.** These represent rare structural variants or borderline cluster assignments, not systematic errors.

### Files and scripts

**Scripts:**
| Script | Description |
|--------|-------------|
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze-community-enrichments.R` | Community gene enrichment testing |
| `/moosefs/guarracino/HPRCv2/scripts/community/plot-community-enrichments.R` | Enrichment result visualization |

**Output files:**
| File | Description |
|------|-------------|
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_summary_table.tsv` | Per-community summary: arms, genes, biotypes |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_enrichment.tsv` | Gene × community presence matrix |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_enrichment_fisher.tsv` | Fisher's exact test results |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_specific_genes.tsv` | Genes found in exactly one community |
| `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_biotypes.tsv` | Biotype distribution per community |


