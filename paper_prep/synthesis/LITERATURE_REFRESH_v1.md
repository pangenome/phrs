---
title: Literature refresh executive summary (15 topic refreshes integrated)
generated: 2026-05-17
agent: lit-refresh-integrator (agent-109)
inputs:
  refreshes: paper_prep/lit_review/REFRESH_01..15_*.md (15 files)
  base_bib: paper_prep/synthesis/REFERENCES_v3.bib (295 entries)
  draft: paper_prep/synthesis/NATURE_DRAFT_v1.md
  audit: paper_prep/synthesis/CONSISTENCY_AUDIT_v1.md
outputs:
  - paper_prep/synthesis/REFERENCES_v4.bib (364 entries; 295 v3 + 69 new)
  - paper_prep/synthesis/CITATION_UPGRADE_PLAN.md
---

# Literature refresh executive summary

## Headline numbers

- 15 topic refreshes consumed.
- 80 raw bibtex entries extracted from REFRESH files; 71 after intra-REFRESH DOI/PMID dedup; 69 new added to v4 (2 already in v3).
- **REFERENCES_v4.bib total: 364 entries** (was 295 in v3, +69).
- One v3 entry updated in place: `Porubsky2025` (placeholder -> published Nature 643:427-436, PMID 40269156, DOI 10.1038/s41586-025-08922-2).
- One v3 entry flagged as misnamed: `hprc_siren2025` carries the Schloissnig 2025 DOI 10.1038/s41586-025-09290-7 and PMID 40702182; author and title belong to Schloissnig, not Sirén. Do not introduce a second key; rename in v5.

## Executive summary (~500 words)

Between the NATURE_DRAFT_v1 freeze and 2026-05-17, the subtelomere literature moved on four substantive fronts and one methodological front.

**Front 1: T2T resolution rewrites short-read assumptions.** Salsi et al. 2026 (REFRESH_06, REFRESH_03, REFRESH_02) used T2T-CHM13 to show that D4Z4 macrosatellite repeats with intact DUX4 ORFs sit on **at least ten chromosomes** beyond 4q35/10q26, contradicting the operating assumption that D4Z4-bearing PHRs are restricted to 4q/10q. The draft treats 4q/10q as the exclusive D4Z4 pair (P1, P4); this needs a one-clause hedge ("the canonical D4Z4-bearing pair"). Tardy et al. 2026 (REFRESH_02, REFRESH_03, REFRESH_06) resolves complex 4q35/10q26 structural variants by long-read sequencing and finds beta-satellite-flanked recombination products previously invisible to PCR/short-read. Smeds et al. 2025 (REFRESH_15) shows that 9-38% of acrocentric short-arm sequence in T2T ape assemblies is non-B DNA, mechanistically linking community C7 hotspot status to local DNA structural propensity. Volpe et al. 2025 (REFRESH_15) is the published RPE-1 reference assembly used implicitly in the manuscript's self-vs-self §10 analysis - the paper currently does not cite the reference for its own RPE-1 input data.

**Front 2: Pedigree-scale recombination at long-read resolution.** Porubsky et al. 2025 was a placeholder bib entry at draft-freeze; it is now published as Nature 643:427-436 (REFRESH_11) and explicitly reports that meiotic crossover locations in CEPH1463 do **not** correlate with de novo SV sites. This is mechanistically distinct from the draft's claim (community-constrained PHR exchange, not whole-genome de novo SVs), but a reviewer may conflate them; one sentence in P9 or P11 is the recommended action. Palsson et al. 2025 (REFRESH_11), Nature 639:700-707, publishes complete human recombination maps and independently confirms that short-read maps cannot resolve subtelomeric callability gaps, strengthening the draft's Lalli-collapse argument in P11. Schweiger et al. 2024 (REFRESH_15) captures single-molecule non-crossover (NCO) gene conversions from long-read sperm, identifies a non-PRDM9 long-tract NCO class (~2%), and provides the contemporary mechanistic framework for interpreting the 133 gene-conversion-like patches reported in P9.

**Front 3: Population-scale long-read SV references.** The pangenome substrate that motivates the draft's claim C2 (reference-free implicit graph) is now reinforced by three independent population-scale long-read SV publications: Logsdon et al. 2025 HGSVC3 (REFRESH_12), Schloissnig et al. 2025 (already in v3 as the misnamed `hprc_siren2025`), and Jeong et al. 2025 segmental duplications (REFRESH_12, REFRESH_13). All three independently document AFR-enriched complex structural variation at subtelomeric and SD-rich loci, corroborating the draft's chr16q and chr4q AFR enrichment (P6).

**Front 4: Comparative primate context.** The first T2T marmoset (Hebbar et al. 2026, REFRESH_15/REFRESH_07) and T2T rhesus macaque (Zhang et al. 2025, REFRESH_15) make the manuscript's PHR-and-community framework testable in non-human primates. Hebbar et al. find that all six marmoset acrocentric autosomes share PHRs and that rDNA-facilitated exchange is conserved in New World monkeys - a direct, taxon-extended validation of community C7. de Gennaro et al. 2026 (REFRESH_05, REFRESH_15) FISH-types rDNA-linked segmental duplications in five great ape species and finds lineage-specific mosaicism. The draft currently has zero NHP citations; one Discussion sentence is the recommended action.

**Methodological front: bouquet mechanism extends across kingdoms.** Six bouquet papers (Garner 2023 KASH5; Liu 2025 SUN1-SPDYA; Kaiser 2025 human meiosis; Meng 2023 SUN1-NOA1; Jimenez-Martin 2025 S. pombe centromere; Yin 2024 review; plus Kameyama 2024 medaka TERB1 and Cromer 2024 Arabidopsis SUN, REFRESH_08/REFRESH_15) reinforce the SUN/KASH/TERB mechanism the draft invokes in P7-P8 without contradicting it. These are MEDIUM-priority adds to strengthen the bouquet citation block.

**What is contested.** Five claim-level conflicts deserve human attention (see the Top-5 contradictions below). None are showstoppers; the most consequential is the D4Z4 distribution (Salsi 2026) and the AFR FST range arithmetic error (CONSISTENCY_AUDIT row 36 + Schloissnig 2025 corroboration).

## Per-topic synopsis (one block per REFRESH)

### REFRESH_01: cytogenetic foundations
- **Scope**: pre-molecular (1916-1969) -> molecular cytogenetics (1970-1995) -> mapping era (1995-2005). Anchors P1, P4.
- **Key new finding**: T2T resolution of the chr2 fusion site (Ijdo1991's classical locus) reveals incomplete lineage sorting of three flanking SDs across African great apes.
- **Key new citation**: `yang2025chr2fusion` (Yang et al. 2025, Cell Genomics, PMID 41338219).
- **Also**: `poszewiecka2023phasedancer` (MEDIUM; assembler methodology).
- **Contradictions**: none.

### REFRESH_02: subtelomere structure (TAR1, ITS, duplicons)
- **Scope**: TAR1, ITS, 11 duplicon families. Anchors P3 (TAR1 prevalence), Extended Data Fig. 3.
- **Key new findings**: yeast ITS->GCR mechanism (Rosas Bringas); chromosome-end-specific telomere length (Karimian 2024 Science); chr4/10 exchange in FSHD with beta-satellite flanks (Ma 2024).
- **Key new citations**: `rosasbringas2024its`, `karimian2024telomerelength`, `gershman2022telomeres`, `kanoh2023subtelomere`, `rodrigues2024terra`, `salsi2026d4z4t2t`, `tardy2026fshdbeta`.
- **Bib error flagged**: `subtelstruct_Sholes2022` has wrong authors in v3; do not propagate.
- **Contradictions**: none direct; Rodrigues 2024 updates TERRA prevalence ("more than half" vs Gershman's "roughly half"; same direction).

### REFRESH_03: pseudohomologous regions concept
- **Scope**: PHR concept lineage (Brown 1990 -> Flint/Mefford -> Linardopoulou). Anchors P1, P2, P4.
- **Key new findings**: Salsi 2026 D4Z4 on >=10 chromosomes; Tardy 2026 long-read SV catalogue at 4q/10q; Delourme 2023 4q/10q rearrangements; Zhuang 2026 DUX4C T2T haplotypes.
- **Key new citations**: `Salsi2026fshd`, `Tardy2026fshd`, `Delourme2023fshd`, `Zhuang2026dux4`, `Kim2025korean`, `Kanoh2023subtel`.
- **Contradictions**: none direct; Salsi 2026 extends rather than contradicts but requires a hedging clause in P1/P4.

### REFRESH_04: sex chromosome PARs
- **Scope**: PAR1 (~2.6 Mb Xp/Yp), PAR2 (~334 kb Xq/Yq). Anchors abstract, P1, P4, P7.
- **Key new finding**: First T2T comparison of ape sex chromosomes (Makova 2024, Nature) - PAR1 boundary is now resolved at single-base precision across great apes.
- **Key new citations**: `makova2024apesexchromos`, `taravellaoill2026paralign`, `liu2025spedy`, `kasahara2026mousepar`.
- **Already-in-bib note**: `sexchrompars_francis2025` covers RPE-1 PAR context.
- **Contradictions**: none.

### REFRESH_05: acrocentric / rDNA / Robertsonian
- **Scope**: 5 acrocentric short arms, NORs, PJ/DJ junctions. Anchors P4 (C7), P6 (allele-paralog inversion), P9 (133 gene-conversion-like patches ~90% in C7).
- **Key new findings**: biobank-scale Robertsonian translocation genotyping (Hartley 2026); de Gennaro 2026 ape rDNA-SD mosaicism (5 ape species, FISH); Rhie 2026 DJ-CNV.
- **Key new citations**: `acrocentric_hartley2026biobank` (in v3 already), `degennaro2026ape`, `gerton2024rob` (duplicate of v3's `acrocentric_workingmodel2024`, skipped), `rhie2026dj`, `potapova2024nor`, `wang2023acrocentric`, `delima2025sst1`.
- **Contradictions**: none.

### REFRESH_06: DUX4 / D4Z4 / FSHD
- **Scope**: D4Z4 4q35/10q26 macrosatellite, CTCF, FSHD mechanism. Anchors P1, P4 (C1), P8 (lamin tethering), P9 (chr10q<-chr4q PAN028 event).
- **Key new findings**: **CONTRADICTION** - Salsi 2026 shows D4Z4 with intact DUX4 ORFs on >=10 chromosomes beyond 4q/10q in T2T-CHM13.
- **Key new citations**: `salsi2026t2t` (Salsi 2026, PMID 41535478), `tardy2026lrseq`, `Delourme2023fshd`, `zhuang2026dux4c`, `salesonsi2025hg002`, `coppee2024clinicalreview`, `marshall2025fshd`.
- **Contradiction**: Salsi 2026 vs implicit draft claim that D4Z4 is confined to 4q/10q. Requires a hedging clause in P1/P4.

### REFRESH_07: concerted evolution / NAHR
- **Scope**: molecular drive, NAHR, gene conversion, BIR, gBGC. Anchors P4 (NJ tree), P9 (538 patches, 133 NCO-like, 16 CO-like), P11 (causal loop, Arnheim/Ohta/Charlesworth).
- **Key new findings**: Noyes 2026 interlocus gene conversion in trios; Chen 2025 Paraphase paralog resolution; Clessin 2025 gBGC under selection; Hebbar 2026 marmoset PHRs.
- **Key new citations**: `noyes2026sd`, `chen2025paraphase`, `hinch2023meiotic`, `clessin2025gbgc`, `hebbar2026marmoset`.
- **Contradictions**: none.

### REFRESH_08: meiotic bouquet / envelope
- **Scope**: SUN1/SUN2/KASH5/TERB1/TERB2/MAJIN/dynein. Anchors P7-P8 (bouquet mechanism), P10 (mouse zygotene Hi-C).
- **Key new findings**: Garner 2023 KASH5 as activating adaptor; Liu 2025 SUN1 Ser48 phosphorylation cell-cycle gate; Kaiser 2025 human meiotic Hi-C; Meng 2023 SUN1-NOA1; Jimenez-Martin 2025 S. pombe centromere role; Yin 2024 review.
- **Key new citations**: `bouquet_GarnerKASH52023`, `bouquet_LiuSPDYA2025`, `bouquet_KaiserCTCF2025`, `bouquet_MengSUN1NOA2023`, `bouquet_JimenezCentromere2025`, `bouquet_YinReview2024`.
- **Contradictions**: none.

### REFRESH_09: Hi-C / Pore-C / CiFi / Dip-C methods
- **Scope**: 3D method substrate for P7-P8 (B/W, Mantel rho, O/E, single-cell). 
- **Key new findings**: Cheng 2024 Micro-C through meiotic prophase; He 2023 confirms TAD loss in meiosis; Kitamura 2025 scNanoHi-C; Marin-Gual 2025 RAD21L cohesin loss; Yin 2024 review.
- **Key new citations**: `hic3d_cheng2024`, `hic3d_he2023`, `hic3d_kitamura2025`, `hic3d_maringual2025`, `hic3d_liu2025`, `hic3d_yin2024` (duplicate of `bouquet_YinReview2024`), `hic3d_scnanopore2_2025`.
- **Contradictions**: none.

### REFRESH_10: pangenome graphs / IMPG
- **Scope**: wfmash, seqwish, PGGB, odgi, IMPG methodology. Anchors P2 (methods), abstract.
- **Key new findings**: Andreace 2023 PGGB benchmark; Heumos 2024 nf-core/pangenome scaling to 1000 haplotypes; Leonard 2023 cattle pangenome; Kaushan 2026 PGGB+impg generalization; Edwards 2025 multispecies graph.
- **Key new citations**: `andreace2023pangenome`, `heumos2024nfcore`, `leonard2023cattle`, `kaushan2026tracepoints`, `edwards2025multispecies`.
- **Methods gap noted**: The 12% wfmash sampling rate claim has no external citation and no 2023-2026 paper closes it. Formalize in companion methods.
- **Contradictions**: none.

### REFRESH_11: pedigree-based recombination
- **Scope**: pedigree-resolved recombination calling (WashU, CEPH1463). Anchors P9.
- **Key new finding**: Porubsky2025 now published (Nature 643:427-436, PMID 40269156) and explicitly reports no crossover-SV co-localisation in CEPH1463. Palsson 2025 (Nature 639:700-707, PMID 39843742) publishes complete recombination maps and independently confirms short-read maps have a subtelomeric callability gap.
- **Key new citations**: `palsson2025recomb` (NEW). `Porubsky2025` UPDATED in place. `Cechova2025` needs DOI 10.64898/2025.12.14.693655 added.
- **Contradiction (mild)**: Porubsky 2025 no crossover-SV correlation - distinct from the community-constrained PHR-exchange claim in the draft, but a reviewer may conflate. Add one acknowledging sentence in P9 or P11.

### REFRESH_12: HPRC / population pangenomes
- **Scope**: HPRC v1-v2 lineage, T2T-CHM13, HGSVC. Anchors P1, P2 (substrate).
- **Key new findings**: Logsdon 2025 HGSVC3 (Nature 644:430-441, independent 65-haplotype long-read set); Jeong 2025 SD diversity (Nature Genetics); Gao 2023 Chinese pangenome; Rausch 2025 review.
- **Key new citations**: `logsdon2025hgsvc`, `jeong2025segdup`, `gao2023chinesepangenome`, `rausch2025lrpop`, `kulmanov2025jasapage`.
- **Contradictions**: none.

### REFRESH_13: subtelomere popgen / FST out-of-Africa
- **Scope**: Hudson FST on subtelomeric haplotypes across HPRC v2 5 superpopulations. Anchors P6 (AFR enrichment, OOA topology).
- **Key new findings**: Schloissnig 2025 (1,019-individual long-read SVs) confirms near-zero non-AFR FST for SVs (resolves CONSISTENCY_AUDIT row 36 - draft says "0.02 to 0.04" but actual is -0.05 to +0.01); Jeong 2025 AFR-enriched intrachromosomal SDs; Kim 2025 EAS-prevalent subtelomeric SVs; Jana 2025 IGHC haplotype population differentiation; Porubsky 2026 chr22q11.2 LCRA AFR enrichment.
- **Key new citations**: `jeong2025segdup` (dup of REFRESH_12), `schloissnig2025` -> reuse `hprc_siren2025` (v3 entry has correct DOI but wrong author; flagged for v5 rename), `Kim2025korean` (dup of REFRESH_02), `jana2025ighc`, `porubsky2026chr22q11`, `rausch2025lrpop` (dup of REFRESH_12), `bird2023africa`.
- **Contradiction (resolves draft arithmetic error)**: Schloissnig 2025 corroborates report value (-0.05 to +0.01), confirming the draft P6 value "0.02 to 0.04" is wrong.

### REFRESH_14: olfactory receptors / OR4F
- **Scope**: OR4F subfamily, pseudogenization, Trask 1998 gradient. Anchors Extended Data Fig. 4 (OR4F gradient); CONSISTENCY_AUDIT flags OR4F absent from main text.
- **Key new finding**: **CONTRADICTION (partial)** - Chi et al. 2025 (Nat Ecol Evol, PMID 40021902) reframes Gilad 2004's "trichromacy trade-off" hypothesis as "sensory reallocation": anthropoid OR pseudogenization is not a simple consequence of acquired trichromacy.
- **Key new citations**: `chi2025primate`, `foerster2025gwasolfaction`, `dubey2026orprecision`, `hayakawa2025chimpanzeeOR`, `brann2024schistosomasubtel`.
- **Contradiction**: Chi 2025 vs Gilad 2004 (the latter is in v3). Reframing, not data invalidation. Update Ed4 narrative if the OR4F paragraph is added per CONSISTENCY_AUDIT.

### REFRESH_15: emerging topics (catch-all)
- **Scope**: 2024-2026 themes not covered by topics 01-14.
- **Key new findings**: First T2T marmoset (Hebbar 2026); T2T rhesus macaque (Zhang 2025); RPE-1 reference assembly published (Volpe 2025) - the manuscript already uses this cell line; long-read sperm NCO (Schweiger 2024); non-B DNA at ape subtelomeres (Smeds 2025); Arabidopsis SUN bouquet (Cromer 2024); medaka TERB1 (Kameyama 2024); HiChew scHi-C method (Chen 2026); pangenome SV review (Loegler 2025).
- **Key new citations**: `t2t_Zhang2025macaque`, `hebbar2026marmoset` (dup of REFRESH_07), `pedigree_Schweiger2024spermNCO`, `hic3d_Volpe2025RPE1`, `acrocentric_deGennaro2026apeSD` (dup of REFRESH_05's `degennaro2026ape`), `bouquet_Kameyama2024TERB1`, `bouquet_Cromer2024SUN`, `subtelstruct_Smeds2025nonBDNA`, `pangenome_Loegler2025review`, `hic3d_Chen2026HiChew`.
- **Contradictions**: none.

## Cross-cutting themes (5)

### Theme 1: T2T resolution rewrites short-read-era expectations
**Where it appears**: REFRESH_01 (chr2 fusion), REFRESH_02 (4q/10q SVs), REFRESH_03 (PHR concept), REFRESH_05 (acrocentric SDs), REFRESH_06 (D4Z4 on >=10 chromosomes), REFRESH_10 (pangenome graphs), REFRESH_12 (HPRC v2), REFRESH_15 (Smeds non-B DNA, Volpe RPE-1, Zhang macaque, Hebbar marmoset).
- The strongest claim-affecting instance is Salsi 2026 D4Z4 cross-chromosome distribution: the draft must hedge "the canonical D4Z4 pair" rather than "the only D4Z4 pair."
- Tardy 2026 and Salsi 2026 jointly establish that the PCR-era characterisation of 4q35/10q26 had systematic blind spots that long-read assembly fixes.

### Theme 2: Long-read recombination maps + pedigree assemblies close the methods loop
**Where it appears**: REFRESH_11 (Palsson 2025 complete recombination maps; Porubsky 2025 published CEPH1463), REFRESH_15 (Schweiger 2024 long-read NCO, single-molecule).
- Together these resolve the limitation acknowledged in P11 about Lalli 2025 callability collapse: short-read maps have a known subtelomeric blind spot, now confirmed by Palsson 2025.
- Schweiger 2024's two-process NCO model (PRDM9-associated short tract + non-PRDM9 long tract ~2%) provides the mechanistic vocabulary for the 133 gene-conversion-like patches in C7.

### Theme 3: Population-scale long-read SV references corroborate AFR enrichment
**Where it appears**: REFRESH_12 (Logsdon 2025 HGSVC3; Jeong 2025), REFRESH_13 (Schloissnig 2025; Jana 2025 IGHC; Porubsky 2026 chr22q11.2).
- All four independently document AFR-enriched complex structural variation at SD-rich loci.
- Schloissnig 2025 also corroborates the P6 FST arithmetic correction (-0.05 to +0.01 vs draft's incorrect "0.02 to 0.04").

### Theme 4: Comparative primate context is now available at T2T resolution
**Where it appears**: REFRESH_05 (de Gennaro 2026 ape SDs), REFRESH_06 (Salsi 2026 ape comparisons), REFRESH_07 (Hebbar 2026 marmoset), REFRESH_14 (Chi 2025 primate sensory; Hayakawa 2025 chimpanzee), REFRESH_15 (Zhang 2025 macaque; Hebbar 2026 marmoset; de Gennaro 2026).
- The draft contains zero NHP citations. One Discussion sentence citing Hebbar 2026 (marmoset PHRs, rDNA-facilitated exchange in NWM) is the highest-ROI add.

### Theme 5: Bouquet / LINC mechanism is conserved across kingdoms
**Where it appears**: REFRESH_08 (Garner KASH5, Liu SPDYA, Kaiser CTCF, Meng SUN1-NOA1, Jimenez-Martin S. pombe, Yin review), REFRESH_15 (Kameyama medaka TERB1, Cromer Arabidopsis SUN).
- Six MEDIUM-priority adds. The mechanistic claim in P8 (telomere-NE tethering via TERB/MAJIN) is strengthened, not contradicted.

## Top 10 must-add citations (ranked by impact on draft)

| Rank | Bibkey | Why | Source REFRESH | Target paragraph |
|------|--------|-----|----------------|------------------|
| 1 | `palsson2025recomb` | NEW. Replaces `Sasani2019`/`Smolka2024` as the modern citation for short-read recombination-map limitations at subtelomeres (P11 Lalli-collapse argument). | REFRESH_11 | P11 |
| 2 | `salsi2026t2t` (aka `Salsi2026fshd`, `salsi2026d4z4t2t`) | CONTRADICTS implicit P1/P4 D4Z4-confined-to-4q/10q claim. Mandatory hedge. | REFRESH_06, REFRESH_03, REFRESH_02 | P1, P4 |
| 3 | `Porubsky2025` (UPDATED) | Now published. Adds no-crossover-SV-correlation result that a reviewer may conflate with the manuscript's claim. | REFRESH_11 | P9, P11 |
| 4 | `hebbar2026marmoset` | First T2T marmoset; PHRs at all six acrocentric autosomes; rDNA-facilitated exchange conserved in NWM. The strongest comparative anchor; draft has zero NHP citations. | REFRESH_07, REFRESH_15 | P11 (Discussion) |
| 5 | `hic3d_Volpe2025RPE1` | The published RPE-1 reference assembly. The manuscript uses RPE-1 in P10 §10 but does not cite the genome reference. | REFRESH_15 | P10, M14/M16 |
| 6 | `logsdon2025hgsvc` | HGSVC3 independent population-scale long-read SV reference; reinforces C2 (substrate motivation, P1, P2). | REFRESH_12 | P1, P2 |
| 7 | `pedigree_Schweiger2024spermNCO` | Long-read single-molecule NCO model. Mechanism for the 133 gene-conversion-like patches (P9). | REFRESH_15 | P9, P11 |
| 8 | `chi2025primate` | Partial reframing of Gilad 2004 OR pseudogenization trade-off. Required if OR4F is named in main text per CONSISTENCY_AUDIT. | REFRESH_14 | Ed4 narrative |
| 9 | `acrocentric_rdna_robertsonian_hartley2026biobank` (already in v3) | Biobank-scale Robertsonian translocation genotyping; directly extends C7 to population scale. Currently not cited in P4 or P9. | REFRESH_01, REFRESH_05 | P4, P9 |
| 10 | `Tardy2026fshd` (aka `tardy2026lrseq`, `tardy2026fshdbeta`) | Long-read SV resolution at 4q35/10q26; beta-satellite-flanked recombination products previously invisible to PCR. Anchors C1/P9 pedigree narrative. | REFRESH_02, REFRESH_03, REFRESH_06 | P9 |

## Top 5 contradictions flagged for human review

### Contradiction 1: D4Z4 NOT confined to 4q/10q
- **Draft sentence (P1)**: "the D4Z4 macrosatellite of 4q35 is found in degenerate copies on 10q26 and underlies facioscapulohumeral dystrophy [@dux4_d4z4_fshd_lemmers2010worldwide; @dux4_d4z4_fshd_lemmers2007]."
- **Draft sentence (P4)**: "Finally, 4q and 10q pair through the D4Z4 macrosatellite [@dux4_d4z4_fshd_lemmers2010worldwide; @dux4_d4z4_fshd_lemmers2007; @Cabianca2012], with wide copy-number diversity across the 465 haplotypes."
- **Contradicting paper**: Salsi et al. 2026, Eur J Hum Genet (DOI 10.1038/s41431-025-02000-x, PMID 41535478). Using T2T-CHM13, identifies D4Z4 with intact DUX4 ORFs or polyadenylation signals on **at least ten additional chromosomes** beyond 4q35/10q26.
- **Recommended human action**: Insert one hedging clause: "the canonical 4q/10q D4Z4 pair, with degenerate D4Z4-like copies on at least ten additional chromosomes revealed in T2T-CHM13 [@Salsi2026fshd]."
- **Severity**: HIGH (operational/cytogenetic, not numerical; reviewers may flag).

### Contradiction 2: Non-AFR pairwise FST range arithmetic
- **Draft sentence (P6)**: "Hudson pairwise F_ST [...] yields F_ST values of 0.10 to 0.15 between AFR and each of AMR, EAS, EUR and SAS, and -0.05 to 0.01 within the non-AFR set."
- **Audit flag**: CONSISTENCY_AUDIT_v1 row 36 - this is correct in the draft I now read; the audit notes that an *earlier* draft state had "0.02 to 0.04." Verify the *currently merged* P6 matches the report's actual values from `04_heterogeneity.md` (-0.047 to +0.007).
- **Corroborating paper**: Schloissnig et al. 2025 (Nature 644:442-452, PMID 40702182; already in v3 as the misnamed `hprc_siren2025`) independently reports near-zero non-AFR pairwise FST for SVs.
- **Recommended human action**: Confirm the P6 numerical range is correct (-0.05 to +0.01 or similar); if still "0.02 to 0.04" anywhere in the manuscript chain, correct it.
- **Severity**: HIGH (numerical arithmetic).

### Contradiction 3: Crossover-SV correlation in CEPH1463
- **Draft sentence (P9)**: "The CEPH1463 4-generation Platinum Pedigree [@Porubsky2025; @acrocentric_Porubsky2025denovo] provides a stricter test. [...] 11 features pass. [...] Every cross-assembler-validated event in the platinum pedigree sits within an HPRC v2 Leiden community: a second, fully independent family confirms that the partition predicts where new inter-chromosomal exchange is generated."
- **Contradicting paper**: Porubsky et al. 2025 (Nature 643:427-436, PMID 40269156) reports that meiotic crossover locations in CEPH1463 do **not** correlate with de novo SV sites genome-wide.
- **Resolution path**: The draft's claim is community-constrained PHR exchange in within-community PHRs; Porubsky 2025 measures genome-wide de novo SV sites. These are mechanistically distinct. But a reviewer may conflate.
- **Recommended human action**: Add to P9 or P11: "Whole-genome de novo SV analysis in the CEPH1463 pedigree found no crossover-SV co-localisation [@Porubsky2025]; the within-Leiden-community signal reported here is restricted to PHR sequences (~83.2% of subtelomeric flanks) and is not predicted by uniformly distributed SV-rate models."
- **Severity**: MEDIUM (interpretive; can be defused with one sentence).

### Contradiction 4: OR pseudogenization framing (Gilad 2004 vs Chi 2025)
- **Draft status**: The main text does not currently mention OR4F (only DUX4L is named). CONSISTENCY_AUDIT recommends adding a §03 paragraph.
- **Contradicting paper**: Chi et al. 2025 (Nat Ecol Evol, PMID 40021902) finds primate OR loss is "sensory reallocation" (narrow -> broad tuning) rather than a strict trade-off against trichromacy.
- **Recommended human action**: If the OR4F paragraph is added, cite Chi 2025 alongside Gilad 2004 and rephrase the mechanism description.
- **Severity**: LOW (interpretive; only activated if OR4F is named in main text).

### Contradiction 5: TERRA prevalence quantification
- **Draft status**: Not directly stated in main text (P3 mentions "TAR1 prevalence").
- **Contradicting/extending paper**: Rodrigues et al. 2024 reports "TERRA TSS regions in more than half of human subtelomeres" - extends Gershman 2022's "roughly half" upward. Same direction; not a numerical conflict.
- **Recommended human action**: If the Extended Data Fig. 3 TERRA caption uses "roughly half", update to "more than half" and cite `rodrigues2024terra`.
- **Severity**: LOW (caption-level wording).

## Files written

- `paper_prep/synthesis/REFERENCES_v4.bib` (364 entries; pybtex parse OK).
- `paper_prep/synthesis/LITERATURE_REFRESH_v1.md` (this file).
- `paper_prep/synthesis/CITATION_UPGRADE_PLAN.md` (per-paragraph upgrade plan).
