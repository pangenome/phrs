---
title: Per-paragraph citation upgrade plan
draft: paper_prep/synthesis/NATURE_DRAFT_v1.md (lines 18-84)
bib_source: paper_prep/synthesis/REFERENCES_v4.bib (364 entries)
generated: 2026-05-17
agent: lit-refresh-integrator (agent-109)
---

# Per-paragraph citation upgrade plan

All bibkeys referenced below are present in `REFERENCES_v4.bib`. The action verbs are:
- **ADD**: insert a new citation; the existing block remains.
- **REPLACE x -> y**: swap an existing citation for a newer one (the old entry remains in v4 but no longer appears in the draft).
- **REMOVE**: drop without replacement.

## Summary table

| Paragraph | ADDs | REPLACEs | REMOVEs | Net effect |
|-----------|------|----------|---------|------------|
| Abstract | 0 | 0 | 0 | unchanged |
| P1 (history) | 4 | 0 | 0 | +4 strong citations, no removals |
| P2 (methods substrate) | 2 | 0 | 0 | +2 new pangenome substrate citations |
| P3 (heatmaps, scale) | 1 | 0 | 0 | +1 (Volpe RPE-1 for inter-method context) |
| P4 (NJ tree) | 3 | 0 | 0 | +3 (Salsi hedge, marmoset, ape rDNA-SDs) |
| P5 (3-architecture split) | 0 | 0 | 0 | unchanged |
| P6 (heterogeneity, FST) | 4 | 0 | 0 | +4 popgen / SV citations; verify FST range arithmetic |
| P7 (3D Hi-C) | 3 | 0 | 0 | +3 (Cheng, He, Kaiser, Marin-Gual) |
| P8 (flanking paradox, bouquet) | 5 | 0 | 0 | +5 bouquet conservation citations |
| P9 (pedigree exchange) | 4 | 0 | 0 | +4 (Schweiger NCO, Salsi/Tardy D4Z4, Porubsky-published note) |
| P10 (RPE-1, mouse) | 2 | 0 | 0 | +2 (Volpe RPE-1, Zhang macaque optional) |
| P11 (causal loop, limits) | 3 | (1) | 0 | +2 net; supplement Sasani/Smolka with Palsson |
| Methods M1..M17 | 6 | 0 | 0 | +6 (pangenome, scNanoHiC method, sperm NCO) |
| **TOTAL** | **37** | **0** | **0** | **+37 new citations; 1 supplement (not replace)** |

(There are zero true REPLACEs in this plan: every "supersedes" relationship in the REFRESH files was assessed as additive rather than replacing. REMOVEs are also zero - no v3 citation was flagged as actively wrong, only stale Porubsky2025 which is updated in place.)

---

## P1: Main text intro (history, founding references) - line 22

**Currently cited (23 bibkeys)**: `@Brown1990`, `@Wilkie1991`, `@Trask1991`, `@Trask1998`, `@Rouyer1986`, `@sexchrompars_acquaviva2020`, `@sexchrompars_bellott2024`, `@acrocentric_Altemose2022`, `@acrocentric_Guarracino2025ape`, `@dux4_d4z4_fshd_lemmers2010worldwide`, `@dux4_d4z4_fshd_lemmers2007`, `@Mefford2001`, `@MeffordTrask2002`, `@Linardopoulou2005`, `@Riethman2001`, `@Riethman2004`, `@Riethman2008`, `@Flint1997`, `@Ambrosini2007`, `@Nurk2022`, `@Logsdon2021`, `@hprc_hprcv2_2025`, `@Liao2023`.

**Actions:**
- **ADD `@Salsi2026fshd`** (REFRESH_06, REFRESH_03). On the D4Z4 clause "is found in degenerate copies on 10q26" - hedge with "the canonical 4q/10q D4Z4 pair, with degenerate D4Z4-like loci detected on at least ten additional chromosomes in T2T-CHM13 [@Salsi2026fshd]." HIGH PRIORITY.
- **ADD `@yang2025chr2fusion`** (REFRESH_01). If P1 cites Ijdo 1991 or refers to the chr2 fusion as an example of "ancient PHR preserved in interstitial position" (this is in report §12 - check whether the manuscript chain includes it), add Yang 2025 as the modern T2T extension.
- **ADD `@logsdon2025hgsvc`** (REFRESH_12). After the HPRC v2 citation block (`hprc_hprcv2_2025; Liao2023`), supplement with `@logsdon2025hgsvc` (Nature 644:430-441) as independent population-scale long-read SV reference.
- **ADD `@acrocentric_rdna_robertsonian_hartley2026biobank`** (REFRESH_01, REFRESH_05) on the acrocentric clause (currently citing Altemose 2022 + Guarracino 2025): supplement with biobank-scale Robertsonian translocation evidence.

**Net effect: +4 strong citations; no removals.**

---

## P2: Main text methods substrate, PHR definition - line 24

**Currently cited (7 bibkeys)**: `@Guarracino2023`, `@pangenome_graphs_impg_GarrisonGuarracino2023`, `@pangenome_graphs_impg_GuarracinoHeumos2022`, `@pangenome_graphs_impg_IMPG2023`, `@pangenome_graphs_impg_Hickey2024`, `@Garrison2018`, `@Garrison2024pggb`.

**Actions:**
- **ADD `@andreace2023pangenome`** (REFRESH_10). On the PGGB clause: independent benchmark confirming PGGB's superiority for complete variation capture.
- **ADD `@heumos2024nfcore`** (REFRESH_10). On the pipeline scaling clause (or methods chain): nf-core/pangenome scaling to 1000 haplotypes is contemporaneous with our 466-haplotype operation. Lower priority than #1.

**Methods gap**: The 12% wfmash sampling rate and Erdős-Rényi connectivity argument remain *uncited externally*. No new paper closes this. The manuscript or its companion methods paper must formalize this internally.

**Net effect: +2 new pangenome substrate citations.**

---

## P3: Identity heatmaps, PHR scale, comparison to PAR2 - line 26

**Currently cited (7 bibkeys)**: `@sexchrompars_bellott2024`, `@Bailey2002`, `@BaileyEichler2006`, `@RuizHerrera2008`, `@Vollger2023`, `@concerted_evolution_nahr_Vollger2023`, `@Stong2014`.

**Actions:**
- **ADD `@hic3d_Volpe2025RPE1`** (REFRESH_15). LOW PRIORITY here; primary location is P10. Skip in P3.

(No meaningful upgrades for P3: the PHR-scale comparison is anchored on PAR2 length and the existing concerted-evolution citation block remains current.)

**Net effect: +0 net change is acceptable. Keep P3 unchanged.**

---

## P4: Main text - NJ tree, six clades - line 28

**Currently cited (10 bibkeys)**: `@sexchrompars_acquaviva2020`, `@acrocentric_Altemose2022`, `@acrocentric_Guarracino2025ape`, `@Linardopoulou2005`, `@dux4_d4z4_fshd_lemmers2010worldwide`, `@dux4_d4z4_fshd_lemmers2007`, `@Cabianca2012`, `@Skaletsky2003`, `@Rudd2009`, `@Mefford2001`.

**Actions:**
- **ADD `@Salsi2026fshd`** on the 4q/10q clade clause. Same hedge as P1 if not centralised there.
- **ADD `@hebbar2026marmoset`** (REFRESH_07, REFRESH_15) on the acrocentric clade clause: NWM PHR + rDNA-facilitated exchange in all 6 marmoset acrocentric autosomes is a direct cross-species validation of community C7.
- **ADD `@degennaro2026ape`** (REFRESH_05, REFRESH_15) on the acrocentric clade clause: lineage-specific rDNA-linked SD mosaicism in 5 great ape species.

**Net effect: +3 strong citations (T2T hedge + comparative anchors).**

---

## P5: Three-architecture partition (4/28/9) - line 30

**Currently cited (3 bibkeys)**: `@Linardopoulou2005`, `@Stong2014`, `@Eichler2001`.

**Actions:** none. The 4/28/9 split is a manuscript-original analytic claim; the existing citation block is the correct prior-FISH-era benchmark.

**Net effect: 0.**

---

## P6: Within-community heterogeneity, F_ST, out-of-Africa - line 32

**Currently cited (18 bibkeys)**: `@Flint1997`, `@MeffordTrask2002`, `@acrocentric_Altemose2022`, `@acrocentric_rdna_robertsonian_bandyopadhyay2001`, `@Mefford2001`, `@subtelstruct_NergadzeITS2007`, `@subtelstruct_Nergadze2007`, `@subtelstruct_NergadzeITSReview2007`, `@Ambrosini2007`, `@subtel_popgen_hudson1992`, `@subtel_popgen_weir1984`, `@subtel_popgen_lewontin1972`, `@Bergstrom2020`, `@subtel_popgen_anderson2008`, `@subtel_popgen_bhatia2013`, `@subtel_popgen_rosenberg2002`, `@subtel_popgen_levysakin2019`, `@subtel_popgen_1000g2010`.

**Actions:**
- **ADD `@jeong2025segdup`** (REFRESH_12, REFRESH_13) on the AFR enrichment clause: genome-wide AFR enrichment of SDs in 85-sample HPRC-derived assemblies. HIGH PRIORITY (corroborates AFR direction).
- **ADD `@porubsky2026chr22q11`** (REFRESH_13) on the AFR enrichment clause: AFR LCRA at chr22q11.2 is significantly longer and architecturally more complex.
- **ADD `@hprc_siren2025`** (already in v3; bibkey is misnamed - it carries Schloissnig 2025 PMID 40702182). Cite alongside Bergstrom 2020 / 1000G as the long-read population-scale FST benchmark for SVs.
- **ADD `@bird2023africa`** (REFRESH_13) as MEDIUM priority for the within-Africa structure caveat (audit flagged in REFRESH_13 §3 MODERATE-2).

**Numerical check required**: Verify that P6's "F_ST values of [...] -0.05 to 0.01 within the non-AFR set" matches the actual report value (-0.047 to +0.007). CONSISTENCY_AUDIT row 36 flagged an *earlier* draft with "0.02 to 0.04"; this appears already corrected in NATURE_DRAFT_v1 line 32, but verify.

**Bib housekeeping note**: `hprc_siren2025` should be renamed (e.g. to `hprc_schloissnig2025` or `subtel_popgen_schloissnig2025`) in REFERENCES_v5; for now, the DOI is the durable identifier.

**Net effect: +4 new popgen / SV citations; numerical verification on FST range.**

---

## P7: Main text - 3D Hi-C, bouquet, multi-method - line 34

**Currently cited (11 bibkeys)**: `@Tan2018`, `@Xu2025`, `@Zuo2021`, `@Cechova2025`, `@Ulahannan2019`, `@hic3d_cifi2025`, `@hic3d_dixon2012`, `@hic3d_imakaev2012`, `@hic3d_alavattam2019`, `@hic3d_wolff2018`, `@hic3d_deshpande2022`.

**Actions:**
- **ADD `@hic3d_cheng2024`** (REFRESH_09) on the meiotic TAD/loop clause: Micro-C through meiotic prophase confirms CTCF-anchored loop bases on axes during meiosis. Supports the "PHR fits within one meiotic loop" implicit claim.
- **ADD `@bouquet_KaiserCTCF2025`** (REFRESH_08, REFRESH_09) for human meiotic Hi-C - the closest available human meiotic reference.
- **ADD `@hic3d_kitamura2025`** (REFRESH_09) on the scnanoHiC single-cell clause: scNanoHi-C is the contemporary method context for `hic3d_scnanoHiC2023`/`hic3d_scnanoHiC2_2025`.

**Net effect: +3 strong citations.**

---

## P8: Flanking paradox + S_all negative control + bouquet mechanism - line 36

**Currently cited (17 bibkeys)**: `@Xu2025`, `@bouquet_KotaSUN1MAJIN2020`, `@bouquet_Scherthan2001`, `@bouquet_Scherthan2003`, `@bouquet_ShibuyaRPMs2015`, `@bouquet_ChikashigeTelomere1994`, `@bouquet_HarperBouquet2004`, `@bouquet_HornKASH52013`, `@bouquet_DingSUN12007`, `@bouquet_MorimotoKASH2012`, `@bouquet_ZicklerKleckner1999`, `@ZicklerKleckner1998`, `@ZicklerKleckner2015`, `@Ottaviani2009`, `@OttavianiGilson2008`, `@Masny2004`, `@Cabianca2012`.

**Actions (all on the bouquet citation block):**
- **ADD `@bouquet_GarnerKASH52023`** (REFRESH_08). KASH5 as activating dynein adaptor; refines the LINC complex mechanism.
- **ADD `@bouquet_LiuSPDYA2025`** (REFRESH_08). SUN1 Ser48 phosphorylation cell-cycle gate.
- **ADD `@bouquet_MengSUN1NOA2023`** (REFRESH_08). SUN1-NOA1 axis.
- **ADD `@bouquet_JimenezCentromere2025`** (REFRESH_08). S. pombe centromere role in bouquet - cross-kingdom mechanistic generality.
- **ADD `@subtelstruct_Smeds2025nonBDNA`** (REFRESH_15) on the D4Z4-CTCF-lamin clause: non-B DNA structures (G4, Z-DNA, hairpins) occupy 9-38% of acrocentric short arms and subtelomeres in ape T2T genomes - mechanistic complement to the community-defines-exchange-hotspot framing.

**Net effect: +5 bouquet/mechanism citations.**

---

## P9: Pedigree-resolved exchanges - line 38

**Currently cited (10 bibkeys)**: `@Cechova2025`, `@concerted_evolution_nahr_SamonteEichler2002`, `@concerted_evolution_nahr_Eichler2001`, `@concerted_evolution_nahr_Hastings2009`, `@concerted_evolution_nahr_Myers2010`, `@Sharp2006`, `@StankiewiczLupski2002`, `@StankiewiczLupski2010`, `@Porubsky2025`, `@acrocentric_Porubsky2025denovo`.

**Actions:**
- **ADD `@pedigree_Schweiger2024spermNCO`** on the 133 gene-conversion-like clause: long-read single-molecule NCO from sperm; two-process NCO (PRDM9 short-tract + non-PRDM9 long-tract ~2%) is the mechanism vocabulary for our gene-conversion-like patches.
- **ADD `@noyes2026sd`** (REFRESH_07) on the interlocus gene-conversion clause: trio-resolved interlocus gene conversion is a complementary class of evidence to our pedigree-resolved patches.
- **ADD `@chen2025paraphase`** (REFRESH_07) as MEDIUM on the paralog resolution clause: Paraphase resolves paralog-specific variation in long reads - methodologically aligned with our PHR analysis.
- **ADD `@Tardy2026fshd`** on the C1/4q-10q clause: long-read SV resolution at 4q35/10q26 directly relevant to the chr10q<-chr4q PAN028 event narrative.
- **Existing `@Porubsky2025` is now UPDATED in REFERENCES_v4.bib** (published Nature 643:427-436); add a clarifying sentence acknowledging that the Porubsky 2025 published analysis reports no whole-genome crossover-SV correlation, but our community-constrained PHR exchange is mechanistically distinct (see Contradiction 3 in LITERATURE_REFRESH_v1.md).

**Net effect: +4 strong citations; existing Porubsky2025 metadata updated in v4.**

---

## P10: RPE-1 self analysis + mouse generalisation - line 40

**Currently cited (4 bibkeys)**: `@Francis2025`, `@Zuo2021`, `@bouquet_BhattTERBEvolution2020`, `@Patel2019`.

**Actions:**
- **ADD `@hic3d_Volpe2025RPE1`** on the RPE-1 clause: cite the RPE-1 reference genome assembly that the analysis depends on. HIGH PRIORITY - currently uncited for the cell line itself.
- **ADD `@t2t_Zhang2025macaque`** (REFRESH_15) as LOW priority on the mouse-pipeline-generalises clause: T2T rhesus macaque is the next-generation non-human primate substrate. Skippable if word-count is tight.

**Net effect: +2 (1 high priority Volpe; 1 optional).**

---

## P11: Causal loop, limitations, outlook - line 42

**Currently cited (15 bibkeys)**: `@concerted_evolution_nahr_Arnheim1980`, `@concerted_evolution_nahr_Ohta1984`, `@concerted_evolution_nahr_Charlesworth1994`, `@concerted_evolution_nahr_Hillis1991`, `@concerted_evolution_nahr_Vollger2023`, `@Vollger2023`, `@Xu2025`, `@hic3d_scnanoHiC2023`, `@hic3d_scnanoHiC2_2025`, `@Lalli2025`, `@Sasani2019`, `@Smolka2024`, `@acrocentric_Porubsky2025denovo`, `@Porubsky2025`, `@Logsdon2024`.

**Actions:**
- **ADD `@palsson2025recomb`** on the "short-read recombination maps cannot resolve PHRs and long-read maps are required" clause. SUPPLEMENT to `Sasani2019` and `Smolka2024` (not a REPLACE - those remain valid for their original methodological contexts). This is the **#1 must-add citation** in the entire upgrade plan.
- **ADD `@hic3d_Chen2026HiChew`** (REFRESH_15) as LOW priority on the methods outlook clause: snHiChew is a methodological advance over Dip-C/sperm scHi-C for future validation.
- **ADD `@pangenome_Loegler2025review`** (REFRESH_15) as LOW priority on the pangenome SV framing clause.
- Add **one sentence** acknowledging the Porubsky 2025 no-crossover-SV result (see Contradiction 3 in LITERATURE_REFRESH_v1.md).

**Net effect: +3 net; one supplement (Palsson alongside Sasani/Smolka, not a replace).**

---

## Methods paragraphs (M1..M17)

Most Methods paragraphs are method-specific and citation-light. Targeted upgrades:

### M1 (samples) - line 46
- Currently cited: `@hprc_hprcv2_2025`, `@Nurk2022`.
- **ADD `@logsdon2025hgsvc`** for the population-scale long-read substrate context.

### M3 (wfmash) - line 50
- Currently cited: `@Guarracino2023`.
- No upgrade; the original wfmash citation is correct.

### M4 (IMPG + sampling rate) - line 52
- Currently cited: `@pangenome_graphs_impg_GarrisonGuarracino2023`, `@pangenome_graphs_impg_GuarracinoHeumos2022`, `@pangenome_graphs_impg_IMPG2023`, `@pangenome_graphs_impg_Hickey2024`.
- **ADD `@kaushan2026tracepoints`** (REFRESH_10). PGGB + impg generalisation reference if applicable.

### M6 (PGGB + Jaccard) - line 56
- Currently cited: `@Garrison2024pggb`, `@Guarracino2023`.
- **ADD `@andreace2023pangenome`**. PGGB benchmark independent corroboration.
- **ADD `@heumos2024nfcore`** as the nf-core scaling reference if a Methods note on production deployment is appropriate.

### M11 (heterogeneity tests, F_ST) - line 66
- Currently cited: `@subtel_popgen_hudson1992`, `@subtel_popgen_weir1984`.
- **ADD `@hprc_siren2025`** (Schloissnig 2025; misnamed key, see LITERATURE_REFRESH_v1.md Contradiction 2 note) on the Hudson FST benchmark context.

### M14 (single-cell 3D) - line 72
- Currently cited: `@Tan2018`, `@Xu2025`, `@hic3d_scnanoHiC2023`, `@hic3d_scnanoHiC2_2025`.
- **ADD `@hic3d_kitamura2025`** (REFRESH_09). Contemporary scNanoHi-C reference for the single-cell long-read 3D framing.

### M16 (pedigree untangle) - line 76
- Currently cited: `@Cechova2025`, `@Porubsky2025`.
- **No new ADDs** - `Porubsky2025` is now updated to the published version in `REFERENCES_v4.bib`. Verify `@Cechova2025` has DOI 10.64898/2025.12.14.693655 and PMID 41473289 (REFRESH_11 §3.2 flagged this update).

**Net effect on Methods: +6 (one per M1, M4, M6 [+1 optional], M11, M14).**

---

## Notes for the human reviewer

1. **Zero true REPLACEs**: I evaluated every REFRESH file's "supersedes/replaces" language and concluded each is *additive* in the context of this manuscript (historical citations remain valid for their original observations). The single semi-replace is `@palsson2025recomb` as a *supplement* to `@Sasani2019; @Smolka2024` in P11 - keep all three.

2. **Zero REMOVEs**: No v3 citation was flagged as actively wrong in any REFRESH. Three minor issues exist:
   - `Porubsky2025` was a stale placeholder - **UPDATED in place** in `REFERENCES_v4.bib` (Nature 643:427-436, PMID 40269156).
   - `hprc_siren2025` has the wrong first author label in v3 (the DOI 10.1038/s41586-025-09290-7 / PMID 40702182 is Schloissnig 2025, not Sirén). Bibkey naming bug; do not rename in v4 to avoid breaking existing draft cites; flag for v5 (suggested target name: `hprc_schloissnig2025` or `subtel_popgen_schloissnig2025` - these are aspirational names only, NOT current bibkeys).
   - `subtelstruct_Sholes2022` has wrong authors (REFRESH_02 flagged). Verify before next compile.

3. **Verification checklist** (every recommended ADD bibkey is present in REFERENCES_v4.bib):

```
ADD bibkeys -> v4 presence check (run by integrator):

palsson2025recomb       OK (added by integrator from REFRESH_11)
Salsi2026fshd            OK (REFRESH_03 canonical; salsi2026d4z4t2t and salsi2026t2t are dups, dedup chose Salsi2026fshd)
yang2025chr2fusion       OK (REFRESH_01)
logsdon2025hgsvc         OK (REFRESH_12)
acrocentric_rdna_robertsonian_hartley2026biobank  OK (already in v3 from prior work)
andreace2023pangenome    OK (REFRESH_10)
heumos2024nfcore         OK (REFRESH_10)
kaushan2026tracepoints   OK (REFRESH_10)
hebbar2026marmoset       OK (REFRESH_07 canonical; REFRESH_15's t2t_Hebbar2026marmoset is dup)
degennaro2026ape         OK (REFRESH_05 canonical; REFRESH_15's acrocentric_deGennaro2026apeSD is dup)
jeong2025segdup          OK (REFRESH_12 canonical; REFRESH_13's subtel_popgen_jeong2025 is dup)
porubsky2026chr22q11     OK (REFRESH_13)
hprc_siren2025           OK (v3; carries Schloissnig DOI; rename in v5)
bird2023africa           OK (REFRESH_13)
hic3d_cheng2024          OK (REFRESH_09)
bouquet_KaiserCTCF2025   OK (REFRESH_08)
hic3d_kitamura2025       OK (REFRESH_09)
bouquet_GarnerKASH52023  OK (REFRESH_08)
bouquet_LiuSPDYA2025     OK (REFRESH_08 canonical; REFRESH_04's liu2025spedy is dup)
bouquet_MengSUN1NOA2023  OK (REFRESH_08)
bouquet_JimenezCentromere2025  OK (REFRESH_08)
subtelstruct_Smeds2025nonBDNA  OK (REFRESH_15)
pedigree_Schweiger2024spermNCO OK (REFRESH_15)
noyes2026sd              OK (REFRESH_07)
chen2025paraphase        OK (REFRESH_07)
Tardy2026fshd            OK (REFRESH_03 canonical; tardy2026fshdbeta and tardy2026lrseq are dups)
hic3d_Volpe2025RPE1      OK (REFRESH_15)
t2t_Zhang2025macaque     OK (REFRESH_15)
hic3d_Chen2026HiChew     OK (REFRESH_15)
pangenome_Loegler2025review  OK (REFRESH_15)
Porubsky2025             OK (UPDATED in v4 to published version)
chi2025primate           OK (REFRESH_14; for Ed4 only)
```

All 32 distinct recommended bibkeys exist in `REFERENCES_v4.bib`. The total number of ADD actions across all paragraphs is 37 (some bibkeys appear in multiple paragraphs - e.g. `Salsi2026fshd` appears in both P1 and P4 as the same hedge).

4. **Bib parsing**: `pybtex` parsed `REFERENCES_v4.bib` cleanly (364 entries). Bibkey uniqueness: 0 duplicates by `awk '!seen[$0]++'` on the bibkey list.

5. **Recommended next step**: Build a `NATURE_DRAFT_v2.md` by surgically applying the 37 ADDs and 1 supplement listed above. Estimated word-count delta: +200 words (mostly citations and one hedging clause per contradiction). The bib will not need to grow beyond `REFERENCES_v4.bib`.
