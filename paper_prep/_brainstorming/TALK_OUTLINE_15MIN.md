# 15-MINUTE TALK — Population-scale subtelomeric communities mirror 3D nuclear organisation

**Audience.** Mixed genomics / structural-biology / population-genetics audience.
**Goal.** In 15 min: motivate the question, deliver four anchoring numbers, and leave one mechanistic model + one falsifiable prediction.
**Allocation.** ~13 talk slides + 1 title + 1 acknowledgement = **14 slides**. Each slide ≤ 70 s except #1 (longer setup) and #11 (longer model slide).
**Keying to figures.** Each slide annotated `[F#]` or `[ED#]` referencing `MANUSCRIPT_SKELETON.md` so figure assets and slides share a build. Slides without figure tag draw on existing per-section survey artefacts.

---

## Slide map

| # | Slide title | Figure tag | Headline number | Take-home |
|---|---|---|---|---|
| 1 | Title + cytogenetic-era continuity | — | f7501 to pangenome | Four-decade arc: cytogenetic banding → flow-sorted FISH → cosmid surveys → pangenome |
| 2 | Why subtelomeres? | [F1a] | 18,827 / 15,668 | Telomere-anchored 500 kb captures the inter-chromosomal sharing zone Flint/Mefford described |
| 3 | The dataset and pipeline at a glance | [F1] | 232 / 465 | Reproducible pipeline; 95 % identity, deliberately captures Ambrosini's high-identity peak |
| 4 | 41 arms → 15 communities | [F1c] | 12 of 15 | Leiden recapitulates known biology (D4Z4, PAR1/2, acrocentric-p, f7501, OR4F); UPGMA agrees on 12/15 |
| 5 | Three architectural categories | [F1d] | 8/41 / 34/41 / 7/41 | Quantitative extension of Mefford & Trask's "patchwork" model |
| 6 | Within-community: arm identity is preserved (almost) everywhere | [F2a] | p < 1e-300, 8/9 | Allele closer than paralog in 8/9 multi-arm communities; C7 is the lone reversal (70.5 % paralog closer) |
| 7 | The two-domain model holds genome-wide | [F2b] | 99.7 % of haplotypes | Flint–Mefford gradient in 39/48 arms; piecewise breakpoints 15–445 kb; ITS within 25 kb in 11/19 arms |
| 8 | f7501 reproduced and extended; out-of-Africa from cross-arm exchange | [F2c, F2d] | Fst 0.10–0.15 vs 0 | Mefford & Trask AFR enrichment confirmed; 3 new AFR-enriched arms; population tree recovered |
| 9 | Sequence communities are physical: 6 technologies, one signal | [F3a, F3b] | p = 3.9 × 10⁻⁸⁵ | Hi-C / Pore-C / CiFi / Dip-C / sperm / mouse meiotic — all 14 effect sizes on the within-closer side |
| 10 | The S_all negative control + the flanking paradox | [F3c, F3d] | 0/16 + 1/20 | Non-sharing arms invert the signal; unique-sequence flanks 100 kb away give *stronger* signal than PHRs (multi-mapping ruled out) |
| 11 | Pedigrees catch the events: 92 % land in known communities | [F4a] | 538 / 494 (92 %) | WashU 3-gen T2T: 16 crossover-like + 133 gene-conversion-like + 229 NAHR-acrocentric; CEPH1463 cross-assembler 11 robust features; RPE-1 t(X;10) rediscovery |
| 12 | A causal feedback loop with a mouse-bouquet origin | [ED8a, ED8b] | 4 links | Similarity → 3D proximity → ectopic exchange → similarity; D4Z4-CTCF-lamin model for C1; bouquet hypothesis for initiation |
| 13 | What we *can't* claim (yet) | [ED8c] | ρ = 0.00, N = 32 | Recombination null after confound removal; no human meiotic Hi-C; LCL somatic-exchange caveat; mouse extrapolation honesty |
| 14 | Three predictions, three missing datasets + acks | — | — | LINC complex (SUN1 W151R Hi-C); CTCF density at PHR; long-read recombination map; PARs as the limit case of a continuum |

---

## Slide-by-slide details

### Slide 1 — Title + cytogenetic-era continuity (~1 min)
- One-line title: "Population-scale subtelomeric communities mirror 3D nuclear organisation."
- One-image set-up: f7501 cosmid FISH pattern (Mefford & Trask 2002 Fig. 3 reference) → "this work, 465 near-complete assemblies, four decades later".
- Four-decade narrative arc in two sentences (cytogenetic banding → flow-sorted FISH → cosmid surveys → pangenome). Source: `SURVEY_FRAMING Part 3`.

### Slide 2 — Why subtelomeres? [F1a] (~50 s)
- "The most distal 500 kb of every chromosome end shares sequence with non-homologous chromosomes — to a degree that has been invisible at single-genome resolution."
- Visual: genome-wide identity heatmap (`p_genome_wide_identity_heatmap.pdf`).
- Four prior anchors named (Flint 1997 two-domain; Mefford & Trask 2002 patchwork; Linardopoulou 2005 hot spots; Ambrosini 2007 duplicons).

### Slide 3 — Dataset & pipeline [F1, ED1a] (~50 s)
- 232 individuals × 2 haplotypes = 464, plus CHM13 = 465 near-complete assemblies → 12,649 classified contigs → 18,827 telomere-anchored 500 kb flanks → 15,668 PHRs (`SURVEY_01 §1`).
- Pipeline schematic (5 stages); 95 % identity threshold deliberately captures the high-identity Ambrosini peak.
- One sentence honesty: "the 91 % older-exchange peak is missed by design".

### Slide 4 — 41 arms → 15 arm-level communities [F1c] (~70 s)
- 41×41 arm Jaccard heatmap with Leiden 15-community blocks (silhouette 0.347).
- Annotate the named communities: C1 D4Z4 (chr4q+chr10q); C2 chr10p+chr18p (Linardopoulou); C3 f7501 (chr3q+chr7p+chr9q+chr11p+chr16q+chr19p); C4 chr7q+chr12q (private, *new*); C5 RPL23A/WASH/DDX11L; C7 acrocentric-p; C11 OR4F; C14 PAR2; C15 PAR1.
- Method robustness (1 line): UPGMA agrees on 12/15.

### Slide 5 — Three architecture categories [F1d] (~50 s)
- Per-arm bar chart split into three blocks.
- 8/41 homogeneous (chr10p, chr12p, chr17q, chr18p, chr1q, chr20q, chr21q, chr9q — 0 % cross-arm).
- 34/41 polymorphic (e.g. chr11_p 62.7 % cross-arm — outlier).
- 7/41 fully interchangeable (chrY_p → chrX_p; chrY_q → chrX_q; chr15p/chr21p/chr13p/chr22p → chr14p; chr10q → chr4q).
- "This is Mefford & Trask's qualitative patchwork model, quantified."

### Slide 6 — Allele vs paralog: arm identity is preserved [F2a] (~70 s)
- Boxplot of allele vs paralog Jaccard distance per community + reversal table.
- 5,946 paired observations; allele closer in 8/9; combined Wilcoxon p < 1e-300.
- C7 reversal is the population-scale signature of complete acrocentric homogenisation: silhouette = −0.029, conversion score = 1.000.
- Optional 1-line: D4Z4 (C1) is the autosomal weak spot — 44.2 % paralog closer.

### Slide 7 — Two-domain model genome-wide [F2b] (~70 s)
- Per-arm gradient + breakpoint composite.
- 39/48 arms significant gradient; 39/41 testable arms prefer the two-segment fit.
- 99.7 % of individual haplotype sequences show the within-sequence gradient.
- Breakpoints arm-specific: chr22q 15 kb → chr3p 445 kb; original FISH-era arms (chr4p 70 kb, chr4q 50 kb, chr16p 295 kb) confirm.
- Internal (TTAGGG)n islands within 25 kb of breakpoint on 11/19 testable arms.

### Slide 8 — Population history from subtelomeres [F2c, F2d] (~70 s)
- f7501 reproduction in two lines: fixed sites confirmed (chr3_q 91.8 %, chr19_p 90.5 %, chr15_q 85.6 %); chr16_q AFR OR = 17.4, p = 6.6e-27; *new*: chr8_p, chr16_p, chr9_q AFR-enriched; chr2_q SAS, chr6_p AMR.
- Cross-arm-affinity Fst panel: AFR vs non-AFR 0.10–0.15; non-AFR/non-AFR ≈ 0.
- Out-of-Africa topology recovered from cross-arm frequencies (the 24-year-old hedged Mefford & Trask 2002 hypothesis, now answered).

### Slide 9 — Three-dimensional convergent evidence [F3a, F3b] (~90 s, key slide)
- HG002 Pore-C contact matrix grouped by community (top half).
- Convergent-evidence forest plot (bottom half): 14 effect sizes from Hi-C / Pore-C / CiFi / Dip-C / sperm / mouse meiotic, all on the within-community-closer side.
- Headline numbers (memorise): HG002 Pore-C **3.9e-85**; sperm 60 % closer **3.9e-51**; GM12878 Dip-C Mantel ρ = 0.296; CHM13 per-pair ρ = 0.674; mouse zygotene ρ = 0.715.
- "Six technologies. Two species. One signal."

### Slide 10 — Two diagnostic controls [F3c, F3d] (~70 s)
- S_all negative control: 7 chromosome arms with zero PHR signal pooled. **0/16 GM12878 + 1/20 sperm cells fall below W/B = 1**; instead 11 % / 40 % *farther* — sequence sharing is necessary, not incidental.
- Flanking paradox: 100 kb centromere-ward of PHR is unique sequence, yet HG002 flanking B/W = 0.002 (vs PHR 0.027) — 13× *stronger*. Multi-mapping is ruled out; the chromosomal-domain effect extends past the duplicated region.
- Confound exclusions strengthen the signal: HG002 Mantel 0.66 → 0.79; HG02148 0.15 → 0.72.

### Slide 11 — Pedigrees + cross-species generalisation [F4a, F4d] (~90 s, key slide)
- WashU PAN028 maternal hap1 untangle ribbon (right panel of `SURVEY_14 §6 P-1`) annotated with chr22p ← chr14p (15.4 kb gene-conversion), chr3q ← chr7p (28 kb crossover), chr1p ← chr8p (14.9 kb crossover).
- Headline: **538 inter-chromosomal exchanges, 92 % within Leiden community**; 16 crossover-like, 133 gene-conversion-like, 229 acrocentric NAHR signatures.
- CEPH1463 corner inset: 11 robust cross-assembler features; chr10/chr18 (Linardopoulou) detected in NA12877 paternal *and* NA12878 maternal independently.
- Mouse zygotene PHR-pair scatter (ρ = 0.715) — same correlation, different species, different cell type, different architecture.
- One-line on RPE-1: pipeline run on a single individual rediscovers t(X;10) translocation from sequence alone.

### Slide 12 — A causal feedback loop with a mouse-bouquet origin [ED8a, ED8b] (~90 s)
- Four-link cycle schematic: similarity → 3D proximity → ectopic exchange → similarity.
- Coloured arrows by support level (3 of 4 supported by present data; 1 inferred).
- D4Z4-CTCF-lamin tethering model for C1: peripheral radial 0.732, 0–15 kb sharing peak at D4Z4, DUX4L median 22 vs outliers 0–2 (p = 5.3e-6).
- Bouquet hypothesis as initiation node (Zuo 2021 mouse meiotic; ~20 % chromosome-end alignment; median PHR 105 kb fits inside one leptotene loop).
- Honest caveat (1 line): causal direction undecidable from present data.

### Slide 13 — What we can't claim (yet) [ED8c] (~70 s)
- Top 4 caveats:
  1. **No human meiotic Hi-C exists** — somatic 3D is interpreted as a Rabl residual, indirectly.
  2. **Recombination null:** ρ = −0.43 across 39 arms collapses to **ρ = 0.00, p = 0.98** when 7 short-read-confounded acrocentric/PAR arms are excluded.
  3. **LCL somatic exchange:** chr4_q/chr10_q (C1) signal could include in-culture exchange; source-stratified validation needed.
  4. **95 % threshold misses the 91 % older-exchange peak** by design (Ambrosini 2007 bimodality).
- Tone: each caveat names the experiment that would resolve it; this is the next-experiment slide, not a hedging slide.

### Slide 14 — Three predictions + acknowledgements (~50 s + 20 s acks)
- LINC complex requirement: re-run Zuo 2021 SUN1 W151R Hi-C through the community framework.
- CTCF density at PHR boundaries: Gershman 2022 T2T-CHM13 ENCODE + Stergachis Fiber-seq (39/46 telomeres).
- Long-read recombination map: re-test the cross-arm/cM correlation at confounded arms.
- Closing line: **"PARs are the special case. PHRs are the general case. The genome has 24 chromosome ends; we can now study them as a system."**
- Acknowledgements: HPRC, Cechova et al. 2025, Porubsky et al. 2025, Xu et al. 2025, Zuo et al. 2021, the ~30 cytogenetics-era papers without which this work would have no scaffold.

---

## Backup slides (have ready, do not show by default)

- **B1.** Cell-cycle modulation (RPE-1 mitotic 3× stronger global W/B but 1.4× weaker per-pair ρ) — for question on synchronisation.
- **B2.** Compartment identity at tips (68 % A but mean e1 = +0.007 — weak) — for question on A/B framework.
- **B3.** Per-resolution + per-window mouse panel (1/2/4 Mb × 5 resolutions) — for question on parameter sensitivity.
- **B4.** Mouse private-pair → human syntenic-net (`SURVEY_08 §6 T-6`) — for question on mouse-human comparability.
- **B5.** Pedigree donor-asymmetry (PAN027 chr22p:h2 dominant donor; CEPH1463 chr10/chr18 maternal/paternal asymmetry) — for question on parent-of-origin effects.

---

## Recall numbers (memorise before delivering)

| Number | Anchor (one-liner if asked) |
|---|---|
| 232 / 465 | individuals / near-complete assemblies (HPRCv2) |
| 18,827 / 15,668 | flanks / signal-bearing PHRs |
| 41 → 15 | arms → arm-level communities |
| 50 | sequence-level communities |
| 15.9 % / 11.1 % | cross-arm at arm / sequence level |
| 47.5 % | maximum discordance (chr22_q, C6) |
| 92 % (494/538) | WashU within-Leiden patches |
| 0.027 / 3.9e-85 | HG002 Pore-C B/W / p-value |
| 0.401 / 60 % closer / 3.9e-51 | sperm W/B |
| 0/16 + 1/20 | S_all reversal (GM12878 + sperm) |
| 13× | HG002 flanking B/W stronger than PHR |
| ρ = 0.715 | mouse zygotene per-PHR-pair |
| ρ = −0.43 → 0.00 | recombination null after confound removal |
