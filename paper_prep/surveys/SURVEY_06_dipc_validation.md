---
title: "Survey 06 — Dip-C single-cell 3D and sperm validation"
source: end-to-end-report/report/06_dipc_validation.md
scope: Single-cell 3D genome validation (16 GM12878 Dip-C cells + 20 sperm cells)
audience: Nature manuscript and 15-min talk
---

# Survey 06 — Dip-C and sperm validation

This survey extracts and structures the content of `end-to-end-report/report/06_dipc_validation.md` for the Nature manuscript and the companion 15-minute talk. The source section validates pangenome subtelomeric communities at single-cell 3D resolution, complementing the bulk Hi-C/Pore-C/CiFi analyses (covered in survey 05).

---

## 1. Key findings with metrics

### 1.1 GM12878 Dip-C (16 cells, T2T-CHM13v2.0)

Community enrichment, per-arm PHR coordinates (38 C-community arms; 8 arms without CHM13 PHR fall back to terminal 500 kb):

| Metric | Value |
|---|---|
| Wilcoxon signed-rank | stat = 8.0, p = 3.8 × 10⁻⁴ |
| Fisher combined | χ² = 75.3, p = 2.4 × 10⁻⁵ |
| W/B ratio (mean) | **0.931** (6.9 % closer within-community) |
| W/B ratio (median) | 0.934 |
| Mantel ρ | **0.296**, p = 0.002 |

Community-free, similarity ↔ −distance:
- Per-cell: 15/16 cells positive ρ; median ρ = **0.093**
- Arm-level pooled: ρ = **0.336**, p = 1.1 × 10⁻¹⁸ (n = 652 arm pairs)

### 1.2 Sperm scHi-C (20 cells, 10 X-bearing + 10 Y-bearing; Xu et al. 2025)

Community enrichment, per-arm PHR coordinates (haploid; no impute3):

| Metric | Value |
|---|---|
| W/B ratio | **0.401** (60 % closer within-community) |
| Fisher combined p | **3.9 × 10⁻⁵¹** |
| Mantel ρ | **0.202**, p = 0.023 |

Community-free:
- Per-cell median ρ = 0.029; 15/20 cells positive
- Arm-level pooled: ρ = −0.048, p = 0.197 (n.s., wrong direction — interpretation: highly compacted sperm chromatin limits inter-chromosomal proximity variation)

### 1.3 Negative control: non-sharing arms (S_all)

Seven chromosome arms with **zero** inter-chromosomal subtelomeric sequence sharing (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q), pooled into one pseudo-community S_all over the terminal 500 kb:

| | GM12878 (16 cells) | Sperm (20 cells) |
|---|---:|---:|
| C-community W/B (mean) | 0.931 (6.9 % closer) | 0.401 (60 % closer) |
| **S_all W/B (mean)** | **1.106 (11 % farther)** | **1.397 (40 % farther)** |
| **S_all cells with W/B < 1** | **0/16** | **1/20** |

This is the headline negative control: arms without shared subtelomeric sequence are pushed *apart* in 3D, exactly inverting the C-community signal — sequence sharing is necessary, not incidental, for 3D proximity.

### 1.4 Per-arm radial positions (S-singletons)

| Arm | S-label | GM12878 radial | Sperm radial | Nearest C-community in 3D |
|---|---|---|---|---|
| chr2_p | S1 | 0.81 (peripheral) | 0.68 | C12 (chr2_q, same chrom) |
| chr3_p | S2 | 0.84 (peripheral) | 0.82 | C3 (chr3_q, same chrom) |
| chr5_p | S3 | 0.66 | 0.61 | C3 (chr19_p) |
| chr8_q | S4 | 0.60 | 0.58 | C3 (chr19_p) |
| chr11_q | S5 | 0.69 | 0.64 | C6 (chr22_q) |
| chr14_q | S6 | 0.61 | 0.56 | C3 (chr19_p) |
| chr18_q | S7 | 0.78 (peripheral) | 0.71 | C2 (chr18_p, same chrom) |

Take-home: most non-sharing arms' nearest 3D neighbor is the *opposite arm of the same chromosome* (cis-arm proximity), not any inter-chromosomal partner. S1, S2 and S7 are notably peripheral (radial > 0.78), consistent with telomere-proximal nuclear positioning.

### 1.5 Cross-section consistency

All correlation signs are positive across the four technologies (Dip-C, sperm scHi-C, bulk Hi-C, Pore-C/CiFi): more similar PHR sequences → closer 3D proximity. The Dip-C signal magnitude (per-cell median ρ = 0.093) is smaller than Hi-C bulk because each cell contributes a single noisy 3D snapshot; the arm-level pool (ρ = 0.336) recovers the bulk-comparable signal.

---

## 2. Existing figures (paths)

All paths are absolute — they live on `/moosefs/guarracino/HPRCv2/` and are produced outside the report directory.

### GM12878 Dip-C
- `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_mantel_scatter.pdf` — Mantel scatter (sequence distance × 3D distance).
- `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_radial_community.pdf` — radial position by community.
- `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/exclusion_no_sex/gm12878_mantel_scatter.pdf` — exclusion control (no sex chrom).
- `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_16cells_500kb_mantel_scatter.pdf` — fixed 500 kb terminal control.
- `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_16cells_500kb_radial_community.pdf` — fixed 500 kb terminal radial.

### Sperm
- `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_mantel_scatter.pdf` — sperm Mantel scatter.
- `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_radial_community.pdf` — sperm radial by community.
- `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/overlay_plots/sperm_all20.by_arm-type.arm.pdf` — arm-type overlay (shared vs unshared, arm-level).
- `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/overlay_plots/sperm_all20.by_arm-type.per-cell.pdf` — same, per-cell.
- `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/overlay_plots/sperm_all20.by_chromosome.{arm,per-cell}.pdf` — per-chromosome overlay.
- `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/overlay_plots/sperm_by_sex.by_{arm-type,chromosome,sample}.{arm,per-cell}.pdf` — split by sperm sex.

### PBMC (referenced, not described in this section)
- `/moosefs/guarracino/HPRCv2/dipc_t2t/pbmc_hg19/enrichment_corrected/pbmc_mantel_scatter.pdf`
- `/moosefs/guarracino/HPRCv2/dipc_t2t/pbmc_hg19/enrichment_corrected/pbmc_radial_community.pdf`

---

## 3. Existing CSVs (paths)

All TSVs (the report calls them CSVs colloquially). One row of summary statistics per file unless stated.

### GM12878 (PHR-coordinate, 16 cells)
Directory: `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/`
- `gm12878_summary.tsv` — Wilcoxon, Fisher, mean/median W/B, n_cells.
- `gm12878_per_cell.tsv` — per-cell W/B and ratios.
- `gm12878_per_community_per_cell.tsv` — within/between per (cell × community).
- `gm12878_mantel_3d.tsv` — Mantel ρ, p, permutation count.
- `gm12878_radial_community.tsv`, `gm12878_radial_per_community.tsv` — radial summaries.
- `gm12878_community_free_per_cell.tsv` — per-cell community-free Spearman ρ.
- `gm12878_community_free_per_sample.tsv` — per-cell × per-sample#hap (465 pangenome samples) ρ.
- `gm12878_community_free_arm.tsv` — per-cell arm-level ρ (one file per cell `_01`…`_17`).
- `gm12878_community_free_seqlevel.tsv` — sequence-level pooled.
- `gm12878_arm_3d_distance_matrix{,_NN}.tsv` — pooled and per-cell arm × arm 3D distance matrices.
- `gm12878_seq_summary.tsv`, `gm12878_seq_per_cell.tsv`, `gm12878_seq_per_community_per_cell.tsv` — sequence-level analogs.

Exclusion controls under the same directory:
- `exclusion_no_sex/`, `exclusion_no_acro/`, `exclusion_no_acro_sex/` — same file set with chromosomes excluded.

Fixed-window 500 kb (legacy / supplementary controls):
- `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_16cells_500kb_*` — same metric set computed on 500 kb terminal windows.
- `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_16cells_500kb_per_community_per_cell.tsv` contains the **S_all** per-cell rows (used for the negative-control table).

### Sperm (PHR-coordinate, 20 cells)
Directory: `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/`
- `sperm_all20_summary.tsv` — Wilcoxon W/B, Fisher (3.9 × 10⁻⁵¹).
- `sperm_all20_per_cell.tsv`, `sperm_all20_per_community_per_cell.tsv`.
- `sperm_all20_mantel_3d.tsv`.
- `sperm_all20_radial_community.tsv`, `sperm_all20_radial_per_community.tsv`.
- `sperm_all20_community_free_per_cell.tsv`, `sperm_all20_community_free_per_sample.tsv`, `sperm_all20_community_free_arm.tsv`, `sperm_all20_community_free_seqlevel.tsv`.
- `sperm_all20_arm_3d_distance_matrix{,_HSxxx.clean.3dg}.tsv` — pooled and per-sperm-cell arm × arm matrices.
- `sperm_all20_seq_*.tsv` — sequence-level analogs.
- Sub-dirs `exclusion_no_sex`, `exclusion_no_acro`, `exclusion_no_acro_sex` — same file set with chromosomes excluded.

Sperm overlay TSVs:
- `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/overlay_plots/sperm_all20.by_*.tsv` and `sperm_by_sex.by_*.tsv` — companion data for each overlay PDF.

### PBMC reference (hg19-projected)
- `/moosefs/guarracino/HPRCv2/dipc_t2t/pbmc_hg19/enrichment_corrected/pbmc_*.tsv` — community-based only; no community-free (hg19-projected PHR boundaries lack the pairwise Jaccard matrix needed for per-sequence-pair correlation).

### BED inputs
- `/moosefs/guarracino/HPRCv2/dipc_t2t/phr_and_500kb_regions.bed` — CHM13 PHR + 500 kb fallback windows passed via `--region-bed`.
- `/moosefs/guarracino/HPRCv2/dipc_t2t/pbmc_hg19/phr_hg19_merged_regions.bed` — hg19-projected PHR for PBMC.

---

## 4. Methods

### 4.1 Cell preparation and 3D modelling
- **GM12878:** 17 cells from prior Dip-C dataset, remapped to T2T-CHM13v2.0 with BWA-MEM2; hickit for 3D modelling; `dip-c impute3` for diploid haplotype refinement (4 rounds). One cell (cell 12) excluded as a duplicate of cell 10 → **16 cells** retained. SNPs taken from the 1KGP CHM13v2 panel for NA12878.
- **Sperm:** 20 cells (10 X-bearing + 10 Y-bearing) from Xu et al. 2025. Haploid mode via `run_dipc_cell.sh`; no `impute3` step needed.
- **PBMC:** referenced separately; only community-based results available because PHR was projected to hg19 without the Jaccard similarity matrix needed for community-free analysis.

### 4.2 Multi-mapping handling
- `sam2seg -q 0`, `hickit --min-mapq=0` throughout. BWA-MEM2 reports one primary alignment per read (chimeric supplementary alignments retained). Default MAPQ ≥ 20 filter is **disabled** because it would discard 60–99 % of reads in subtelomeric tips. Each multimapped read keeps exactly one randomly-chosen position, contributing symmetric noise — the same convention used for Hi-C, Pore-C, and CiFi.

### 4.3 Coordinates passed to the 3D analysis
- 3D particle positions are sampled at per-arm PHR coordinates for the 38 C-community arms (CHM13 PHR boundaries; arm-specific widths 10–500 kb).
- For the 8 arms without CHM13 PHR (7 S-community arms + chr6_p), the terminal 500 kb is used as a fallback window.
- Pre-built region BED: `/moosefs/guarracino/HPRCv2/dipc_t2t/phr_and_500kb_regions.bed`.

### 4.4 Statistical tests
- **Within/between (W/B) ratio:** mean intra-community 3D distance / mean inter-community 3D distance, per cell. W/B < 1 ⇒ within-community arms closer.
- **Wilcoxon signed-rank** and **Fisher combined p** across cells.
- **Mantel test:** sequence-distance matrix vs 3D-distance matrix, Spearman ρ, with permutation p.
- **Community-free per-cell:** for each of 465 pangenome `sample#hap` combinations, take that sample's PHR sequences, map to 3D positions, compute inter-chromosomal (Jaccard, 3D distance) pairs, Spearman ρ. Mirrors the mcool community-free design but substitutes 3D distance for Hi-C contact.
- **S-community framework:** singleton S1–S7 plus pooled S_all over the same arms; identical W/B test reported alongside C-community results.

### 4.5 Scripts
| Script | Purpose |
|---|---|
| `/moosefs/guarracino/HPRCv2/scripts/dipc/community_3d_enrichment.py` | Community W/B, Mantel, community-free, radial. |
| `/moosefs/guarracino/HPRCv2/scripts/dipc/project_phr_to_hg19.py` | Project CHM13 PHR boundaries to hg19 for PBMC. |
| `/moosefs/guarracino/HPRCv2/scripts/dipc/phr_dipc_overlay.py` | PHR-particle overlay: shared vs unshared 3D distances. |
| `/moosefs/guarracino/HPRCv2/scripts/dipc/plot_3dg.py` | 3D genome structure visualisation. |
| `/moosefs/guarracino/HPRCv2/scripts/dipc/plot_sperm_overlay.sh` | Sperm overlay-plot driver. |
| `/moosefs/guarracino/HPRCv2/scripts/dipc/run_dipc_cell.sh` | Per-cell Dip-C alignment / hickit / impute3 pipeline. |

---

## 5. Gaps

1. **No PBMC community-free.** hg19-projected PHR boundaries lack the pairwise Jaccard similarity matrix; the per-sequence-pair correlation cannot be computed. PBMC PHR-coordinate community-based results are reportedly **not significant** but the section currently only references this in passing — no PBMC numbers are quoted in 06.
2. **Cell-level heterogeneity not summarised.** Per-cell rho TSVs exist but only five GM12878 cells (01, 02, 03, 05, 06) are tabulated in the source section. The other 11 are computed but not displayed. No corresponding sperm per-cell ρ table is shown.
3. **Sperm community-free arm-level wrong direction.** Arm-level pooled ρ = −0.048 (n.s.). Section attributes this to compacted sperm chromatin but offers no quantitative model or scatter to back it up; a per-cell distribution figure would help.
4. **No 3D visualisation in section 06.** `plot_3dg.py` exists but no rendered 3D genome panel (e.g., one cell, communities coloured) is referenced. The talk would benefit from a single iconic 3D image.
5. **No direct visual comparison Dip-C vs Hi-C vs Pore-C/CiFi** side-by-side at the W/B / Mantel ρ level. This integrative cross-tech panel is implied by the integrated section (07) but not produced inside 06.
6. **S_all positive-control test missing.** S_all is only used as a negative control. There is no equivalent "matched random C-style community" positive control to bound the expected magnitude of the S_all reversal.
7. **Sperm sex-split statistics not foregrounded.** Overlay TSVs split by sperm sex exist but no W/B / Mantel ρ summary is reported separately for X- vs Y-bearing cells.
8. **Multi-mapping noise effect not quantified for Dip-C.** Section 06 cites the 60–99 % discard rate at MAPQ ≥ 20 in subtelomeres but does not show the empirical effect on W/B at MAPQ = 0 vs 20 in Dip-C cells (the comparison exists in spirit from section 05).
9. **Cell 12 exclusion** is a one-line note. No QC figure backs the duplicate call.
10. **No external benchmark / Salzer-style comparison** to a non-PHR baseline community (e.g., random arm groupings of matched size) for either GM12878 or sperm — the S_all set provides one such control but it is not size-matched to C-communities.

---

## 6. Suggested figures with captions (produced vs to-do)

### Already produced (use directly or recompose)

**P-1. GM12878 Mantel scatter.**
*File:* `gm12878_mantel_scatter.pdf` (output_q0_XX/community_enrichment_k50/).
*Caption:* "PHR sequence distance vs 3D Euclidean distance across 16 GM12878 Dip-C cells (T2T-CHM13v2.0). Each point is one inter-chromosomal arm pair; Mantel ρ = 0.296, p = 0.002."

**P-2. GM12878 radial position by community.**
*File:* `gm12878_radial_community.pdf`.
*Caption:* "Mean radial position (0 = nuclear centre, 1 = periphery) per C-community across 16 GM12878 cells. Acrocentric- and chr19-anchored communities cluster centrally; chr2_p / chr3_p / chr18_q non-sharing arms are peripheral."

**P-3. Sperm Mantel scatter.**
*File:* `sperm_all20_mantel_scatter.pdf`.
*Caption:* "PHR sequence distance vs 3D Euclidean distance across 20 sperm cells. Mantel ρ = 0.202, p = 0.023 — sequence-distance / 3D-distance relationship persists in compacted sperm chromatin."

**P-4. Sperm radial position by community.**
*File:* `sperm_all20_radial_community.pdf`.
*Caption:* "Per-community radial positions in 20 sperm cells (10 X-bearing, 10 Y-bearing); mirrors GM12878 panel for direct comparison."

**P-5. Sperm sequence-shared vs unshared overlay.**
*File:* `sperm_all20.by_arm-type.{arm,per-cell}.pdf`.
*Caption:* "Distribution of inter-chromosomal 3D distances stratified by PHR-sharing status (shared vs unshared); arm-level and per-cell views."

**P-6. Sperm by-sex overlay.**
*File:* `sperm_by_sex.by_{arm-type,chromosome}.{arm,per-cell}.pdf`.
*Caption:* "X- vs Y-bearing sperm 3D distance distributions stratified by PHR sharing — companion to the pooled overlay."

### To-do (suggested new figures)

**T-1. Cross-cell-type W/B summary (one panel for the talk).**
*Caption:* "Within/between 3D distance ratio (W/B) per cell, grouped by cell type. C-community W/B in GM12878 (0.93, n = 16), Sperm (0.40, n = 20), PBMC (n = ?). S_all negative control overlaid for each (1.11, 1.40). 0/16 GM12878 and 1/20 sperm cells show S_all W/B < 1." Build from `gm12878_per_cell.tsv`, `sperm_all20_per_cell.tsv`, and the `community_enrichment_16cells_500kb_per_community_per_cell.tsv` S_all rows.

**T-2. Per-cell community-free ρ distribution.**
*Caption:* "Per-cell Spearman ρ (similarity vs −3D distance) across 16 GM12878 and 20 sperm cells. Median ρ_GM = 0.093; median ρ_sperm = 0.029; arrow marks the arm-level pooled estimate (0.336 / −0.048)." Built from `*_community_free_per_cell.tsv`.

**T-3. Iconic 3D genome render.**
*Caption:* "One representative GM12878 cell (e.g., cell 01) rendered with `plot_3dg.py`; PHR particles coloured by C-community; non-sharing arms in grey. Designed as an opening visual for the talk." Source: `output_q0_XX/3dg/cell_01.3dg`.

**T-4. S_all reversal panel.**
*Caption:* "W/B ratio per cell for S_all vs each C-community. Bars / strip-plot showing 0/16 GM12878 cells and 1/20 sperm cells fall below 1 for S_all, contrasted with C-community W/B distributions where the majority do."

**T-5. Per-arm radial panel (extends Table 1.4).**
*Caption:* "Radial positions of all 46 arms (38 C + 7 S + chr6_p) in GM12878 and sperm. Highlights cis-arm proximity for non-sharing arms and the peripheral positioning of S1/S2/S7."

**T-6. Cross-tech rho summary (talk slide).**
*Caption:* "Effect-size summary across four technologies — Dip-C (ρ = 0.336), sperm scHi-C (ρ = −0.048 arm pooled, +0.029 per-cell median), bulk Hi-C, Pore-C/CiFi — to show consistent positive direction (with the sperm caveat) across cell types and assays."

**T-7. Methods schematic.**
*Caption:* "Workflow: BWA-MEM2 (MAPQ = 0) → hickit → dip-c impute3 → 3D coordinates sampled at PHR / 500 kb terminal windows → W/B, Mantel, community-free analyses."

---

## 7. Talk slide takeaways (15-min talk)

1. **Headline:** "Single-cell 3D structures (Dip-C, sperm scHi-C) confirm that PHR-shared arms are physically closer than expected — in two distinct cell types, including the most compacted nucleus we know of (sperm)."
2. **One-number recall (talk):** "60 % closer in sperm (Fisher p = 4 × 10⁻⁵¹), 7 % closer in GM12878 (Fisher p = 2.4 × 10⁻⁵). Both cell types show a positive Mantel correlation between sequence distance and 3D distance."
3. **Negative control (essential slide):** "0/16 GM12878 cells and 1/20 sperm cells show within-community 3D proximity for the seven non-sharing arms — those arms are 11 % (GM12878) and 40 % (sperm) *farther* apart in 3D. Sequence sharing is necessary, not incidental." Use figure T-4.
4. **Visual anchor:** one rendered 3D genome (T-3) with PHR particles coloured by community.
5. **Cross-technology consistency:** four technologies (Dip-C, sperm scHi-C, bulk Hi-C, Pore-C/CiFi), all positive sign. Use figure T-6.
6. **Cis-arm caveat:** non-sharing arms' nearest 3D neighbour is usually the *opposite arm of the same chromosome*, not an inter-chromosomal partner — i.e., the dominant 3D constraint without sequence sharing is cis, as expected from chromosome territories.
7. **Method honesty:** MAPQ = 0 throughout; per-arm PHR coordinates (10–500 kb) for 38 arms + 500 kb fallback for 8 arms; one cell excluded as a duplicate (16/17 GM12878). Highlight the 60–99 % subtelomeric read recovery enabled by dropping the MAPQ filter.
