# SURVEY 05 — Hi-C / Pore-C 3D validation

Survey of `end-to-end-report/report/05_hic_validation.md` for the Nature manuscript and 15-min talk on
subtelomeric sequence sharing across HPRCv2.

Source document: `end-to-end-report/report/05_hic_validation.md` (lines 1–689).

---

## 1. Key findings with metrics

The chapter establishes that arms placed in the same Leiden subtelomeric community by sequence sharing
also contact each other more often in 3D, and shows this with three orthogonal statistical tests on
8 datasets (CHM13 + 5 diploid HPRC Hi-C, HG002 Pore-C, HG002 CiFi) plus 3 RPE-1 datasets and dozens of
control comparisons. All major numbers (test statistics, p-values, ARI, B/W) are gathered below.

### 1.1 Community enrichment — within/between contact ratio (PHR, 50 kb)
339 bootstrap tests across 7 samples + 7 global Mann–Whitney tests + 7 Mantel tests. **All 8
datasets significant (all p < 0.01)**:

| Sample / Tech | B/W ratio | Global p |
|---|---|---|
| CHM13 Hi-C | 0.071 | 6.0e-18 |
| HG002 Hi-C | **0.027** (strongest) | 4.0e-66 |
| HG002 Pore-C | 0.056 | 3.9e-85 |
| HG002 CiFi | 0.036 | 2.0e-74 |
| HG00658 Hi-C | 0.056 | 7.6e-12 |
| HG02148 Hi-C | 0.050 | 9.1e-05 |
| HG02559 Hi-C | 0.074 | 9.4e-03 |
| NA19036 Hi-C | 0.049 | 1.9e-07 |

### 1.2 Mantel test (arm-level Jaccard similarity × Hi-C contact, PHR, 50 kb)
**7 / 8 datasets significant** (HG02148 borderline, p=0.085):

| Sample / Tech | Mantel ρ | p | n shared arms |
|---|---|---|---|
| CHM13 Hi-C | 0.656 | <1e-4 | 38 |
| HG002 Hi-C | 0.657 | <1e-4 | 41 |
| HG002 Pore-C | 0.486 | <1e-4 | 41 |
| HG002 CiFi | 0.308 | <1e-4 | 41 |
| HG02559 Hi-C | 0.397 | 4.0e-4 | 37 |
| HG00658 Hi-C | 0.276 | 9.3e-3 | 37 |
| NA19036 Hi-C | 0.266 | 0.018 | 34 |
| HG02148 Hi-C | 0.152 | 0.085 (ns) | 37 |

### 1.3 Independent Hi-C community detection — Adjusted Rand Index (50 kb)
All 8 samples have **ARI > 0** vs the 15 sequence-defined communities:

| Sample / Tech | ARI vs seq |
|---|---|
| CHM13 Hi-C | 0.539 |
| HG002 Hi-C | 0.296 |
| HG02559 Hi-C | 0.128 |
| HG00658 Hi-C | 0.133 |
| HG02148 Hi-C | 0.123 |
| NA19036 Hi-C | 0.165 |
| HG002 Pore-C | 0.264 |
| HG002 CiFi | 0.056 |

### 1.4 Per-arm-pair Spearman ρ (mean Jaccard × Hi-C contact)
All 8 samples significant; ρ ≈ 0.66 in best (CHM13, HG002 Hi-C), down to 0.16 in HG02148
(p = 4.0e-5).

### 1.5 Community-free per-sequence-pair Spearman ρ (PHR, 50 kb / 10 kb)
**All 8 datasets significant; p < 3.7e-24 throughout.** Hi-C ρ = 0.64–0.83 at 10 kb (highest in
NA19036 0.83, HG02148 0.81). Pore-C ρ ≈ 0.38, CiFi ρ ≈ 0.19. Includes RPE-1 (Async CiFi 0.32, Async
Pore-C 0.43, Mitotic CiFi 0.18) and HG002 Hi-C with 2 544 sequence pairs.

### 1.6 Multi-resolution robustness (5 kb–100 kb)
B/W ratios remain well below 1.0 at all 5 mcool resolutions for every sample. HG002 Hi-C: 0.029 / 0.028 / 0.028 / 0.027 / 0.021 (5–100 kb). Confirms resolution invariance.

### 1.7 Acrocentric / sex / strong-community exclusion controls
Mantel ρ and community-free ρ **persist or strengthen** in every exclusion. HG002 Hi-C Mantel: 0.657
→ 0.790 (no acro pq + sex). HG02148 Hi-C: 0.152 → 0.720. CHM13: 0.656 → 0.850. The signal is
**not** driven by acrocentric nucleolar association or PAR.

### 1.8 Flanking control (100 kb centromere-ward of PHR)
Flanking is unique sequence (no multi-mapping). Mantel ρ comparable to PHR (CHM13 0.522 vs 0.656;
HG002 0.520 vs 0.657). All datasets significant at all resolutions; sequence-level community
enrichment 27/27 multi-arm communities significant for HG002 Hi-C (BH q<0.05).

### 1.9 p-arm vs q-arm 3D enrichment
Mann–Whitney U (q-dominated vs p-dominated communities) U=812, **p=0.90 (ns)**. Median enrichment:
q-dom 12.0×, p-dom 10.2×, mixed 7.9×. No systematic p/q asymmetry.

### 1.10 O/E inter-chromosomal normalization
After O/E normalization (controlling for chromosome size & marginal contact), within-community
contacts are still **8.6× to 34.4×** higher than between-community contacts (HG002 Hi-C 34.4×,
CHM13 12.9×, HG02148 17.3×, NA19036 19.2×, HG002 Pore-C 8.6×).

### 1.11 RPE-1 cell-type validation (HPRC communities)
3 RPE-1 datasets, all 5 resolutions, 92 arms (diploid), HPRC 15-community labels.
- Community enrichment: all **p < 1.6e-65**; 50 kb B/W = 0.024 (Async CiFi), 0.031 (Async Pore-C),
  0.008 (Mitotic CiFi).
- Mantel ρ: 0.457 (CiFi), 0.611 (Pore-C), 0.340 (Mitotic CiFi). All p < 1e-300.
- ARI 0.27–0.28 — comparable to LCL.
- Mitotic vs async: global W/B 128.6× vs 41.8× (3× stronger), but Mantel attenuated ~1.3× and per-PHR
  ρ attenuated ~1.4× — interphase signal partially "smoothed" by mitotic condensation (caveat:
  synchronization efficiency unknown).

### 1.12 Cross-platform reproducibility (RPE-1 Async, same cells)
PacBio CiFi vs ONT Pore-C: W/B 41.8× vs 31.9×, Mantel ρ 0.457 vs 0.611, PHR-pair ρ 0.538 vs 0.681,
ARI 0.174 vs 0.216 — Pore-C slightly stronger, consistent with multi-way contact amplification.

### 1.13 Per-individual cross-arm affinity
54 sample × community pairs: Spearman ρ = **−0.31 (p = 0.024)** between mean discordance and W/B
ratio. Significant communities have lower mean discordance (0.087 vs 0.139, Wilcoxon p = 0.013) —
intact duplicated sequence drives the 3D contact.

### 1.14 Per-community reproducibility (10 datasets)
- **C1 (D4Z4, chr4q+chr10q)**: 7.9–45.5× across all 8 datasets — most reproducible.
- **C2 (chr10p+chr18p)**: 8.6–35.3× across all 8 datasets.
- **C7 (acrocentric p-arms)**: 5/8 datasets — nucleolar.
- **C14/C15 (PARs)**: 71–165× in male diploid; zero / absent elsewhere — sex-controlled.
- **C13 (chr4p)**: only significant in deep datasets — depth-limited.

### 1.15 PBMC Dip-C cell-type split (negative)
18 PBMC cells (15 lymphocytes + 3 monocytes/neutrophils, hg19): W/B = 0.983 (1.7% closer),
Wilcoxon p = 0.305, Fisher combined p = 0.217. Not significant — attributed to hg19 projection
noise + mixed cell types + small N. GM12878 T2T (6.9% closer, p = 2.4e-5) remains the primary
Dip-C finding.

### 1.16 Lamina cross-reference
C1 (D4Z4) radial = 0.732 (peripheral, consistent with lamin A/C tethering — Masny 2004,
Ottaviani 2009). C14 (PAR2) most peripheral (0.840). C10 (chr17p) most interior (0.474).
Lamina contributes but is not the primary correlate.

---

## 2. Existing figures referenced (paths)

The chapter does not embed figure files inline, but it points to results directories that already
contain produced PDFs. Concrete paths verified to exist on `/moosefs/guarracino/HPRCv2/PHR_III/`:

### 2.1 Per-sample figure suites (exist for all 8 datasets at 50 kb)
Pattern: `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/<sample>_*`.
- `*_hic_bootstrap_distributions.pdf` — 10 000-permutation null vs observed within-community contact.
- `*_hic_community_heatmap.pdf` — arm × arm contact heatmap, blocked by community.
- `*_hic_mds_comparison.pdf` — MDS embedding of arms colored by sequence community vs Hi-C community.
- `*_hic_similarity_contact_scatter.pdf` — Mantel scatter (Jaccard similarity vs Hi-C contact).
- `*_phr_pair_scatter.pdf` — per-arm-pair Spearman scatter.

Verified samples: `chm13`, `hg002`, `hg002_cifi`, `hg002_porec`, `hg00658`, `hg02148`, `hg02559`,
`na19036` (40 PDFs total at 50 kb).

### 2.2 No-acrocentric exclusion PDFs
`/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/no_acrocentric/50000bp/<sample>_hic_*.pdf` —
matched suite for the acro-exclusion control.

### 2.3 Community-free flanking scatters
`/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_free/<sample>_flanking_<res>bp_seqlevel_scatter.pdf` (CHM13, HG002, HG002 cifi/porec, HG00658, HG02148, HG02559, NA19036; 5–100 kb).

### 2.4 RPE-1 figure suites
`/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/RPE1/{5000,10000,20000,50000,100000}bp/`
and `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/flanking/RPE1/{res}bp/` for Async CiFi, Async
Pore-C, Mitotic CiFi.

### 2.5 Pedigree-based figures (sibling section)
`/moosefs/erikg/phrs/.wg-worktrees/agent-636/end-to-end-report/pedigree-plots/` — repository-tracked,
but cross-referenced in the report rather than directly in chapter 05.

**No paper-style composite figure (multi-panel "main figure 4 / 5 / 6") yet exists** — every
asset above is a per-sample scratch PDF.

---

## 3. Existing CSVs / data files (paths)

Verified output trees on `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/`:

### 3.1 Community-based (8 samples × 5 resolutions)
`community_based/{5000,10000,20000,50000,100000}bp/<sample>_*`
- `*_community_bootstrap_tests.tsv` — per-community W/B + bootstrap p (raw + BH q).
- `*_global_test.tsv` — Mann–Whitney global p, Mantel ρ + p.
- `*_seq_community_bootstrap_tests.tsv` / `*_seq_global_test.tsv` — sequence-level (50-community)
  variant.
- `*_hic.arm-leiden.communities.tsv` — Leiden partition from Hi-C contact matrix.
- `*_hic.arm-leiden.leiden_scan.tsv` — resolution-parameter scan.
- `*_hic.dist_matrix.tsv` — arm-level Hi-C distance matrix.
- `*_contact_matrix.tsv` — raw arm × arm contact matrix.
- `*_oe_matrix.tsv` — O/E-normalized arm × arm matrix.
- `*_seq_vs_hic_ari.tsv` — ARI vs sequence partition.
- `*_phr_pair_correlation.tsv` — per-arm-pair Spearman.
- `*_subtelomeric_regions.bed` — PHR coords used.
- `*_singleton_C{8,9,10,13}_neighbors.tsv` — singleton-community contextual neighbors (CHM13).

### 3.2 Flanking (7 samples × 5 resolutions)
`flanking/{res}bp/` — same file pattern as 3.1 but on 100 kb centromere-ward windows.

### 3.3 Acrocentric-exclusion control (7 samples × 5 resolutions)
`no_acrocentric/{res}bp/` — same file pattern; B/W and Mantel after excluding chr13–22 + chr X/Y.

### 3.4 Community-free sequence-level
`community_free/<sample>_flanking_{res}bp_seqlevel.tsv` (full pair table) and `*_seqlevel_summary.tsv`
(Spearman ρ + p). Also `human_<sample>_flanking_{res}bp_seqlevel*` variants.

### 3.5 RPE-1 (3 datasets × 5 resolutions)
`community_based/RPE1/{res}bp/` and `flanking/RPE1/{res}bp/` — RPE-1 self-discovered + HPRC
community labels.

### 3.6 PHR sequence inventory (input)
`all-vs-all.1Mb.p95.id95.len.tsv` — 15 666 sequences, per-haplotype FASTAs, 0 truncated. Used for
PHR coordinate extraction.

### 3.7 PBMC Dip-C
`/moosefs/guarracino/HPRCv2/dipc_t2t/pbmc_hg19/` — 18 PBMC cells, hg19 coords, PHR boundaries
projected via impg.

---

## 4. Methods used

1. **Multi-mapping policy.** All three technologies disable MAPQ and multi-mapper filters
   (`MIN_MAPQ=0`, `RM_MULTI=0`); each multi-mapped read keeps **exactly one randomly chosen
   primary alignment** (Bowtie2 default for Hi-C; minimap2 default for Pore-C / CiFi). Adds
   symmetric noise; aggregate community-level enrichments hold; flanking unique-sequence is the
   clean-mapping control.
2. **Per-haplotype analysis.** Maternal / paternal arms kept separate (~75 arms for diploid
   samples, 92 for RPE-1, 38 for haploid CHM13).
3. **Five mcool resolutions:** 5 / 10 / 20 / 50 / 100 kb.
4. **W/B + bootstrap permutation test** — per community; 10 000 random label permutations;
   Benjamini–Hochberg FDR control at q < 0.05.
5. **Global Mann–Whitney U** on pooled within vs between contacts.
6. **Mantel test** — Spearman rank ρ on upper-triangle of arm × arm Jaccard similarity vs Hi-C
   contact matrix; 10 000 row/column-label permutations.
7. **Independent Leiden community detection** on the **O/E-normalized inter-chromosomal contact
   matrix** at 50 kb (per-haplotype), then **Adjusted Rand Index (ARI)** vs the 15 sequence-based
   communities.
8. **Per-arm-pair Spearman** ρ between mean Jaccard (averaged across all underlying sequence pairs)
   and Hi-C contact.
9. **Community-free per-sequence-pair Spearman** — direct cooler queries at PHR sequence
   coordinates; balanced contact sum normalized by region size (bp).
10. **O/E normalization:** E_ij = (row_sum_i × col_sum_j) / total_inter; observed/expected ratio
    used to control for chromosome size and marginal mappability.
11. **Acrocentric / sex / strong-community exclusion controls** (re-running Mantel + community-free
    with five exclusion sets: no acro p, no sex, no acro p+sex, no acro pq+sex, no strong [C1, C7,
    C14, C15]).
12. **Cell-type generalization (RPE-1)** — same 5-step workflow as human samples, on a near-diploid
    non-transformed cell line; cross-platform (PacBio CiFi vs ONT Pore-C); cell-cycle modulation
    (async vs mitotic).
13. **Sequence-level (50-community) variant** — same tests using the finer 50-community partition.
14. **Lamina / radial cross-reference** — Dip-C radial position joined to per-community signal.

Scripts (all under `/moosefs/guarracino/HPRCv2/scripts/`):
- `community/analyze_hic_communities.py` — community-based B/W, Mantel, ARI, per-PHR-pair, flanking.
- `community/sequence_hic_correlation.py` — community-free per-sequence-pair Jaccard × Hi-C.
- `community/plot_seqlevel_overlay.py` — sequence-level overlay scatters.
- `community/parm_qarm_3d_enrichment.py` — p-arm vs q-arm comparison.
- `similarity/export_arm_dist_matrix.R` — arm-level Jaccard distance from RDS.
- `hic/analyzer.py` — NOR clustering via MDS.
- `verify_rpe1_results.py` — global summary extraction.

---

## 5. Gaps

Items the chapter does **not** answer but the manuscript / talk likely needs:

1. **No multi-panel composite figure.** Every existing PDF is a per-sample scratch artifact.
   A main figure (e.g. "3D validation of subtelomeric communities") combining a contact heatmap, a
   Mantel scatter, an ARI bar, and a per-community W/B forest plot does not yet exist as a single
   asset.
2. **No headline summary table.** No single TSV/CSV that joins B/W, Mantel ρ, ARI, per-PHR-pair ρ,
   community-free ρ across all 8 + 3 datasets in one file. Currently scattered across
   `*_global_test.tsv` per sample.
3. **HG02148 Mantel marginal (p=0.085)** at the full-arm level is not explicitly addressed in the
   results table — only the exclusion control rescues it (0.152 → 0.720). The reason is plausibly
   chromosome fragmentation (chr1, 4, 6 as random contigs) but the report does not quantitatively
   tie HG02148's NaN flanking ρ to its marginal Mantel.
4. **Mitotic synchronization efficiency unknown** — flagged as caveat for the cell-cycle
   modulation result but not bounded.
5. **Pore-C / CiFi multi-way contact resolution** is suggestive (Pore-C consistently
   stronger than Hi-C on RPE-1) but the chapter never separates 2-way from 3+-way contacts
   numerically.
6. **PBMC Dip-C result is reported as null** but the lymphocyte-only subset (n = 15) is not
   tested in isolation in the chapter — only the all-18-cells combined number is cited.
7. **No effect size for the lamina / radial correlation.** Per-community radial values are
   given individually but no global Spearman of (radial × community W/B) is reported.
8. **No FDR correction across the 8 cross-confound exclusions** (Mantel and community-free are
   reported per-cell rather than as a joint inference).
9. **No comparison to randomized null Hi-C maps** at the per-PHR-pair scale (only the bootstrap of
   community labels is reported, not contact-matrix randomization).
10. **Chromosome-size confound for O/E** is addressed at the arm-pair level but the contact
    matrix could also be **distance-decay-normalized** for cis effects bleeding into trans —
    not done here (intentional, since the focus is inter-chromosomal).
11. **Multi-resolution Mantel "no strong" rho** is given for the 6 Hi-C samples — Pore-C / CiFi
    not in the multi-resolution exclusion table.
12. **No effect size relative to non-PHR inter-chromosomal contact** as a formal control beyond
    flanking — i.e. random matched-length unique-sequence regions away from the subtelomere are
    not used as a third negative control.

---

## 6. Suggested figures (main + SI) with captions — produced vs to-do

**Convention:** "produced" = the underlying data and per-sample PDFs already exist; "to-do" = the
multi-panel composite still needs to be assembled.

### 6.1 Main figure candidates

**Main Fig 5A — Hi-C concordance overview (TO-DO composite; data PRODUCED).**
*Caption.* Subtelomeric sequence communities have a 3D counterpart. (a) HG002 Hi-C inter-chromosomal
contact matrix at 50 kb, arms ordered by sequence community (15 communities) — diagonal blocks
visibly enriched. (b) Mantel test scatter: arm-level Jaccard similarity (x) vs Hi-C contact (y),
ρ = 0.66, p < 1e-4 (10 000 permutations). (c) ARI vs sequence partition for all 8 datasets — bar
plot, all > 0. (d) Per-community W/B ratio, log-scale, faceted by sample. Source PDFs:
`hg002_hic_community_heatmap.pdf`, `hg002_hic_similarity_contact_scatter.pdf`,
`<sample>_seq_vs_hic_ari.tsv` (×8), per-community `*_community_bootstrap_tests.tsv`.

**Main Fig 5B — Multi-resolution + acrocentric robustness (TO-DO composite; data PRODUCED).**
*Caption.* The signal is resolution- and confound-invariant. (a) B/W ratio across 5 mcool
resolutions for all 8 datasets — all far below 1.0. (b) Mantel ρ before vs after "no acro pq + sex"
exclusion: HG02148 0.152→0.720; CHM13 0.656→0.850. (c) Community-free per-sequence-pair Spearman ρ
across exclusions. (d) O/E within vs between violin plot, 8.6–34.4× separation.

### 6.2 Supplementary figure candidates

**SI Fig S5.1 — Bootstrap null distributions, 8 datasets (50 kb).** Per-sample 10 000-permutation
null vs observed within-community mean contact. PRODUCED:
`<sample>_hic_bootstrap_distributions.pdf` (×8).

**SI Fig S5.2 — Per-sample MDS comparison.** Arms colored by sequence community vs by Hi-C
community. PRODUCED: `<sample>_hic_mds_comparison.pdf` (×8).

**SI Fig S5.3 — Per-arm-pair Spearman scatter, 8 datasets.** PRODUCED: `<sample>_phr_pair_scatter.pdf` (×8).

**SI Fig S5.4 — Flanking-region community-free overlays.** Spearman ρ at 100 kb centromere-ward
unique sequence. PRODUCED:
`community_free/<sample>_flanking_<res>bp_seqlevel_scatter.pdf` (multi-resolution).

**SI Fig S5.5 — Acrocentric / sex / strong-community exclusion grid.** Mantel ρ matrix [sample ×
exclusion]; community-free ρ matrix. TO-DO composite from existing TSVs.

**SI Fig S5.6 — RPE-1 cell-type validation panel.** B/W across resolutions for Async CiFi, Async
Pore-C, Mitotic CiFi (HPRC + self-discovered communities); Mantel scatters; ARI. PRODUCED per-dataset
PDFs in `community_based/RPE1/{res}bp/`.

**SI Fig S5.7 — Async vs mitotic cell-cycle modulation.** RPE-1 paired comparison: W/B 41.8× vs
128.6×, Mantel ρ 0.457 vs 0.340, PHR-pair ρ 0.538 vs 0.389. TO-DO composite.

**SI Fig S5.8 — Per-community reproducibility heatmap.** 15 communities × 10 datasets, cell color =
log enrichment ratio, asterisks for q < 0.05 / 0.001. PRODUCED data: per-sample
`*_community_bootstrap_tests.tsv` × 8 + RPE-1 ×3.

**SI Fig S5.9 — Lamina / radial position cross-reference.** Per-community Dip-C radial (x) vs Hi-C
W/B (y) — TO-DO; needs joining lamina/radial table to community W/B.

**SI Fig S5.10 — Multi-mapping symmetry diagram.** Schematic showing that randomly assigning a
multi-mapper to chr4q vs chr10q adds **symmetric** noise, so aggregate community contacts hold but
individual pair-level repetitive contacts do not. TO-DO; conceptual diagram.

### 6.3 Tables (suggested)

- **Main Table 5.** Joined headline statistics across 8 + 3 datasets: B/W, Mantel ρ + p, ARI,
  per-arm-pair ρ, community-free ρ — currently scattered, needs single summary CSV (TO-DO).
- **SI Table 5.1.** Multi-resolution B/W (PHR + flanking + no-acro), 5 mcool resolutions × 8
  samples — PRODUCED data, TO-DO compile.
- **SI Table 5.2.** Per-community enrichment table (15 × 10) with q-values — PRODUCED data, TO-DO
  compile.
- **SI Table 5.3.** Cross-confound Mantel ρ (full / no acro p / no sex / no acro p+sex / no acro
  pq+sex / no strong) — PRODUCED, mostly already laid out as text tables (lines 329–351).
- **SI Table 5.4.** RPE-1 vs LCL comparison (Mantel, ARI, PHR-pair, B/W) — PRODUCED, needs
  consolidation.

---

## 7. Talk slide takeaways (15-min talk)

Framed for a Nature-style 15-minute talk on subtelomeric sequence sharing across HPRCv2:

1. **"Sequence communities are 3D-real."** Arms placed in the same Leiden subtelomeric community by
   pangenome graph similarity also touch each other in the nucleus — across 8 datasets, 3
   technologies (Hi-C, Pore-C, CiFi), and 5 resolutions (5–100 kb).

2. **"Three orthogonal tests, same conclusion."** Within/between contact ratio (B/W = 0.027–0.074,
   all p < 0.01), Mantel arm-similarity × contact (ρ up to 0.66), and **independent Leiden on the
   Hi-C O/E matrix recovers the sequence partition (ARI = 0.06–0.54, all > 0)**. Strongest:
   HG002 Hi-C — Mantel ρ = 0.66, B/W = 0.027 (p = 4e-66).

3. **"Not driven by acrocentrics or sex."** Excluding chr13–22 and chrX/Y *strengthens* the signal
   in 7/8 datasets (HG002 Mantel 0.66 → 0.79; HG02148 0.15 → 0.72). This rules out nucleolar
   association and PAR sharing as the sole drivers.

4. **"Not a chromosome-size artifact."** O/E normalization preserves an 8.6–34.4× within/between
   contact gap. The signal is community-specific, not size-confounded.

5. **"Even unique-sequence flanks show it."** 100 kb of unique sequence centromere-ward of the PHR
   boundary still shows significant community structure (CHM13 Mantel 0.52, HG002 0.52). The
   contact signal is a **chromosomal-domain** property, not a multi-mapping artifact.

6. **"Sub-arm-pair resolution."** Community-free Spearman on 627–2 830 individual sequence pairs:
   ρ = 0.66–0.83 in Hi-C, all p < 3.7e-24. The proximity signal is graded — more shared
   subtelomeric sequence ⇒ more 3D contact, no community labels needed.

7. **"Cell-type general."** RPE-1 (non-transformed retinal) shows the same enrichment as LCLs
   (Mantel ρ 0.34–0.61, ARI 0.27–0.28). **Cross-platform** (PacBio CiFi vs ONT Pore-C, same cells)
   gives matched results — Pore-C slightly stronger via multi-way contacts.

8. **"Cell-cycle modulated."** Mitotic-arrested CiFi: global W/B 3× stronger but per-arm-pair ρ ~1.4×
   weaker — condensation amplifies the bulk segregation but smooths fine-grained arm-pair
   correspondence.

9. **"Mechanism hint."** Per-individual: communities with **more concordant** (intact-duplicated)
   sequence sharing have stronger 3D contact (ρ = −0.31, p = 0.024 between discordance and W/B).
   Suggests subtelomeric exchange both **leaves** and **requires** physical proximity.

10. **"D4Z4 (C1) and PARs (C14/C15) are textbook controls."** D4Z4 community is peripheral
    (radial 0.732 — lamin A/C), reproduced in all 8 datasets. PARs reach 71–165× enrichment in male
    diploids and zero / absent in female / haploid — exactly as expected.

**One-slide summary line:**
*"Sequence-defined subtelomeric communities are physical: 8 samples × 3 technologies × 5 resolutions
all show the same arms touching, with effect sizes that survive every confound we threw at them
(acrocentrics, sex, chromosome size, mappability, cell type, cell cycle)."*

---

**Inputs cited.**
- `end-to-end-report/report/05_hic_validation.md` (lines 1–689)
- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/{5,10,20,50,100}000bp/`
- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/no_acrocentric/{5,10,20,50,100}000bp/`
- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_free/`
- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/RPE1/{res}bp/`
- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/flanking/RPE1/{res}bp/`
- `/moosefs/guarracino/HPRCv2/dipc_t2t/pbmc_hg19/`
- `/moosefs/guarracino/HPRCv2/scripts/community/{analyze_hic_communities.py, sequence_hic_correlation.py, plot_seqlevel_overlay.py, parm_qarm_3d_enrichment.py}`
- `/moosefs/guarracino/HPRCv2/scripts/similarity/export_arm_dist_matrix.R`
- `/moosefs/guarracino/HPRCv2/scripts/hic/analyzer.py`
- `/moosefs/guarracino/HPRCv2/scripts/verify_rpe1_results.py`
