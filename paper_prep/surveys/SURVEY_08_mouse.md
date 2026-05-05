---
title: "Survey 08 — Mouse cross-species T2T subtelomeric validation"
source: end-to-end-report/report/08_mouse.md
scope: Mouse T2T (B6 + CAST) subtelomeric PHR detection, communities, multi-window optimization, Hi-C across 4 meiotic stages, flanking Hi-C, sequence-level community-free
audience: Nature manuscript and 15-min talk
---

# Survey 08 — Mouse cross-species section

This survey extracts and structures the content of `end-to-end-report/report/08_mouse.md` for the Nature manuscript and the companion 15-minute talk. The source section asks: does the subtelomeric **sequence homology → 3D proximity** link generalise beyond human? It applies the human pipeline (contig classification → community detection → Hi-C validation) to two mouse T2T assemblies (C57BL/6J = B6 and CAST/EiJ from Francis et al. 2025), exploits mouse-specific telocentric architecture (one subtelomere per chromosome), runs a 1Mb→2Mb→4Mb window-size sweep up to 33Mb, and adds a flanking-Hi-C control. It is the cross-species generalisation arm of the paper.

---

## 1. Key findings with metrics

### 1.1 Pipeline executes end-to-end on mouse T2T (B6 + CAST)

- Two assemblies: B6 = `GCA_964188535`, CAST = `GCA_964188545` (Francis et al. 2025).
- 78 total subtelomeric flanks (B6 19 autosomes + CAST 19 autosomes + chrX); flank length 474–498 kb after telomere trimming (telomeric tracts: B6 1,695–26,101 bp; CAST 3,329–21,637 bp).
- 49,911 inter-chromosomal alignments (`wfmash -p 95`, do_not_overfilter branch cc60cd8 — k-mer-frequency filtering no longer needs the manual `-F 0.1` override; `--min-count 1`).

### 1.2 Subtelomeric end is uniquely subtelomeric (centromere-distal only)

- 39/78 flanks carry inter-chromosomal signal.
- **All 39 `_parm` flanks (centromere-distal / subtelomeric)** have signal.
- **All 39 `_qarm` flanks (centromere-proximal)** have **zero** signal — satellite-dominated; confirms centromere coordinates are unnecessary for mouse and matches the human "subtelomeric ends only" finding.

### 1.3 Mouse subtelomeric sharing is broad and uniform — far beyond the TLC repeat

- Most subtelomeric flanks share inter-chromosomal sequence over the **entire 500 kb window** (440–495 kb regions match all 20 chromosomes).
- The TLC repeat itself is only 6–12 kb in B6 (Francis et al. 2025); the shared region therefore extends ~50× beyond the TLC motif into flanking repeat-rich (L1-LINE / LTR) sequence.
- Notable B6 exceptions (consistent with Francis et al. structural exceptions):
  - **B6 chr11/chr18** form a private group with chr4 (chr11 80 kb, chr18 30 kb signal). B6 chr11 lacks the standard L1-LINE → TLC → minor-satellite end motif.
  - **B6 chr7** forms a small private group with chrX/chr11 (10 kb).
  - **B6 chr4** has only 15 kb signal (vs typical ~490 kb), missing chr7 and chrX.

### 1.4 Pairwise Jaccard (1Mb p-arm flanks; PGGB `-p 95 -n 2`, `odgi similarity --all`)

| Pair | Jaccard |
|---|---:|
| CAST chr1 ↔ CAST chr2 | **0.987** |
| B6 chr3 ↔ B6 chr15 | 0.944 |
| B6 chr10 ↔ B6 chr12 | 0.934 |
| CAST chr7 ↔ CAST chr16 | 0.914 |
| Most B6 cross-chr p-arm pairs | 0.88–0.94 |
| Most CAST cross-chr p-arm pairs | 0.70–0.91 |
| B6 chr11 ↔ B6 chr18 (outlier) | **0.292** |

**Take-home:** subtelomeric content is broadly shared (J = 0.7–0.99) across most p-arm pairs — the mouse subtelomeric landscape is much more uniform than human; this directly explains why Leiden returns only 2 communities.

### 1.5 Two-community Leiden structure (arm level, B6+CAST collapsed)

| Window | Communities | Notes |
|---|---:|---|
| 1 Mb | 2 | C1 = 16 arms (mostly p-arms); C2 = 11 arms (7 q-arms + chr4_p, chr7_p, chr11_p, chr18_p) |
| 2 Mb | 2 | identical to 1 Mb |
| 4 Mb | 2 | chr11_p moves C2 → C1 (17 vs 10 arms); 2-community structure preserved |

**Contrast with human:** 15 multi-arm communities across 41 arms with overlapping duplicon structure — mouse is fundamentally simpler.

### 1.6 Mouse private pairs are NOT human subtelomeres (mm39 → hg38 syntenic net)

| Mouse chr (distal end) | Human syntenic region | Human subtelomeric community? |
|---|---|---|
| CAST chr1 (J=0.987 pair) | chr8: 51.8–55.6 Mb | interior — no |
| CAST chr2 (J=0.987 pair) | chr10: 5.9–15.4 Mb | interior — no |
| B6 chr11 (J=0.292 pair) | chr22: 28.8–31.6 Mb | interior — no |
| B6 chr18 (J=0.292 pair) | chr10: 35.0–35.2 Mb | interior — no |

Mouse subtelomeric sharing is driven by **repeat architecture (TLC, L1-LINE)**, not by syntenic conservation of a duplicon system — exactly as expected when mouse is telocentric and human is metacentric.

### 1.7 Hi-C validation across 4 meiotic stages (Zuo et al. 2021)

**Community-free, sequence-level (1Mb window, 50 kb resolution):**

| Stage | Seq pairs | ρ (PHR) | p (PHR) | ρ (flanking) | p (flanking) |
|---|---:|---:|---:|---:|---:|
| Leptotene | 1,088 | 0.372 | 4.0e-37 | 0.740 | 3.4e-74 |
| Zygotene  | 1,135 | **0.425** | 4.2e-51 | 0.604 | 1.2e-68 |
| Pachytene | 1,135 | **0.428** | 9.6e-52 | 0.766 | 7.2e-88 |
| Diplotene | 999   | **0.416** | 4.9e-43 | 0.715 | 1.1e-62 |

**Community-based, per arm-pair (1Mb, 50 kb resolution):**

| Stage | rho | p | Mantel ρ | Mantel p |
|---|---:|---:|---:|---:|
| Leptotene | 0.680 | 5.0e-48 | 0.687 | < 0.0001 |
| Zygotene  | **0.715** | 4.4e-55 | 0.718 | < 0.0001 |
| Pachytene | 0.677 | 1.6e-47 | 0.683 | < 0.0001 |
| Diplotene | 0.574 | 1.8e-31 | 0.577 | < 0.0001 |

- All four meiotic stages: significant positive correlation between sequence similarity and Hi-C contact, both PHR and flanking.
- Flanking ρ (0.60–0.77) > PHR ρ (0.37–0.43) — the **human flanking paradox** reproduces in mouse: unique-sequence flanks beat repeat-rich PHR centres because flanks are not multi-mapping-suppressed.
- Compared to old 500 kb pipeline (rho 0.08–0.17, n.s.), 1Mb captures the full PHR extent.

### 1.8 Window-size optimization (1Mb–33Mb)

| Window | n PHR | Median | Mean | Max | Saturated (≥ 90% window) |
|---|---:|---:|---:|---:|---|
| 1Mb  | 49 | 0.98 Mb | 0.69 Mb | 1.0 Mb  | 30/49 (61 %) |
| 2Mb  | 49 | 1.85 Mb | 1.30 Mb | 2.0 Mb  | 25/49 (51 %) |
| 4Mb  | 49 | 2.53 Mb | 2.30 Mb | 4.0 Mb  | 19/49 (39 %) |
| 10Mb | 49 | 3.42 Mb | 3.92 Mb | 10.0 Mb |  8/49 (16 %) |
| 15Mb | 50 | 3.44 Mb | 4.32 Mb | 15.0 Mb |  4/50 ( 8 %) |
| 33Mb | 50 | 2.11 Mb | 3.37 Mb | 18.9 Mb |  0/50 ( 0 %) |

- Largest mouse PHR: **CAST chr1_p = 18.9 Mb**.
- Median stabilises between 10–33 Mb at ≈ 2 Mb — most PHRs are 1–4 Mb.
- 2 Mb is the practical sweet spot (resolves 6 arms truncated at 1 Mb at modest extra cost).

**Identity-threshold scan at 33 Mb:**

| Identity | n PHR | Median | Mean | Max |
|---|---:|---:|---:|---:|
| ≥95 % | 50 | 2.11 Mb | 3.37 Mb | 18.9 Mb |
| ≥96 % | 48 | 1.69 Mb | 2.64 Mb | 18.9 Mb |
| ≥97 % | 42 | 1.71 Mb | 2.60 Mb | 18.9 Mb |
| ≥98 % | 38 | 0.60 Mb | 1.84 Mb | 18.7 Mb |

- Raising threshold 95→98 % cuts PHR count by 24 % and median size to less than a third.
- Outer mouse subtelomeres are **moderately diverged (95–98 %) repeats**; a high-identity (≥ 98 %) core extends 0.5–2 Mb from the telomere.

### 1.9 Window-scaling Hi-C signal (community-based, 50 kb resolution)

| Window | Stage | B/W ratio | p | Mantel ρ |
|---|---|---:|---:|---:|
| 1 Mb | Lep | 0.055 | 4.7e-35 | 0.687 |
| 1 Mb | Zyg | 0.093 | 4.0e-93 | 0.718 |
| 1 Mb | Pac | 0.071 | 9.3e-61 | 0.683 |
| 1 Mb | Dip | 0.046 | 2.7e-23 | 0.577 |
| 2 Mb | Lep | 0.042 | 7.2e-38 | 0.700 |
| 2 Mb | Zyg | 0.072 | 1.2e-103 | 0.685 |
| 2 Mb | Pac | 0.056 | 2.4e-60 | 0.680 |
| 2 Mb | Dip | 0.038 | 3.3e-26 | 0.609 |
| 4 Mb | Lep | 0.026 | 6.4e-100 | 0.727 |
| 4 Mb | Zyg | 0.038 | 1.1e-150 | 0.693 |
| 4 Mb | Pac | 0.030 | 1.3e-110 | 0.711 |
| 4 Mb | Dip | 0.034 | 1.8e-50 | 0.650 |

- B/W ratio (lower = stronger enrichment) **decreases monotonically** with window size: 1Mb ≈ 0.05–0.09 → 4Mb ≈ 0.03–0.04. Bigger window → more sharing signal captured.
- Mantel ρ is high and stable (0.58–0.73) across all windows.

### 1.10 Mouse flanking Hi-C (1 Mb window, 100 kb centromere-ward of PHR)

**B/W ratios — flanking is even more enriched than PHR centres:**

| Stage | 5kb | 10kb | 20kb | 50kb | 100kb |
|---|---:|---:|---:|---:|---:|
| Lep | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| Zyg | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| Pac | 0.003 | 0.002 | 0.002 | 0.003 | 0.003 |
| Dip | 0.001 | 0.001 | 0.001 | 0.001 | 0.002 |

**Mantel ρ:**

| Stage | 5kb | 10kb | 20kb | 50kb | 100kb |
|---|---:|---:|---:|---:|---:|
| Lep | 0.656 | 0.662 | 0.597 | 0.588 | 0.589 |
| Zyg | 0.739 | 0.750 | 0.722 | 0.643 | 0.566 |
| Pac | 0.712 | 0.683 | 0.662 | 0.612 | 0.622 |
| Dip | 0.644 | 0.622 | 0.617 | 0.495 | **0.100 (n.s.)** |

- Significant Mantel at all stages × resolutions except diplotene at 100 kb.
- Flanking B/W (0.001–0.003) much lower than PHR B/W (0.03–0.12) — extreme enrichment, consistent with unique-sequence freeing from multi-mapping suppression.
- **Caveat:** mouse "flanking" still contains 138–297 nonzero-Jaccard pairs (vs human flanks which are ~unique). The flanking negative-control concept applies *less cleanly* to mouse because TLC + L1-LINE extends well past the PHR boundary.

### 1.11 Multi-resolution community-based B/W ratios

**1 Mb window:**

| Stage | 5kb | 10kb | 20kb | 50kb | 100kb |
|---|---:|---:|---:|---:|---:|
| Lep | 0.073 | 0.055 | 0.057 | 0.055 | 0.029 |
| Zyg | 0.122 | 0.112 | 0.105 | 0.093 | 0.108 |
| Pac | 0.119 | 0.098 | 0.091 | 0.071 | 0.076 |
| Dip | 0.055 | 0.050 | 0.040 | 0.046 | 0.061 |

**4 Mb window:**

| Stage | 5kb | 10kb | 20kb | 50kb | 100kb |
|---|---:|---:|---:|---:|---:|
| Lep | 0.030 | 0.022 | 0.024 | 0.026 | 0.020 |
| Zyg | 0.057 | 0.049 | 0.044 | 0.038 | 0.042 |
| Pac | 0.052 | 0.037 | 0.041 | 0.030 | 0.036 |
| Dip | 0.024 | 0.019 | 0.023 | 0.034 | 0.043 |

All B/W ≪ 1 across the entire grid — strong community-structured 3D enrichment in every meiotic stage at every resolution.

---

## 2. Existing figures (paths)

All paths absolute, on the shared `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/` tree.

### Community-based Hi-C validation (per stage × resolution × window)
Pattern: `community_analysis_{1,2,4}Mb/{5000,10000,20000,50000,100000}bp/zuo2021_{lep,zyg,pac,dip}tene_*.pdf`

For each meiotic stage there are four PDFs:
- `*_hic_community_heatmap.pdf` — Hi-C contact matrix grouped by community
- `*_hic_mds_comparison.pdf` — MDS embedding comparing sequence-distance vs Hi-C-distance arrangement
- `*_hic_similarity_contact_scatter.pdf` — community-based per-pair similarity vs contact scatter
- `*_phr_pair_scatter.pdf` — community-free per-PHR-pair scatter (the paper-headline 16-panel design analog)

Total: 4 stages × 5 resolutions × 3 windows × 4 figure types = **240 PDFs** at `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_{1,2,4}Mb/{res}bp/`.

Representative anchors (1 Mb, 50 kb resolution):
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp/zuo2021_zygotene_phr_pair_scatter.pdf`
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp/zuo2021_zygotene_hic_community_heatmap.pdf`
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp/zuo2021_zygotene_hic_mds_comparison.pdf`
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp/zuo2021_zygotene_hic_similarity_contact_scatter.pdf`

### Sequence-level community-free correlation
Pattern: `seqlevel_correlation_{1,2,4}Mb/mouse_{stage}_{phr,flanking}_{5000,10000,20000,50000,100000}bp_seqlevel_scatter.pdf`

Examples:
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/seqlevel_correlation_1Mb/mouse_pachytene_phr_50000bp_seqlevel_scatter.pdf`
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/seqlevel_correlation_1Mb/mouse_pachytene_flanking_50000bp_seqlevel_scatter.pdf`

### Flanking Hi-C
Pattern: `flanking_hic_{1,2,4}Mb/{res}bp/zuo2021_{stage}_*.pdf` (same 4 figure types as community_analysis):
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/flanking_hic_1Mb/50000bp/zuo2021_zygotene_phr_pair_scatter.pdf`
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/flanking_hic_1Mb/50000bp/zuo2021_zygotene_hic_community_heatmap.pdf`

### Pangenome graph (similarity)
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/flanking_1Mb/pggb/mouse_flank1M_combined.paf.11fba48.e3aa42b.smooth.final.gfa` (graph)
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/flanking_1Mb/pggb/mouse_flank1M_combined.paf.11fba48.e3aa42b.smooth.final.og` (odgi)

> **Note:** the PGGB build for mouse subtelomeres lives under `flanking_1Mb/pggb/`. No similarity-matrix PDF is rendered in the source section — only the per-pair Jaccard tables (see §3).

---

## 3. Existing CSVs (paths)

### Community assignments and Leiden scans (per window)
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/similarity_{1,2,4}Mb/`
- `mouse.communities.tsv` — final per-sequence community labels.
- `mouse.arm-leiden.communities.tsv` — arm-collapsed assignments (the ones quoted in the source).
- `mouse.leiden_scan.tsv`, `mouse.arm-leiden.leiden_scan.tsv` — resolution scan 0.1–3.0.
- `mouse.dist_matrix.tsv` — pairwise Jaccard distance matrix.
- `similarity_1Mb/mouse.communities.tsv.bak.3comm` — pre-fix backup (3-community variant; superseded).

### Pairwise similarity (PGGB / odgi)
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/flanking_1Mb/pggb/mouse_flank1M_combined.paf.11fba48.e3aa42b.smooth.final.similarity.tsv.gz` — full pairwise.
- `…similarity.tsv.gz.arm_pair_jaccard.tsv` — arm-pair-collapsed Jaccard (source of §1.4 table).

### PHR detection per window
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/subtelo_{1,2,4,10,15,33}Mb/`
- `mouse.all-vs-all.p95.id95.len.tsv` — PHR length per detected region (source of §1.8).
- `mouse.all-vs-all.p95.id95.len.per_window.tsv` (2Mb, 33Mb only).
- `B6.telo.tsv`, `CAST.telo.tsv` — telomere coordinates (source of §1.1 telomeric tract sizes).
- `mouse_subtelo_{N}Mb_trimmed.fa.gz` — telomere-trimmed subtelomere flanks (input to wfmash).

Identity-threshold scan: `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/subtelo_33Mb_p{96,97,98}/` (same files as 33 Mb base).

### Community-free per-pair correlation tables
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/seqlevel_correlation_{1,2,4}Mb/`
- `mouse_{lep,zyg,pac,dip}tene_{phr,flanking}_{5,10,20,50,100}kb_seqlevel.tsv` — per-pair raw points.
- `mouse_{stage}_{kind}_{res}bp_seqlevel_summary.tsv` — Spearman ρ, p, n (source of §1.7 §1.10).

### Community-based Hi-C tests
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_{1,2,4}Mb/{res}bp/`
- `zuo2021_{stage}_community_bootstrap_tests.tsv` — B/W ratios + bootstrap p (source of §1.9 §1.11).
- `zuo2021_{stage}_global_test.tsv` — Mantel ρ + permutation p.
- `zuo2021_{stage}_contact_matrix.tsv` — community × community contact matrix.
- `zuo2021_{stage}_phr_pair_correlation.tsv` — per-pair correlations.
- `mouse_{stage}_subtelomeric_regions.bed` — region BEDs used by the analysis.

### Flanking Hi-C tables
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/flanking_hic_{1,2,4}Mb/{res}bp/` — same file set as community_analysis (source of §1.10).

### Hi-C inputs
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/HiC/zuo2021_{leptotene,zygotene,pachytene,diplotene}.mcool` — per-stage mouse meiotic Hi-C (Zuo et al. 2021).

### Other
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/annotations/` — repeat / GFF annotations (input for the open repeat-annotation question, §5).

---

## 4. Methods

### 4.1 Inputs and assemblies
- B6 (C57BL/6J): `GCA_964188535` → `/moosefs/guarracino/HPRCv2/data/mouse_T2T/B6.PanSN.fa`
- CAST (CAST/EiJ): `GCA_964188545` → `/moosefs/guarracino/HPRCv2/data/mouse_T2T/CAST.PanSN.fa`
- PanSN-prefixed (B6 = PATERNAL, CAST = MATERNAL — kept separate per haplotype throughout).
- Hi-C: 4 meiotic stages from Zuo et al. 2021 (`zuo2021_{leptotene,zygotene,pachytene,diplotene}.mcool`).

### 4.2 Subtelomeric flank extraction
- 500 kb pipeline was deleted (PHRs filled the entire window; 30/49 PHRs truncated). Replaced by a 1 Mb pipeline (0 truncated).
- Flanks: 500 kb (legacy) and {1, 2, 4, 10, 15, 33} Mb (current sweep) extracted from both chromosome ends.
- Telomeres trimmed in-place by `trim_mouse_flanks.py` (preserves PanSN naming). B6 telomere tracts 1,695–26,101 bp; CAST 3,329–21,637 bp.
- Mouse-specific advantage: telocentric chromosomes have one subtelomere per chromosome; the centromere-proximal end (`_qarm`) yields zero inter-chromosomal signal in 39/39 cases — no centromere mask required.

### 4.3 Inter-chromosomal alignment + PHR detection
- `wfmash -p 95` on the **`do_not_overfilter` branch (commit `cc60cd8`)** — k-mer-frequency filtering no longer requires the `-F 0.1` override that was needed on the human pipeline. Total alignments at 1 Mb: 49,911.
- `impg v0.4.0` for projection.
- `find-multichr-regions-incremental.py` with `--min-count 1` — same script as human, no modification.
- PHR identity sweep at 33 Mb at p {95, 96, 97, 98}.

### 4.4 Pangenome graph and similarity
- `pggb -p 95 -n 2` on 39 PHR sequences (number 39 = signal-bearing flanks).
- `odgi similarity --all` → pairwise Jaccard.
- Arm-level collapse (B6+CAST per arm) for the Leiden input — same convention as human.

### 4.5 Leiden community detection
- `detect_communities.R --organism mouse --level arm`, resolution scan 0.1–3.0.
- Input = arm-level Jaccard distance matrix from §4.4.
- Output: 27 chromosomal arms partition into 2 communities at 1, 2 and 4 Mb (chr11_p flips C2 → C1 between 1 Mb and 4 Mb, but the 2-community structure is preserved).

### 4.6 Hi-C validation (community-based)
- Driver: `analyze_hic_communities.py` — per-haplotype (B6 PATERNAL, CAST MATERNAL kept separate), arm-level communities, PHR-specific coordinates from the per-window-size PHR TSV.
- Resolutions: 5, 10, 20, 50, 100 kb.
- Stages: leptotene, zygotene, pachytene, diplotene.
- Outputs: B/W bootstrap (within-community vs between-community Hi-C contact), Mantel ρ vs sequence distance, contact matrix, MDS comparison.

### 4.7 Hi-C validation (community-free, sequence-level)
- Driver: `sequence_hic_correlation.py` (mouse mode) — analog of the human individual sequence-pair correlation.
- For each PHR-pair on different chromosomes, correlate Jaccard similarity ↔ Hi-C contact at the exact PHR coordinates; no community labels.
- Run separately for **PHR-region** (repeat-rich centre) and **flanking** (100 kb centromere-ward, see §4.8).

### 4.8 Flanking Hi-C
- Extraction: `extract_flanking_sequences.py` → 100 kb immediately centromere-ward of each PHR boundary.
- Analyzed at all 5 resolutions × 4 stages × 3 windows (1, 2, 4 Mb).
- Mouse caveat: 138–297 nonzero-Jaccard flanking pairs persist (vs human ~unique flanks); the "negative control" interpretation is therefore weaker for mouse.

### 4.9 Cross-species synteny (open question)
- mm39 → hg38 syntenic net from UCSC; used to map mouse private-pair distal ends to human positions (B6 chr11/18 and CAST chr1/2 → human chromosome interiors, never subtelomeres).

### 4.10 Scripts (canonical paths)
| Script | Purpose |
|---|---|
| `/moosefs/guarracino/HPRCv2/scripts/trim_mouse_flanks.py` | Trim telomeres from mouse flanks (preserves PanSN). |
| `/moosefs/guarracino/HPRCv2/scripts/find-multichr-regions-incremental.py` | PHR detection (same as human). |
| `/moosefs/guarracino/HPRCv2/scripts/community/detect_communities.R` | Leiden community detection (`--organism mouse --level arm`). |
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze_hic_communities.py` | Community-based Hi-C validation (per-haplotype, PHR coords). |
| `/moosefs/guarracino/HPRCv2/scripts/community/sequence_hic_correlation.py` | Community-free per-pair correlation. |
| `/moosefs/guarracino/HPRCv2/scripts/extract_flanking_sequences.py` | Flanking-region extraction. |

---

## 5. Gaps

1. **Repeat annotation not yet intersected.** RepeatMasker / GFF3 from Ensembl (`https://projects.ensembl.org/mouse_genomes/`) for both T2T assemblies are downloaded into `mouse_T2T/annotations/` but not yet intersected with PHR regions. The 440–495 kb shared block dwarfs the 6–12 kb TLC repeat — what fraction is L1-LINE vs LTR vs minor satellite vs other? This is the single biggest open question and determines whether the "subtelomeric duplicon" framing holds in mouse.
2. **CAST repeat annotation availability uncertain.** Source notes Ensembl as the source for both T2T assemblies, but CAST repeat tracks may be less complete than B6 — needs verification before §1 figures can quote per-class fractions.
3. **No 3D / CCAN / single-cell mouse data.** Source has only bulk meiotic Hi-C (Zuo 2021); no mouse Dip-C / Pore-C / scHi-C analog is shown. The cross-tech consistency claim (positive correlation across modalities) is human-only.
4. **No statistical comparison of B/W slope across windows.** Source quotes the monotone B/W decrease (1Mb → 4Mb) but does not formally test whether window scaling is significant or asymptotic — useful for the "2 Mb is the sweet spot" recommendation.
5. **No haplotype-resolved Hi-C (B6 vs CAST) breakdown.** Pipeline runs "per-haplotype" (kept separate during alignment) but the final correlations are pooled across B6 + CAST; no B6-only or CAST-only stage table is provided, so we can't say whether the signal differs between strains.
6. **F1-hybrid Hi-C not used here.** `analyze_hic_communities.py` advertises "F1 hybrid support" in its docstring, but the source uses Zuo 2021 (non-F1, B6-mapped Hi-C). An F1 (B6 × CAST) phased Hi-C dataset would test whether trans contacts remain enriched between haplotypes of *different* parental origin.
7. **Synteny analysis is qualitative.** Only 4 mouse private-pair endpoints mapped to hg38; no genome-wide enrichment test of "do mouse subtelomeric PHRs preferentially map to human subtelomeric communities?" (expected null, but worth quantifying).
8. **Flanking control caveat under-quantified.** Section says mouse flanking is not truly unique (138–297 nonzero Jaccard pairs) but doesn't tabulate the residual PHR-vs-flanking sequence-divergence distribution — a key gap because the human flanking paradox depends on flanks being clean.
9. **Diplotene 100 kb flanking outlier.** Mantel ρ collapses to 0.100 (n.s.) at diplotene 100 kb only; not explained — could be Hi-C coverage drop-off or genuine biological compaction at diplotene exit.
10. **No talk-ready figure synthesising the cross-species message.** Source has 240+ PDFs distributed by stage × resolution × window but no single-panel "human vs mouse" comparison panel.
11. **Window-saturation test is one-sided.** PHR sizes peak at CAST chr1_p = 18.9 Mb (33 Mb scan) but minimum-PHR detection threshold not stress-tested — could very small PHRs (< 100 kb) be missed at small windows?
12. **2 Mb communities table not shown explicitly.** Source says "2 Mb identical to 1 Mb" but no per-arm table for 2 Mb is rendered; users must infer from 1 Mb table.

---

## 6. Suggested figures with captions (produced vs to-do)

### Already produced (use directly or recompose)

**P-1. Mouse PHR-pair similarity-vs-contact scatter (zygotene, 1 Mb, 50 kb).**
*File:* `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp/zuo2021_zygotene_phr_pair_scatter.pdf`.
*Caption:* "Mouse community-free per-PHR-pair correlation, zygotene meiotic Hi-C at 50 kb resolution. Each point = one inter-chromosomal PHR pair; Spearman ρ = 0.425, p = 4.2 × 10⁻⁵¹ (n = 1,135). Direct mouse analog of the human zygotene panel — same sign and similar effect size."

**P-2. Mouse flanking similarity-vs-contact scatter (pachytene, 1 Mb, 50 kb).**
*File:* `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/seqlevel_correlation_1Mb/mouse_pachytene_flanking_50000bp_seqlevel_scatter.pdf`.
*Caption:* "Flanking-region (100 kb centromere-ward of PHR) similarity-vs-Hi-C-contact at pachytene, 50 kb resolution. ρ = 0.766, p = 7.2 × 10⁻⁸⁸. Flanking ρ exceeds PHR ρ — same flanking paradox as human."

**P-3. Mouse community heatmap (zygotene, 1 Mb, 50 kb).**
*File:* `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp/zuo2021_zygotene_hic_community_heatmap.pdf`.
*Caption:* "Mouse meiotic Hi-C contact matrix grouped by Leiden community (C1, C2). Within-community contact density visibly exceeds between-community — mirrored by B/W = 0.093 (zygotene, 1 Mb)."

**P-4. Mouse Hi-C MDS comparison (pachytene, 1 Mb, 50 kb).**
*File:* `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp/zuo2021_pachytene_hic_mds_comparison.pdf`.
*Caption:* "MDS embedding of mouse arms by sequence distance (left) vs Hi-C distance (right) at pachytene; concordant clustering of C1 vs C2 confirms the sequence ↔ 3D mapping at the arm level."

**P-5. Mouse flanking Hi-C scatter (zygotene, 5 kb, 1 Mb window).**
*File:* `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/flanking_hic_1Mb/5000bp/zuo2021_zygotene_phr_pair_scatter.pdf`.
*Caption:* "Mouse flanking-region community-free correlation at the highest resolution (5 kb). Mantel ρ = 0.739, B/W = 0.001 — strongest enrichment in the dataset."

### To-do (suggested new figures)

**T-1. Cross-species headline panel.**
*Caption:* "Sequence similarity ↔ Hi-C contact correlation in human (T2T-CHM13, lymphoblast / Pore-C / CiFi) vs mouse (B6 + CAST, 4 meiotic stages). Bars = Spearman ρ at 50 kb resolution. The same positive correlation holds across species, despite mouse subtelomeric architecture being telocentric and TLC-driven rather than syntenic-duplicon-driven." Built from `seqlevel_correlation_1Mb/*_summary.tsv` + the human equivalents.

**T-2. Window-size scaling of Hi-C signal.**
*Caption:* "B/W ratio (lower = stronger enrichment) and Mantel ρ vs window size (1 Mb / 2 Mb / 4 Mb) for each meiotic stage. Monotone improvement in B/W; Mantel ρ saturates ~window-independently. 2 Mb is the practical compromise: 6 of 30 truncated arms recovered relative to 1 Mb." Built from `community_analysis_{1,2,4}Mb/50000bp/zuo2021_*_global_test.tsv`.

**T-3. PHR-size saturation curve (1 Mb → 33 Mb).**
*Caption:* "Distribution of mouse PHR sizes at each window. Saturation fraction (≥ 90 % window) drops 61 % → 0 % from 1 Mb to 33 Mb. Median stabilises at ≈ 2 Mb, identifying the natural mouse subtelomeric scale (max = 18.9 Mb at CAST chr1_p)." Built from `subtelo_{1..33}Mb/mouse.all-vs-all.p95.id95.len.tsv`.

**T-4. Identity-threshold scan at 33 Mb.**
*Caption:* "PHR length distribution at p ≥ 95, 96, 97, 98 % identity. The 95 → 98 % step removes 24 % of PHRs and cuts median size from 2.11 Mb → 0.60 Mb — outer mouse subtelomeres are moderately diverged repeats with a high-identity (≥ 98 %) 0.5–2 Mb core." Built from `subtelo_33Mb_p{96,97,98}/mouse.all-vs-all.p95.id95.len.tsv`.

**T-5. Per-arm Leiden community map (mouse vs human).**
*Caption:* "Side-by-side: 27 mouse arms in 2 communities vs 41 human arms in 15 communities. Highlights the architectural simplicity of mouse subtelomeres relative to the human duplicon system." Built from `similarity_1Mb/mouse.arm-leiden.communities.tsv` plus human equivalent.

**T-6. Mouse private-pair synteny panel.**
*Caption:* "Karyotype schematic of mouse showing private-pair pairs (CAST chr1↔chr2 at J = 0.987; B6 chr11↔chr18 at J = 0.292) connected by ribbons to their hg38 syntenic positions — all in human chromosome interiors, none in subtelomeres. Reinforces that mouse subtelomeric sharing is repeat-architectural, not syntenic." Built from mm39 → hg38 syntenic net + the §1.6 table.

**T-7. PHR vs flanking effect-size paradox (mouse).**
*Caption:* "Spearman ρ for PHR (repeat-rich) vs flanking (unique-ish) regions across 4 meiotic stages × 5 resolutions. Flanking ρ exceeds PHR ρ at all stages × resolutions — the human flanking paradox replicates in mouse, supporting the multi-mapping-suppression interpretation." Built from `seqlevel_correlation_1Mb/*_summary.tsv`.

**T-8. Repeat-class composition of mouse PHRs (open question §5).**
*Caption:* "Stacked bar of repeat-class fractions (TLC, L1-LINE, LTR, minor satellite, other) inside mouse PHR vs in flanking 100 kb, derived from RepeatMasker tracks in `mouse_T2T/annotations/`. Not yet rendered." To do — depends on resolving Gap 1 + 2.

**T-9. F1-hybrid haplotype-resolved Hi-C (open question §5).**
*Caption:* "If/when F1 (B6 × CAST) phased Hi-C is added: B6-haplotype × B6-haplotype, CAST × CAST, and B6 × CAST trans contact correlations vs sequence similarity. Tests whether subtelomeric proximity is intra-strain or trans-strain." To do.

**T-10. Methods schematic — mouse cross-species pipeline.**
*Caption:* "Mouse pipeline: PanSN B6+CAST FASTAs → telomere trim → 1 / 2 / 4 / 10 / 15 / 33 Mb flank extraction → wfmash p95 (do_not_overfilter) → impg → PHR detection → PGGB + odgi similarity → Leiden (arm) → community-based + community-free Hi-C across 4 meiotic stages × 5 resolutions × 3 windows." To do.

---

## 7. Talk slide takeaways (15-min talk)

1. **Headline.** "The sequence-homology → 3D-proximity link generalises across mammals: in mouse meiotic Hi-C (4 stages, B6 + CAST T2T) we recover the same positive correlation as in human — ρ = 0.37–0.43 at sequence level, 0.57–0.72 with community structure, p < 10⁻³⁰ at every stage."

2. **One-number recall.** "CAST chr1 ↔ CAST chr2 share 98.7 % subtelomeric Jaccard over ~990 kb — mouse's analog of the human acrocentric p-arm uniformity. B/W = 0.03 at 4 Mb means within-community contacts are 30× more enriched than between."

3. **Architectural contrast (essential slide).** "Human = 15 multi-arm communities, complex overlapping duplicon structure. Mouse = 2 communities (1 broad p-arm, 1 q-arm-enriched). Mouse subtelomeric sharing is **repeat-architecture-driven (TLC + L1-LINE)**, not syntenic — confirmed by mm39 → hg38 syntenic net showing all mouse private-pair endpoints land in human chromosome interiors." Use figure T-5 + T-6.

4. **Visual anchor.** Mouse zygotene PHR-pair scatter (P-1) next to its human counterpart — same shape, same sign, smaller community structure.

5. **Window-size sweep is informative, not arbitrary.** "We swept 1 Mb → 33 Mb. 1 Mb truncates 61 % of PHRs; 4 Mb still saturates 39 %; CAST chr1_p reaches 18.9 Mb. 2 Mb is the practical sweet spot. The B/W ratio drops monotonically as the window grows, confirming that mouse subtelomeric repeat fields really are huge." Use figure T-2 + T-3.

6. **Identity-threshold scan reveals a high-identity core.** "Raising the PHR-detection identity threshold from 95 % to 98 % cuts median PHR size from 2.11 Mb → 0.60 Mb. Mouse subtelomeres have a 0.5–2 Mb high-identity (≥ 98 %) core embedded in 1–4 Mb of moderately diverged (95–98 %) repeats." Use figure T-4.

7. **Flanking control replicates the human paradox — with caveats.** "Mouse flanking (100 kb centromere-ward of PHR) gives even higher ρ (0.60–0.77) and lower B/W (0.001–0.003) than the PHR centres themselves — same multi-mapping-suppression argument as human. But mouse flanks aren't truly unique (138–297 nonzero-Jaccard pairs persist), so the negative-control reading is weaker than in human."

8. **Method honesty.** Per-haplotype throughout (B6 = PATERNAL, CAST = MATERNAL kept separate), `wfmash do_not_overfilter` branch (no `-F 0.1` override), arm-collapsed Leiden, PHR-coord (not fixed-window) Hi-C analysis. The 500 kb pipeline was deleted because PHRs filled it; 1 Mb was the minimum that worked, 2 Mb is the recommended default.

9. **Open question to flag at end.** Repeat-annotation intersection (RepeatMasker / GFF3 from Ensembl) is the next step to quantify what fraction of the 440–495 kb shared block is L1-LINE vs LTR vs other — see Gap 1 in the survey. This will turn "TLC-architectural" into a quantitative repeat-class breakdown.
