## Cross-species validation: mouse T2T subtelomeric analysis

### Motivation

**What it does.** Tests whether the subtelomeric sequence homology → 3D proximity link is conserved across mammals.

**How.** The same pipeline (the contig classification section–the community detection section) applied to two mouse T2T assemblies: C57BL/6J (B6, GCA_964188535) and CAST/EiJ (GCA_964188545) from Francis et al. (2025). Mouse chromosomes are telocentric — centromere at one end, telomere at the other — giving one subtelomeric region per chromosome.

### Flank extraction

**What it does.** Extracts 500 kb flanks from both chromosome ends.

**Key metrics.** B6 19 autosomes, CAST 19 autosomes + chrX. 78 total flanks. Telomeric tracts: B6 1,695–26,101 bp, CAST 3,329–21,637 bp.

**Result.** 78 flanks of 474–498 kb.

### Inter-chromosomal region detection

**What it does.** Identifies inter-chromosomal signal in mouse subtelomeric flanks.

**How.** wfmash `-p 95` (do_not_overfilter branch, cc60cd8 — k-mer frequency filtering no longer needs manual `-F 0.1` override), impg v0.4.0, same detection script with `--min-count 1`. Total alignments: 49,911.

**Result.** 39/78 flanks have signal:
- All 39 `_parm` flanks (centromere-distal/subtelomeric end) have signal.
- All 39 `_qarm` flanks (centromere-proximal end) have zero signal — satellite-dominated, confirming centromere coordinates are unnecessary.

Most subtelomeric flanks show inter-chromosomal signal spanning the **entire 500 kb window** (440–495 kb regions matching all 20 chromosomes). This is much larger than the TLC repeat itself (6–12 kb in B6 per Francis et al. 2025), indicating that the shared subtelomeric content extends well beyond the TLC motif into flanking repeat-rich sequence.

**Notable exceptions in B6:**
- B6 chr11 and chr18 form a private group with chr4 (signal only with chr4/chr11/chr18; chr11: 80 kb, chr18: 30 kb). Francis et al. (2025) noted that B6 chr11 is an exception to the standard TLC end structure — it lacks the conserved L1-LINE → TLC → minor satellite motif.
- B6 chr7 forms a small private group (chr7/chrX/chr11, 10 kb). Francis et al. also identified chr7 as a structural exception.
- B6 chr4: 15 kb signal (much smaller than the typical ~490 kb), matching 18 chromosomes (missing chr7 and chrX).

### Pangenome graph and similarity

**What it does.** Computes pairwise similarity among mouse PHR sequences.

**How.** PGGB (`-p 95 -n 2`) on 39 PHR sequences; Jaccard similarity via `odgi similarity --all`.

**Result.** Key pairwise similarities (1Mb p-arm flanks):

| Pair | Jaccard |
|------|---------|
| CAST chr1 ↔ CAST chr2 | 0.987 |
| B6 chr3 ↔ B6 chr15 | 0.944 |
| B6 chr10 ↔ B6 chr12 | 0.934 |
| CAST chr7 ↔ CAST chr16 | 0.914 |
| Most B6 cross-chr p-arm pairs | 0.88–0.94 |
| Most CAST cross-chr p-arm pairs | 0.70–0.91 |
| B6 chr11 ↔ B6 chr18 | 0.292 (outlier: B6 chr11 lacks standard TLC end) |

**Conclusion.** At the 1Mb window scale, most mouse p-arm subtelomeric flanks share extensive sequence content (Jaccard 0.7–0.99), driven by the L1-LINE + TLC repeat architecture that extends far beyond the TLC motif itself. CAST chr1/chr2 (Jaccard = 0.987) shows the highest pair identity. B6 chr11/chr18 is a low-Jaccard outlier (0.292), consistent with Francis et al.'s observation that B6 chr11 lacks the standard TLC end structure. The broadly high Jaccard across most pairs explains why Leiden finds only 2 communities — the mouse subtelomeric landscape is much more uniform than human.

### Open questions

1. **Repeat annotation**: RepeatMasker and GFF3 annotations are available for both T2T assemblies from Ensembl (https://projects.ensembl.org/mouse_genomes/). These should be downloaded and intersected with the PHR regions to determine what repeat content drives the inter-chromosomal signal — the large shared regions (440–495 kb) extend far beyond the TLC repeat (6–12 kb) and may include L1-LINE, LTR, and other repeat classes.
2. **Mouse-human synteny** (from mm39→hg38 syntenic net, UCSC): The distal (telomeric) ends of the key mouse chromosomes map to the following human regions:

| Mouse chr (distal end) | Human syntenic region | Human subtelomeric community |
|---|---|---|
| CAST chr1 (pair, J=0.987) | chr8:51.8–55.6 Mb | interior, not subtelomeric |
| CAST chr2 (pair, J=0.987) | chr10:5.9–15.4 Mb | interior, not subtelomeric |
| B6 chr11 (pair, J=0.292) | chr22:28.8–31.6 Mb | interior, not subtelomeric |
| B6 chr18 (pair, J=0.292) | chr10:35.0–35.2 Mb | interior, not subtelomeric |

The mouse private pairs do NOT correspond to human subtelomeric regions — the syntenic human positions are in the interior of human chromosomes, far from any telomere. This is expected: mouse chromosomes are telocentric (centromere at one end) while human chromosomes are metacentric, so the mouse distal telomeric end often maps to human chromosome interiors. The mouse subtelomeric inter-chromosomal sharing is driven by repeat architecture (TLC, L1-LINE), not by syntenic conservation of a subtelomeric duplicon system.
3. **Community detection**: Leiden community detection (`detect_communities.R`, `--organism mouse --level arm`, resolution scan 0.1–3.0) on the arm-level Jaccard distance matrix (collapsing B6+CAST to chromosomal arms, same as human) finds **2 communities** from 27 chromosomal arms at all three window sizes.

| Window | Communities | Structure |
|--------|------------|-----------|
| 1Mb | 2 | C1 (16 arms): chr1_p–chr3_p, chr5_p–chr6_p, chr8_p–chr10_p, chr12_p–chr17_p, chr19_p, chrX_p; C2 (11 arms): chr4_p, chr7_p, chr10_q, chr11_p, chr11_q, chr17_q, chr18_p, chr19_q, chr8_q, chr9_q, chrX_q |
| 2Mb | 2 | Same as 1Mb |
| 4Mb | 2 | C1 (17 arms): 1Mb C1 + chr11_p; C2 (10 arms): 1Mb C2 − chr11_p |

**C1** contains 16 p-arms at 1Mb. **C2** contains 11 arms enriched for q-arms (chr10_q, chr11_q, chr17_q, chr19_q, chr8_q, chr9_q, chrX_q = 7 q-arms out of 11), plus chr4_p, chr7_p, chr11_p, chr18_p. At 4Mb, chr11_p moves from C2 to C1, but the overall 2-community structure is preserved.

The 2-community structure reflects the uniform TLC-based mouse subtelomeric architecture: most arms share broadly similar repeat content (C1), while a subset with distinct q-arm sharing patterns clusters separately (C2). This is fundamentally different from human (15 multi-arm communities with complex overlapping duplicon structure across 41 arms).

### Pipeline update: 500kb → 1Mb windows

**What it does.** The 500kb subtelomeric pipeline was deleted because PHRs filled the entire 500kb window, causing flanking regions to be truncated on 30/49 per-haplotype PHR sequences. A 1Mb pipeline was completed instead.

**Key metrics.** 1Mb subtelomeric windows (`subtelo_1Mb` directory): 49 PHR regions detected. Flanking extraction: 0 truncated (the larger window provides sufficient room for both PHR and 100kb flanking regions).

**Result.** Community-free sequence-level correlation (1Mb, 50kb resolution). This is the mouse equivalent of the individual sequence-pair correlation section — each pair of PHR sequences on different chromosomes is correlated directly with Hi-C contact at their exact coordinates, without community labels:

| Meiotic stage | Seq pairs | ρ (PHR) | p (PHR) | ρ (flanking) | p (flanking) |
|--------------|-----------|---------|---------|-------------|-------------|
| Leptotene | 1,088 | 0.372 | 4.0e-37 | 0.740 | 3.4e-74 |
| Zygotene | 1,135 | **0.425** | 4.2e-51 | 0.604 | 1.2e-68 |
| Pachytene | 1,135 | **0.428** | 9.6e-52 | 0.766 | 7.2e-88 |
| Diplotene | 999 | **0.416** | 4.9e-43 | 0.715 | 1.1e-62 |

All four stages show significant positive correlation between sequence similarity and Hi-C contact for both PHR and flanking regions. The flanking correlations (ρ=0.60–0.77) are stronger than PHR (ρ=0.37–0.43), consistent with the human flanking paradox (the flanking paradox section) — unique-sequence flanking regions show stronger signal because they are not affected by multi-mapping. The 1Mb windows capture the full PHR extent (vs old 500kb: rho=0.08–0.17, non-significant).

Community-based per-arm-pair correlation (1Mb, 50kb resolution):

| Meiotic stage | rho | p-value | Mantel rho | Mantel p |
|--------------|-----|---------|------------|----------|
| Leptotene | 0.680 | 5.0e-48 | 0.687 | <0.0001 |
| Zygotene | **0.715** | 4.4e-55 | 0.718 | <0.0001 |
| Pachytene | **0.677** | 1.6e-47 | 0.683 | <0.0001 |
| Diplotene | 0.574 | 1.8e-31 | 0.577 | <0.0001 |

Flanking region control (100kb unique sequence centromere-ward of PHR, 50kb resolution):

| Meiotic stage | rho | p-value | nonzero Jaccard pairs |
|--------------|-----|---------|----------------------|
| Leptotene | 0.740 | 3.4e-74 | 155 |
| Zygotene | 0.604 | 1.2e-68 | 297 |
| Pachytene | 0.766 | 7.2e-88 | 156 |
| Diplotene | 0.715 | 1.1e-62 | 138 |

**Note on mouse flanking:** Unlike human (where flanking regions are truly unique and show weak/no correlation), mouse flanking regions retain substantial inter-chromosomal similarity (138–297 nonzero Jaccard pairs). This reflects the extent of mouse subtelomeric repeats, which span nearly the entire 1Mb window — the 100kb flanking region adjacent to the PHR boundary is still within the repeat zone. The flanking "negative control" concept does not apply straightforwardly to mouse, where TLC and L1-LINE repeats extend far beyond the PHR boundary.

### Window size optimization: 1Mb → 2Mb → 4Mb

**What it does.** Tests whether larger subtelomeric windows capture the full extent of mouse PHR regions. At 1Mb, 61% of PHRs saturate the window (>=900kb); at 2Mb, 51% still saturate (>=1.8Mb). Mouse acrocentric p-arms have massive subtelomeric repeats that extend beyond 4Mb.

**Result.** PHR size comparison across window sizes:

| Window | n PHR | Median | Mean | Max | Saturated (>=90% window) |
|--------|-------|--------|------|-----|--------------------------|
| 1Mb | 49 | 980 kb | 688 kb | 995 kb | 30/49 (61%) |
| 2Mb | 49 | 1,845 kb | 1,298 kb | 1,995 kb | 25/49 (51%) |
| 4Mb | 49 | 2,525 kb | 2,299 kb | 3,990 kb | 19/49 (39%) |

Even at 4Mb, 19/49 PHRs (39%) saturate the window. The p-arms of acrocentric chromosomes have truly massive subtelomeric repeat regions — many exceeding 4Mb. The 2Mb window represents a practical balance: it resolves 6 of the 30 arms that were truncated at 1Mb, while keeping computational costs manageable.

**Community-based Hi-C validation across window sizes (50kb resolution):**

| Window | Stage | B/W ratio | p-value | Mantel rho |
|--------|-------|-----------|---------|------------|
| 1Mb | Leptotene | 0.055 | 4.7e-35 | 0.687 |
| 1Mb | Zygotene | 0.093 | 4.0e-93 | 0.718 |
| 1Mb | Pachytene | 0.071 | 9.3e-61 | 0.683 |
| 1Mb | Diplotene | 0.046 | 2.7e-23 | 0.577 |
| 2Mb | Leptotene | 0.042 | 7.2e-38 | 0.700 |
| 2Mb | Zygotene | 0.072 | 1.2e-103 | 0.685 |
| 2Mb | Pachytene | 0.056 | 2.4e-60 | 0.680 |
| 2Mb | Diplotene | 0.038 | 3.3e-26 | 0.609 |
| 4Mb | Leptotene | 0.026 | 6.4e-100 | 0.727 |
| 4Mb | Zygotene | 0.038 | 1.1e-150 | 0.693 |
| 4Mb | Pachytene | 0.030 | 1.3e-110 | 0.711 |
| 4Mb | Diplotene | 0.034 | 1.8e-50 | 0.650 |

All analyses use per-haplotype treatment (B6 = PATERNAL, CAST = MATERNAL kept separate), chromosomal-arm-level communities (collapsing B6+CAST per arm, same as human), and PHR-specific coordinates from the per-window-size PHR TSV. B/W ratios decrease (stronger enrichment) with window size (1Mb: 0.05–0.09; 2Mb: 0.04–0.07; 4Mb: 0.03–0.04), consistent with larger windows capturing more of the sharing signal. The Mantel rho is strong and consistent across all windows (1Mb: 0.58–0.72; 2Mb: 0.61–0.70; 4Mb: 0.65–0.73).

**Window-size optimization (PHR saturation at p95 identity):**

| Window | n PHR | Median | Mean | Max | Saturated (≥90%) |
|--------|-------|--------|------|-----|-------------------|
| 1Mb | 49 | 0.98 Mb | 0.69 Mb | 1.0 Mb | 30/49 (61%) |
| 2Mb | 49 | 1.85 Mb | 1.30 Mb | 2.0 Mb | 25/49 (51%) |
| 4Mb | 49 | 2.53 Mb | 2.30 Mb | 4.0 Mb | 19/49 (39%) |
| 10Mb | 49 | 3.42 Mb | 3.92 Mb | 10.0 Mb | 8/49 (16%) |
| 15Mb | 50 | 3.44 Mb | 4.32 Mb | 15.0 Mb | 4/50 (8%) |
| 33Mb | 50 | 2.11 Mb | 3.37 Mb | 18.9 Mb | 0/50 (0%) |

At 33Mb (maximum feasible = smallest chromosome / 2), zero PHRs saturate. The largest PHR is CAST chr1_p at 18.9 Mb. The median (2.11 Mb at 33Mb) stabilizes between 10-33 Mb, confirming that most PHRs are 1-4 Mb in extent.

**Identity threshold effect at 33Mb:**

| Identity | n PHR | Median | Mean | Max |
|----------|-------|--------|------|-----|
| ≥95% | 50 | 2.11 Mb | 3.37 Mb | 18.9 Mb |
| ≥96% | 48 | 1.69 Mb | 2.64 Mb | 18.9 Mb |
| ≥97% | 42 | 1.71 Mb | 2.60 Mb | 18.9 Mb |
| ≥98% | 38 | 0.60 Mb | 1.84 Mb | 18.7 Mb |

Raising the identity threshold from 95% to 98% reduces the number of detected PHRs by 24% (50→38) and reduces the median size to less than a third (2.11→0.60 Mb), showing that the outer portions of mouse subtelomeric PHRs consist of moderately diverged (95-98%) repeats, while a core of high-identity (≥98%) sharing extends 0.5-2 Mb from the telomere.

### Mouse flanking Hi-C validation (1Mb window)

**What it does.** Tests whether unique-sequence flanking regions (100kb centromere-ward of mouse PHR boundaries) also show community-structured 3D clustering, controlling for multi-mapping artifacts in the repeat-rich mouse subtelomeres.

**How.** Same pipeline as human flanking (the flanking analysis section): extract 100kb regions immediately centromere-ward of each PHR boundary, run community-based W/B bootstrap and Mantel test at all 5 resolutions. Note: unlike human flanking regions, mouse flanking regions retain substantial inter-chromosomal similarity (see the mouse pipeline section note), so the "unique sequence negative control" interpretation is weaker for mouse.

**Result.** Flanking B/W ratios and Mantel rho across all resolutions (1Mb window):

| Stage | 5kb B/W | 10kb B/W | 20kb B/W | 50kb B/W | 100kb B/W |
|-------|---------|----------|----------|----------|-----------|
| Leptotene | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 |
| Zygotene | 0.002 | 0.002 | 0.002 | 0.002 | 0.002 |
| Pachytene | 0.003 | 0.002 | 0.002 | 0.003 | 0.003 |
| Diplotene | 0.001 | 0.001 | 0.001 | 0.001 | 0.002 |

| Stage | 5kb Mantel | 10kb Mantel | 20kb Mantel | 50kb Mantel | 100kb Mantel |
|-------|-----------|-------------|-------------|-------------|--------------|
| Leptotene | 0.656 | 0.662 | 0.597 | 0.588 | 0.589 |
| Zygotene | 0.739 | 0.750 | 0.722 | 0.643 | 0.566 |
| Pachytene | 0.712 | 0.683 | 0.662 | 0.612 | 0.622 |
| Diplotene | 0.644 | 0.622 | 0.617 | 0.495 | 0.100 (ns) |

**Conclusion.** The Mantel test shows significant positive correlations at all stages and resolutions (except diplotene at 100kb), confirming that the sequence similarity-contact relationship extends into flanking regions. The flanking B/W ratios (0.001–0.003) are even lower than PHR B/W ratios (0.03–0.12), indicating extremely strong enrichment — consistent with unique-sequence flanking regions being free from multi-mapping suppression. The signal is most consistent at fine resolutions (5-20kb). As noted in the mouse pipeline section, mouse flanking regions are not truly unique sequence — they still contain inter-chromosomal repeats — so the flanking control is less clean than in human.

**Multi-resolution mouse community-based B/W ratios (1Mb window, all 5 resolutions):**

| Stage | 5kb | 10kb | 20kb | 50kb | 100kb |
|-------|-----|------|------|------|-------|
| Leptotene | 0.073 | 0.055 | 0.057 | 0.055 | 0.029 |
| Zygotene | 0.122 | 0.112 | 0.105 | 0.093 | 0.108 |
| Pachytene | 0.119 | 0.098 | 0.091 | 0.071 | 0.076 |
| Diplotene | 0.055 | 0.050 | 0.040 | 0.046 | 0.061 |

**Multi-resolution mouse community-based B/W ratios (4Mb window, all 5 resolutions):**

| Stage | 5kb | 10kb | 20kb | 50kb | 100kb |
|-------|-----|------|------|------|-------|
| Leptotene | 0.030 | 0.022 | 0.024 | 0.026 | 0.020 |
| Zygotene | 0.057 | 0.049 | 0.044 | 0.038 | 0.042 |
| Pachytene | 0.052 | 0.037 | 0.041 | 0.030 | 0.036 |
| Diplotene | 0.024 | 0.019 | 0.023 | 0.034 | 0.043 |

All B/W ratios are well below 1.0 at every resolution for both window sizes, confirming strong community-structured 3D enrichment.

### Implications

**Conclusion.**
1. Mouse subtelomeres share sequence across chromosomes over large regions (440–495 kb), much larger than the TLC repeat alone.
2. With per-haplotype treatment and proper k-mer filtering (do_not_overfilter branch), Leiden clustering finds 2 communities at all window sizes — a broad p-arm-enriched community (C1) and a smaller q-arm-enriched community (C2) with cross-strain members.
3. The CAST chr1/chr2 pair (J=0.987) shows nearly complete subtelomeric identity across ~990 kb, analogous to human acrocentric p-arms.
4. The mouse private pairs do NOT correspond to human subtelomeric regions — syntenic human positions are in chromosome interiors. Mouse sharing is driven by repeat architecture (TLC, L1-LINE), not syntenic conservation of a duplicon system.
5. Flanking region analysis confirms the signal extends into adjacent regions, though the "unique sequence" control is less clean in mouse than in human due to the massive extent of mouse subtelomeric repeats.

### Files and scripts

**Scripts:**
| Script | Description |
|--------|-------------|
| `/moosefs/guarracino/HPRCv2/scripts/trim_mouse_flanks.py` | Trim telomeres from mouse flanks (preserves PanSN) |
| `/moosefs/guarracino/HPRCv2/scripts/find-multichr-regions-incremental.py` | PHR detection (same as human) |
| `/moosefs/guarracino/HPRCv2/scripts/community/detect_communities.R` | Leiden community detection (mouse) |
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze_hic_communities.py` | Community-based Hi-C validation (per-haplotype, PHR coords, F1 hybrid support) |
| `/moosefs/guarracino/HPRCv2/scripts/community/sequence_hic_correlation.py` | Community-free sequence-level (mouse) |
| `/moosefs/guarracino/HPRCv2/scripts/extract_flanking_sequences.py` | Flanking region extraction |

**Input data:**
| File | Description |
|------|-------------|
| `/moosefs/guarracino/HPRCv2/data/mouse_T2T/B6.PanSN.fa` | B6 (C57BL/6J) T2T assembly |
| `/moosefs/guarracino/HPRCv2/data/mouse_T2T/CAST.PanSN.fa` | CAST (CAST/EiJ) T2T assembly |
| `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/HiC/zuo2021_{stage}.mcool` | Mouse meiotic Hi-C (4 stages) |

**Output directories:**
| Directory | Description |
|-----------|-------------|
| `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/subtelo_{1,2,4}Mb/` | Subtelomeric flanks + PHR detection per window size |
| `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/subtelo_{10,15,33}Mb/` | PHR saturation scanning |
| `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/subtelo_33Mb_p{96,97,98}/` | Identity threshold scanning |
| `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/similarity_{1,2,4}Mb/` | Communities per window size |
| `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_{1,2,4}Mb/` | Hi-C validation (all res, 4 stages, per-haplotype, PHR coords) |
| `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/flanking_hic_1Mb/{res}bp/` | Mouse flanking Hi-C (all res, 4 stages) |
| `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/seqlevel_correlation_{1,2,4}Mb/` | Community-free sequence-level |
| `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/flanking_{1,2,4}Mb/` | Flanking sequences + pggb similarity |


---

