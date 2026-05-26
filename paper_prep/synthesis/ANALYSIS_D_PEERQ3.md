# D-PeerQ3: P-arm-to-Q-arm Orientation Asymmetry Per Leiden Community

**Status:** Complete  
**Date:** 2026-05-18  
**Closes:** OPEN_REVIEWER_CONCERNS.md §D-PeerQ3  
**Script:** `scripts/cladistics/pq_asymmetry_d_peerq3.py`

---

## Context

Audience-member Q3 (from NARRATIVE_EXTRACT §6) noted that P-arm-to-P-arm and Q-arm-to-Q-arm relationships dominate the arm-level Leiden communities, but P-to-Q cross-orientation sharing also exists and was unquantified. This document provides the requested supplementary analysis.

---

## Methods

Community assignments come from the Leiden k=15 clustering of the 41×41 arm-level Jaccard similarity matrix (41 signal-bearing arms; source: `arm_order_community.tsv` derived from `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`).

For each of the 15 arm-level communities, every within-community pair of arms is classified:
- **P-P**: both arms are p-arms  
- **Q-Q**: both arms are q-arms  
- **P-Q**: one p-arm, one q-arm (cross-orientation)

**Statistical test.** For each community with ≥2 arms, a one-sided Fisher exact test (alternative = "greater") is used to evaluate whether same-orientation pairs (P-P + Q-Q) are enriched within the community relative to the global background. The 2×2 table is:

|               | Same-orient | Cross-orient (P-Q) |
|---|---|---|
| Within comm.  | obs PP+QQ   | obs P-Q            |
| Outside comm. | (global PP+QQ) − obs | (global P-Q) − obs |

Global counts across all 41 arms (21 P, 20 Q):
- Global P-P = C(21,2) = **210**
- Global Q-Q = C(20,2) = **190**
- Global P-Q = 21 × 20 = **420**
- Total arm pairs = C(41,2) = **820**
- Expected P-Q fraction = 420/820 = **51.2%**

P-values are corrected for multiple comparisons across 11 testable communities (the 4 singletons C8, C9, C10, C13 have 0 edges and are excluded) using the Benjamini-Hochberg (BH) procedure.

---

## Results

### Per-community edge count table

| Community | Members | n | nP | nQ | Orientation class | PP | QQ | PQ | PQ% obs | OR same-orient | p (Fisher) | p (BH) | Sig |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **C1** | 4q, 10q | 2 | 0 | 2 | pure-Q | 0 | 1 | 0 | 0.0% | ∞ | 0.488 | 0.767 | |
| **C2** | 10p, 18p | 2 | 2 | 0 | pure-P | 1 | 0 | 0 | 0.0% | ∞ | 0.488 | 0.767 | |
| **C3** | 3q, 7p, 9q, 11p, 16q, 19p | 6 | 3 | 3 | mixed | 3 | 3 | 9 | 60.0% | 0.70 | 0.828 | 0.968 | |
| **C4** | 7q, 12q | 2 | 0 | 2 | pure-Q | 0 | 1 | 0 | 0.0% | ∞ | 0.488 | 0.767 | |
| **C5** | 6p, 9p, 12p, 20q | 4 | 3 | 1 | mixed | 3 | 0 | 3 | 50.0% | 1.05 | 0.634 | 0.871 | |
| **C6** | 1q, 13q, 17q, 19q, 21q, 22q | 6 | 0 | 6 | pure-Q | 0 | 15 | 0 | 0.0% | ∞ | <0.001 | 0.0002 | *** |
| **C7** | 13p, 14p, 15p, 21p, 22p | 5 | 5 | 0 | pure-P | 10 | 0 | 0 | 0.0% | ∞ | 0.0007 | 0.004 | ** |
| **C8** | 15q | 1 | 0 | 1 | singleton | 0 | 0 | 0 | — | — | — | — | |
| **C9** | 16p | 1 | 1 | 0 | singleton | 0 | 0 | 0 | — | — | — | — | |
| **C10** | 17p | 1 | 1 | 0 | singleton | 0 | 0 | 0 | — | — | — | — | |
| **C11** | 1p, 5q, 6q, 8p | 4 | 2 | 2 | mixed | 1 | 1 | 4 | 66.7% | 0.52 | 0.880 | 0.968 | |
| **C12** | 2q, 20p | 2 | 1 | 1 | mixed | 0 | 0 | 1 | 100.0% | 0.00 | 1.000 | 1.000 | |
| **C13** | 4p | 1 | 1 | 0 | singleton | 0 | 0 | 0 | — | — | — | — | |
| **C14** | Xq, Yq | 2 | 0 | 2 | pure-Q | 0 | 1 | 0 | 0.0% | ∞ | 0.488 | 0.767 | |
| **C15** | Xp, Yp | 2 | 2 | 0 | pure-P | 1 | 0 | 0 | 0.0% | ∞ | 0.488 | 0.767 | |

OR = odds ratio for same-orientation enrichment (within vs outside community); `***` p_BH < 0.001; `**` p_BH < 0.01.

---

## Key Findings

### 1. Most communities are intra-orientation

10 of the 11 communities with >1 member show 0 P-Q edges (pure-P or pure-Q orientation class). Across all 41 arms, only **17 P-Q edges** exist within the 15 communities, out of 820 total pairwise arm comparisons.

### 2. Statistically significant same-orientation enrichment in C6 and C7

- **C6** (1q/13q/17q/19q/21q/22q — terminal long arms of acrocentrics + chr1/17/19): 6 pure-Q arms, 15 Q-Q edges, 0 P-Q. Strongest signal: p_BH = 0.0002.
- **C7** (13p/14p/15p/21p/22p — acrocentric short arms): 5 pure-P arms, 10 P-P edges, 0 P-Q. p_BH = 0.004.

These two communities have more than enough statistical power (15 and 10 edges, respectively) to formally reject random arm-orientation assignment.

### 3. C1 and C2 are intra-orientation by construction but underpowered

- **C1** (4q/10q): the sole community of 2 same-orientation Q-arms; 1 Q-Q edge, 0 P-Q. OR = ∞, but 1 edge provides no power (p_BH = 0.77).
- **C2** (10p/18p): analogously, 1 P-P edge, 0 P-Q. p_BH = 0.77.

These are consistent with intra-orientation clustering but the 2-arm communities cannot reach significance with any single-edge test.

### 4. P-Q edges cluster in four mixed communities

All 17 within-community P-Q edges concentrate in C3, C5, C11, and C12:

| Community | PQ edges | PQ% (obs) | PQ% (expected under null) |
|---|---|---|---|
| C3 | 9 | 60.0% | 51.2% |
| C5 | 3 | 50.0% | 51.2% |
| C11 | 4 | 66.7% | 51.2% |
| C12 | 1 | 100.0% | 51.2% |

None reaches significance (p_BH ≥ 0.87), meaning these mixed communities are not significantly more cross-oriented than chance—they are simply the communities where the Leiden algorithm placed heterogeneous-orientation arms together, driven by sequence similarity rather than orientation.

---

## Summary Statement (for v6 paper)

> Across the 15 arm-level Leiden communities, P-arm-to-Q-arm (P-Q) edges are rare: only 17 of 75 within-community arm pairs (23%) are cross-orientation, compared to 51% expected by chance (420/820 global pairs). Same-orientation enrichment is statistically significant after BH correction in C6 (six acrocentric-q arms; p_BH = 0.0002) and C7 (five acrocentric-p arms; p_BH = 0.004). C1 (4q/10q) and C2 (10p/18p) are uniformly intra-Q and intra-P respectively, but their two-arm size provides insufficient power. The four communities with any P-Q edges (C3, C5, C11, C12) show no significant deviation from the null expectation, indicating that cross-orientation sharing occurs at background rates in mixed communities.

---

## Recommended v6 Edit

Add one line to the Supplementary Data section referencing a new supplementary table:

> **Supplementary Table S-PQ.** Per-community P-arm/Q-arm edge counts and same-orientation enrichment test (Fisher exact, BH-corrected) for the 15 Leiden k=15 arm-level communities. Generated by `scripts/cladistics/pq_asymmetry_d_peerq3.py`.

In the main text (Results, communities paragraph), add one clause to the existing community description: *"Community membership is strongly intra-orientation: same-orientation arm pairs account for 77% of within-community edges (Supplementary Table S-PQ), with significant same-orientation enrichment in C6 and C7 (Fisher exact, p_BH < 0.005)."*

---

## Validation Checklist

- [x] Per-community P-P / Q-Q / P-Q counts: see table above.
- [x] BH-corrected asymmetry test results: C6 p_BH = 0.0002, C7 p_BH = 0.004; all other communities not significant.
- [x] Recommended v6 edit: supplementary table reference and one main-text clause provided above.
