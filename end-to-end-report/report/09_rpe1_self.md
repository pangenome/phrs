## RPE-1 self-vs-self subtelomeric pipeline

*For RPE-1 tested with HPRC population-level communities, see [Hi-C validation — RPE-1](05_hic_validation.md#cell-type-validation-rpe-1).*

### Motivation

**What it does.** Discovers RPE-1's own subtelomeric community structure from its diploid assembly, then validates it with Hi-C/Pore-C — testing whether a single individual's subtelomeric similarities predict its own 3D nuclear organization, without relying on population-level communities.

**How.** Same pipeline as mouse (the mouse section): telomere detection → 500kb flank extraction → wfmash all-vs-all → impg/PHR detection → pggb → odgi similarity → Leiden community detection → Hi-C validation. RPE-1 is a near-diploid human cell line with a known t(X;10) translocation (chrX_HAP1 carries chr10q material — a natural positive control for cross-chromosome sequence sharing).

### Flank extraction and alignment

**Key metrics.** 46 chromosomes (chr1_HAP1..chrX_HAP2), 92 flanks (46 × 2 ends), PanSN naming (RPE1#1#chr*_parm, RPE1#2#chr*_qarm). Telomeres: 0–6,329 bp trimmed per flank. wfmash `-p 95` (do_not_overfilter branch), 6,410 total alignments.

**Result.** 68/91 flanks (75%) have inter-chromosomal signal at ≥95% identity, `--min-count 1` (2 haplotypes per arm). One sequence (RPE1#2#chrX_qarm) is absent from the PHR TSV (91 of 92 flanks processed). Flanks with no signal are primarily centromere-proximal q-arm regions.

### Pangenome graph and similarity

**What it does.** Computes pairwise Jaccard similarity among all 92 flanks.

**How.** pggb on full flanks (not PHR subregions). Because RPE-1 has a single sample prefix (`RPE1#`), pggb's internal wfmash `-T RPE1#` combined with `-Y #` (exclude same-sample self-mapping) produces zero mappings. **Workaround**: wfmash run externally with `-T RPE1#1#` and `-T RPE1#2#` separately (one per haplotype), then combined PAF fed to pggb via `-a`.

**Key metrics.** 1,267 wfmash alignments, pangenome graph: 191,005 nodes, 262,658 edges, 92 paths. Similarity matrix: 8,464 pairs, 1,540 non-zero Jaccard.

### Community detection

**What it does.** Discovers arm-level communities from RPE-1's own subtelomeric similarity.

**How.** `detect_communities.R --organism human --level arm`. Leiden clustering on the arm-level Jaccard distance matrix (46 arms from 23 chromosomes × 2 arms).

**Result.** 37 communities from 46 arms. 5 multi-arm communities:

| Community | Members | Interpretation |
|---|---|---|
| C2 | **chr10_q, chrX_q** | **t(X;10) translocation detected**: chrX_HAP1 carries translocated chr10q |
| C9 | chr14_p, chr15_p, chr21_p, chr22_p | Acrocentric p-arms (NOR-bearing), same as human C7 |
| C20 | chr1_p, chr5_q, chr6_q | Cross-chromosome sharing |
| C18 | chr3_q, chr9_q, chr19_p | Cross-chromosome sharing |
| C13 | chr7_p, chr16_q | Same as human C3 (f7501 arms) |

**Key finding.** The t(X;10) translocation is independently discovered by the pipeline: chrX_q and chr10_q share subtelomeric sequence because chrX_HAP1 physically carries chr10q material. This serves as a positive control validating the method.

### Hi-C validation against self-discovered communities

**What it does.** Tests whether RPE-1's self-discovered community structure predicts 3D contacts in its own Hi-C/Pore-C data.

**How.** Three datasets (async CiFi, async Pore-C, mitotic CiFi), all 5 resolutions (5 kb–100 kb), per-haplotype (92 arms), PHR-specific coordinates. Communities from the RPE-1 self community detection. Same 3-test framework: W/B ratio, Mantel test, bootstrap permutation.

**Result.**

| Dataset | W/B ratio | Global p | Mantel rho | Mantel p | Sig communities (BH q<0.05) |
|---|---|---|---|---|---|
| Async CiFi | 83.0x | 8.2e-113 | 0.548 | <1e-300 | 19/39 |
| Async Pore-C | 45.2x | 5.5e-86 | 0.684 | <1e-300 | 20/39 |
| Mitotic CiFi | 212.8x | 2.8e-95 | 0.409 | <1e-300 | 14/39 |

**Interpretation.** The W/B ratios (45–213x) far exceed those from HPRC-community-based RPE-1 analysis (the RPE-1 validation section). This is because self-discovered communities are more specific to RPE-1's own genome — the HPRC communities are population-level averages that may not perfectly match RPE-1's specific subtelomeric architecture.

**Caveats.** Most of the 37 communities contain only 2 arms (one per haplotype from the same chromosome). For these 2-arm communities, within-community contacts are intra-chromosomal homolog contacts (e.g., chr1_HAP1_p ↔ chr1_HAP2_p), which conflates homolog pairing with subtelomeric community structure. The biologically meaningful signal comes from the 5 multi-arm communities (C2, C9, C13, C18, C20) that group arms from different chromosomes. The Mantel test, which operates on the continuous distance matrix rather than binary community labels, is not affected by this issue.

### Comparison: self-discovered vs HPRC communities

| Metric | HPRC communities (the RPE-1 validation section) | Self-discovered (the RPE-1 self Hi-C validation) |
|---|---|---|
| Communities | 15 (population-level) | 37 (RPE-1-specific) |
| Multi-arm | 15 (all) | 5 |
| W/B ratio (async CiFi) | 41.8x | 83.0x |
| Mantel rho (async CiFi) | 0.457 | 0.548 |
| Mantel rho (async Pore-C) | 0.611 | 0.684 |
| t(X;10) detected? | No (population averages) | Yes (C2: chrX_q + chr10_q) |

**Conclusion.** Self-discovered communities capture individual-specific features (t(X;10) translocation) invisible to population-level communities. The Mantel correlations are comparable, confirming that the continuous similarity-contact relationship is robust to community definition. The inflated W/B ratio in self-discovered mode reflects the dominance of intra-chromosomal 2-arm communities rather than stronger inter-chromosomal signal.

### Flanking region control (100kb centromere-ward)

**What it does.** Tests whether unique-sequence regions immediately centromere-ward of PHR boundaries also show similarity-contact correlation. This controls for multi-mapping artifacts: flanking sequences are unique (not duplicated), so any signal must reflect genuine 3D proximity, not alignment ambiguity.

**How.** For each of the 92 RPE-1 flanks, the 100kb immediately centromere-ward was extracted (`extract_flanking_sequences.py`), yielding 68 flanking sequences (24 skipped: too short or at sequence boundary). These were processed through the same pipeline: wfmash → pggb → odgi similarity → `sequence_hic_correlation.py`. Graph: 13,519 nodes, 18,391 edges, 68 paths. Only 94 non-zero Jaccard pairs (vs 1,540 for PHR flanks) and only 5 inter-chromosomal pairs with nonzero Jaccard.

**Result.**

| Dataset | Technology | Seq pairs | ρ (50kb) | p (50kb) | ρ (10kb) | p (10kb) |
|---|---|---|---|---|---|---|
| Async CiFi | CiFi | 1,177 | -0.011 | 0.717 | -0.006 | 0.828 |
| Async Pore-C | Pore-C | 1,177 | 0.127 | 1.3e-5 | 0.232 | 8.7e-16 |
| Mitotic CiFi | CiFi | 1,177 | -0.010 | 0.742 | -0.006 | 0.828 |

**Conclusion.** Flanking regions show near-zero sequence similarity across chromosomes (only 5 inter-chromosomal pairs with Jaccard > 0), so the similarity-contact correlation is absent or very weak. The Pore-C dataset shows a weak positive correlation (ρ = 0.13–0.23, driven by the few nonzero pairs), while CiFi shows none. This confirms that the strong PHR correlation (ρ = 0.30–0.44) is driven by shared subtelomeric sequence, not by a general property of chromosome-end regions. The comparison mirrors the arm-level flanking analysis (the flanking analysis section): flanking regions show community-structured 3D clustering (from the Mantel/W/B tests using arm-level aggregation), but at the sequence level the similarity signal is too sparse to drive a correlation.

### Flanking results (updated)

**Key metrics.** Flanking pggb: 4,625 similarity pairs (68 sequences). Flanking community-free correlations:

| Dataset | rho |
|---------|-----|
| Async CiFi | -0.008 |
| Async Pore-C | 0.136 |
| Mitotic CiFi | -0.007 |

**Conclusion.** Flanking community-free rho ≈ 0 across all datasets, confirming that the PHR community-free signal (the individual sequence-pair correlation section) is driven by shared subtelomeric sequence content rather than generic chromosome-end proximity effects.

### Files and scripts

**Scripts:**
| Script | Description |
|--------|-------------|
| `/moosefs/guarracino/HPRCv2/scripts/community/detect_communities.R` | Leiden community detection on RPE-1 self-similarity |
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze_hic_communities.py` | Hi-C validation (PHR coords from RPE-1 self PHR TSV) |
| `/moosefs/guarracino/HPRCv2/scripts/community/sequence_hic_correlation.py` | Community-free sequence-level (RPE-1) |

**Key files:**
| File | Description |
|------|-------------|
| `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/rpe1.all-vs-all.p95.id95.len.tsv` | RPE-1 self-discovered PHR boundaries |
| `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/rpe1_subtelo.communities.tsv` | RPE-1 self-discovered communities (37) |
| `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/rpe1_subtelo.dist_matrix.tsv` | RPE-1 arm-level distance matrix |
| `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/pggb/*similarity.tsv.gz` | RPE-1 pairwise Jaccard similarity |
| `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/hic_validation/{res}bp/` | Hi-C validation output (all resolutions) |
| `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/rpe1_*_seqlevel_summary.tsv` | Community-free correlation (all resolutions) |


---

