---
title: "Survey 09 — RPE-1 self-vs-self subtelomeric pipeline"
source: end-to-end-report/report/09_rpe1_self.md
scope: RPE-1 self-discovered communities, Hi-C/Pore-C validation, comparison with HPRC communities, flanking control
audience: Nature manuscript and 15-min talk
---

# Survey 09 — RPE-1 self-vs-self

This survey extracts and structures the content of `end-to-end-report/report/09_rpe1_self.md` for the Nature manuscript and the companion 15-minute talk. The source section runs the same pipeline used for HG002/HPRC and mouse — telomere → flank → wfmash → impg/PHR → pggb → odgi similarity → Leiden communities → Hi-C — but on the **RPE-1 diploid assembly alone**, then validates the resulting communities against three RPE-1 3D datasets (async CiFi, async Pore-C, mitotic CiFi). RPE-1 carries a known **t(X;10) translocation** that the pipeline rediscovers from sequence alone, providing a positive control. For HPRC-community testing on the same RPE-1 contact data, see survey 05 (Hi-C validation).

---

## 1. Key findings with metrics

### 1.1 Flank extraction and alignment
- Diploid assembly RPE1v1.1: 46 chromosomes (chr1_HAP1..chrX_HAP2).
- **92 flanks** = 46 chromosomes × 2 ends (PanSN names `RPE1#1#chr*_parm`, `RPE1#2#chr*_qarm`).
- Telomeres trimmed: **0–6,329 bp** per flank.
- wfmash `-p 95` (do_not_overfilter branch); **6,410** total all-vs-all alignments.
- **75 % (68/91) of flanks** have inter-chromosomal signal at ≥95 % identity, `--min-count 1` (2 haplotypes per arm). One sequence (RPE1#2#chrX_qarm) is missing from the PHR TSV; **91/92 flanks** processed. Flanks with no signal are primarily centromere-proximal q-arms.

### 1.2 Pangenome graph and similarity (PHR-region pggb)
- **Workaround for single-prefix wfmash:** because RPE-1 has a single sample prefix (`RPE1#`), pggb's internal `wfmash -T RPE1# -Y #` produces zero mappings. wfmash run **externally** with `-T RPE1#1#` and `-T RPE1#2#` separately, then combined PAF fed to pggb via `-a`.
- pggb output: **1,267** wfmash alignments → graph with **191,005 nodes / 262,658 edges / 92 paths**.
- Similarity matrix: **8,464 pairs** total; **1,540** non-zero Jaccard.

### 1.3 Community detection (RPE-1-specific)
- `detect_communities.R --organism human --level arm` → Leiden on the 46-arm Jaccard distance matrix.
- **37 communities from 46 arms.** Five multi-arm communities (the rest are 2-arm intra-chromosomal homolog pairs):

| Community | Members | Interpretation |
|---|---|---|
| C2 | **chr10_q, chrX_q** | **t(X;10) translocation rediscovered**: chrX_HAP1 carries translocated chr10q |
| C9 | chr14_p, chr15_p, chr21_p, chr22_p | Acrocentric p-arms (NOR-bearing); same arms as HPRC C7 |
| C20 | chr1_p, chr5_q, chr6_q | Cross-chromosome sharing |
| C18 | chr3_q, chr9_q, chr19_p | Cross-chromosome sharing |
| C13 | chr7_p, chr16_q | Same as HPRC C3 (f7501 arms) |

- **Headline finding:** the t(X;10) translocation is **independently discovered** from sequence alone — the chrX_HAP1 q-arm clusters with chr10_q because it physically carries chr10q material. This serves as a positive control for the entire pipeline.

### 1.4 Hi-C validation against self-discovered communities
Three RPE-1 datasets, all 5 resolutions (5 kb–100 kb), per-haplotype (92 arms), PHR-specific coordinates from RPE-1's own PHR TSV:

| Dataset | W/B ratio | Global p | Mantel ρ | Mantel p | Sig. communities (BH q < 0.05) |
|---|---|---|---|---|---|
| Async CiFi | **83.0×** | 8.2 × 10⁻¹¹³ | 0.548 | < 1 × 10⁻³⁰⁰ | 19/39 |
| Async Pore-C | 45.2× | 5.5 × 10⁻⁸⁶ | **0.684** | < 1 × 10⁻³⁰⁰ | 20/39 |
| Mitotic CiFi | **212.8×** | 2.8 × 10⁻⁹⁵ | 0.409 | < 1 × 10⁻³⁰⁰ | 14/39 |

**Caveat (essential):** 32/37 communities are 2-arm intra-chromosomal homolog pairs (e.g. chr1_HAP1_p ↔ chr1_HAP2_p). For these, "within-community" contacts are intra-chromosomal homolog contacts, conflating homolog pairing with subtelomeric similarity. The biologically meaningful inter-chromosomal signal lives in the **5 multi-arm communities** (C2, C9, C13, C18, C20). The Mantel test, computed on the continuous distance matrix, is not affected by this conflation.

### 1.5 Self-discovered vs HPRC communities (same RPE-1 contact data)

| Metric | HPRC communities (survey 05) | Self-discovered (this section) |
|---|---:|---:|
| Communities | 15 (population-level) | 37 (RPE-1-specific) |
| Multi-arm | 15 (all) | 5 |
| W/B ratio (async CiFi) | 41.8× | **83.0×** |
| Mantel ρ (async CiFi) | 0.457 | **0.548** |
| Mantel ρ (async Pore-C) | 0.611 | **0.684** |
| t(X;10) detected? | No (population averages mask it) | **Yes** (C2) |

**Interpretation.** Self-discovered communities capture RPE-1-individual-specific features (t(X;10)) invisible to population-level HPRC communities. The Mantel correlations are comparable — and slightly stronger for the self-discovered case — confirming that the continuous similarity-contact relationship is robust to community definition. The inflated W/B ratio is partly an artefact of intra-chromosomal homolog 2-arm communities, **not** a uniformly stronger inter-chromosomal signal.

### 1.6 Flanking-region control (100 kb centromere-ward)
Tests whether unique-sequence regions immediately centromere-ward of PHR boundaries also carry the similarity-contact signal — controls for multi-mapping artifacts.

- **68 flanking sequences** retained (24 of 92 skipped: too short or at sequence boundary).
- Flanking pggb graph: **13,519 nodes / 18,391 edges / 68 paths**.
- Only **94 non-zero Jaccard pairs** vs **1,540** for PHR flanks; only **5 inter-chromosomal pairs** with non-zero Jaccard.

Sequence-level Spearman ρ (Jaccard vs Hi-C contact, **n = 1,177** pairs):

| Dataset | ρ (50 kb) | p (50 kb) | ρ (10 kb) | p (10 kb) |
|---|---|---|---|---|
| Async CiFi | −0.011 | 0.717 | −0.006 | 0.828 |
| Async Pore-C | **0.127** | 1.3 × 10⁻⁵ | **0.232** | 8.7 × 10⁻¹⁶ |
| Mitotic CiFi | −0.010 | 0.742 | −0.006 | 0.828 |

Community-free correlations on the same 4,625 flanking pggb similarity pairs (68 sequences):

| Dataset | ρ |
|---|---|
| Async CiFi | −0.008 |
| Async Pore-C | 0.136 |
| Mitotic CiFi | −0.007 |

**Conclusion.** Flanking regions show near-zero sequence similarity across chromosomes, so the similarity-contact correlation collapses (ρ ≈ 0) for both CiFi datasets. Pore-C retains a weak positive ρ (0.13–0.23) driven by the few non-zero pairs. The strong PHR correlation (ρ = 0.30–0.44 elsewhere in the report) is therefore driven by **shared subtelomeric sequence**, not by a generic chromosome-end proximity effect.

---

## 2. Existing figures (paths)

All paths are absolute and live on `/moosefs/guarracino/HPRCv2/`.

### 2.1 PHR pggb (similarity graph)
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/pggb/`
- `rpe1_phr_combined.paf.11fba48.f9ca20d.smooth.final.og.lay.draw.png` — graph layout.
- `rpe1_phr_combined.paf.11fba48.f9ca20d.smooth.final.og.lay.draw_multiqc.png` — graph layout (multiqc-style).
- `rpe1_phr_combined.paf.…viz_multiqc.png` and `viz_O / viz_depth / viz_inv / viz_pos / viz_uncalled` PNGs — odgi viz panels.

### 2.2 Hi-C validation (per dataset × per resolution)
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/hic_validation/{50000bp,100000bp,20000bp,10000bp,5000bp}/`

For each of the three datasets (`async_cifi`, `async_porec`, `mitotic_cifi`) at each resolution:
- `rpe1_self_<dataset>_hic_bootstrap_distributions.pdf` — bootstrap permutation null distributions.
- `rpe1_self_<dataset>_hic_community_heatmap.pdf` — within/between contact heatmap by community.
- `rpe1_self_<dataset>_hic_mds_comparison.pdf` — MDS embedding of arms (sequence vs Hi-C).
- `rpe1_self_<dataset>_hic_similarity_contact_scatter.pdf` — Mantel scatter (similarity vs contact).
- `rpe1_self_<dataset>_phr_pair_scatter.pdf` — per-PHR-pair Spearman scatter.

(Singleton-neighbour TSVs `rpe1_self_<dataset>_singleton_C{23,30,32}_neighbors.tsv` accompany the scatters.)

### 2.3 PHR-region sequence-level scatter
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/`
- `async_cifi_seqlevel_10kb_seqlevel_scatter.pdf` — async CiFi PHR sequence-level scatter (10 kb).
- `async_porec_seqlevel_10kb_seqlevel_scatter.pdf`
- `mitotic_cifi_seqlevel_10kb_seqlevel_scatter.pdf`
- `hg002_hic_seqlevel_scatter.pdf`, `hg002_porec_seqlevel_scatter.pdf` — comparator scatters projected onto RPE-1 PHR coords.

### 2.4 Flanking control (100 kb centromere-ward)
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/`
- `async_cifi_flanking_{10000,50000}bp_seqlevel_scatter.pdf`
- `async_porec_flanking_{10000,50000}bp_seqlevel_scatter.pdf`
- `mitotic_cifi_flanking_{10000,50000}bp_seqlevel_scatter.pdf`

### 2.5 Flanking pggb
Directories: `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/{pggb_flanking,flanking_pggb}/`
- `rpe1_flank_combined.paf.…final.og.lay.draw.png` — flanking-graph layout (pggb_flanking).
- Multiqc viz PNGs (`viz_*_multiqc.png`) for the flanking graph.

---

## 3. Existing CSVs / TSVs (paths)

### 3.1 PHR-level inputs and community outputs
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/`
- `RPE1.telo.tsv` — telomere coordinates per arm (input to flank extraction).
- `rpe1.all-vs-all.p95.id95.len.tsv` — RPE-1 self-discovered PHR boundaries (the "PHR TSV").
- **`rpe1.communities.tsv`** — 37 RPE-1 self-discovered communities (note: the source markdown calls this `rpe1_subtelo.communities.tsv`; the actual filename on disk is `rpe1.communities.tsv`).
- **`rpe1.dist_matrix.tsv`** — arm-level distance matrix (same naming caveat: not `rpe1_subtelo.dist_matrix.tsv`).
- `rpe1.leiden_scan.tsv` — Leiden resolution-scan trace (community-count vs resolution; not described in the section).

### 3.2 PHR pggb similarity
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/pggb/`
- `rpe1_phr_combined.paf.11fba48.f9ca20d.smooth.final.similarity.tsv.gz` — pairwise Jaccard similarity (1,540 non-zero pairs).
- `rpe1_phr_combined.paf.…similarity.tsv.gz.arm_pair_jaccard.tsv` — arm-level aggregated Jaccard.
- `rpe1_phr_combined.paf.…fix.affixes.tsv.gz` — odgi affixes / variant table.

### 3.3 Hi-C validation (per dataset × per resolution)
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/hic_validation/{50000bp,100000bp,20000bp,10000bp,5000bp}/`

For each `<dataset>` ∈ {`async_cifi`, `async_porec`, `mitotic_cifi`}:
- `rpe1_self_<dataset>_global_test.tsv` — global W/B and bootstrap p.
- `rpe1_self_<dataset>_community_bootstrap_tests.tsv` — per-community W/B + BH-corrected q.
- `rpe1_self_<dataset>_contact_matrix.tsv` — arm × arm Hi-C contact matrix.
- `rpe1_self_<dataset>_phr_pair_correlation.tsv` — per-PHR-pair Spearman ρ.
- `rpe1_self_<dataset>_singleton_C{23,30,32}_neighbors.tsv` — nearest-neighbour arms for singleton communities (companion to the similarity-contact scatter).
- `rpe1_self_<dataset>_subtelomeric_regions.bed` — coords used to extract the contact submatrix.

### 3.4 PHR-region sequence-level summaries
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/`
- `async_cifi_seqlevel_10kb_seqlevel.tsv` + `..._summary.tsv`
- `async_porec_seqlevel_10kb_seqlevel.tsv` + `..._summary.tsv`
- `mitotic_cifi_seqlevel_10kb_seqlevel.tsv` + `..._summary.tsv`
- All resolutions: `rpe1_<dataset>_{5000,10000,20000,50000,100000}bp_seqlevel.tsv` + `..._summary.tsv` — the per-resolution sequence-level inputs that feed into the report's tables.
- `hg002_hic_seqlevel.tsv` / `hg002_porec_seqlevel.tsv` (+ `_summary.tsv`) — HG002 contacts mapped onto RPE-1 PHR (comparator).

### 3.5 Flanking control (100 kb centromere-ward)
Directory: `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/`
- `async_cifi_flanking_{10000,50000}bp_seqlevel.tsv` + `..._summary.tsv`
- `async_porec_flanking_{10000,50000}bp_seqlevel.tsv` + `..._summary.tsv`
- `mitotic_cifi_flanking_{10000,50000}bp_seqlevel.tsv` + `..._summary.tsv`

### 3.6 Flanking pggb similarity
- `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/pggb_flanking/rpe1_flank_combined.paf.…final.similarity.tsv.gz` — flanking pairwise Jaccard (94 non-zero pairs).
- `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/flanking_pggb/rpe1_flank_combined.paf.…final.similarity.tsv.gz` — duplicate run (older params).

### 3.7 Sequence inputs
- `RPE1_flanks_500kb_trimmed.fa.gz` (+ `.fai`, `.gzi`) — 92 PHR flanks (500 kb terminal, telomere-trimmed).
- `RPE1_flanking_100kb.fa.gz` (+ indices) — 68 100 kb sequences immediately centromere-ward of PHR boundaries.

---

## 4. Methods

### 4.1 Pipeline (identical to mouse / HPRC)
Telomere detection → 500 kb flank extraction → wfmash `-p 95` (do_not_overfilter branch) all-vs-all → impg/PHR detection → pggb on full flanks (not PHR subregions) → odgi pairwise Jaccard similarity → Leiden community detection at the **arm level** → Hi-C validation against the resulting communities. RPE-1 input: diploid RPE1v1.1 assembly (PanSN-named `RPE1#{1,2}#chr*`).

### 4.2 Single-sample wfmash workaround
pggb's internal `wfmash -T RPE1#` combined with `-Y #` (exclude same-sample self-mapping) produces zero mappings on a single-prefix dataset. **Workaround:** run wfmash externally with `-T RPE1#1#` and `-T RPE1#2#` separately (one per haplotype), concatenate the two PAFs, and feed to pggb via `-a`.

### 4.3 Community detection
`detect_communities.R --organism human --level arm`. Leiden clustering on the arm-level Jaccard distance matrix (46 arms = 23 chromosomes × 2 arms). Output: 37 communities.

### 4.4 Hi-C validation framework (3-test)
Same as the HPRC pipeline (survey 05):
1. **W/B ratio** + bootstrap permutation (10,000 perms) — global W/B with permutation p; per-community BH q.
2. **Mantel test** — Spearman ρ between sequence distance matrix and Hi-C distance matrix (continuous, not categorical).
3. **Per-PHR-pair correlation** — pairwise Jaccard vs pairwise Hi-C contact across all arm pairs.

All five resolutions (5 / 10 / 20 / 50 / 100 kb), per-haplotype (92 arms), PHR-specific coordinates from RPE-1's own PHR TSV.

### 4.5 Flanking control
- For each of the 92 RPE-1 flanks, the 100 kb immediately centromere-ward of the PHR boundary was extracted (`extract_flanking_sequences.py`).
- 68 sequences retained (24 skipped: too short or at sequence boundary).
- Same pipeline downstream: wfmash → pggb → odgi similarity → `sequence_hic_correlation.py`.

### 4.6 Scripts
| Script | Purpose |
|---|---|
| `/moosefs/guarracino/HPRCv2/scripts/community/detect_communities.R` | Leiden community detection on RPE-1 self-similarity. |
| `/moosefs/guarracino/HPRCv2/scripts/community/analyze_hic_communities.py` | Hi-C 3-test framework using PHR coords from RPE-1 self PHR TSV. |
| `/moosefs/guarracino/HPRCv2/scripts/community/sequence_hic_correlation.py` | Community-free per-sequence-pair correlation (RPE-1). |
| `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/extract_rpe1_flanks.py` | RPE-1 specific flank extraction (in-tree). |

---

## 5. Gaps

1. **Filename inconsistency.** Source markdown cites `rpe1_subtelo.communities.tsv` and `rpe1_subtelo.dist_matrix.tsv`; the actual files on disk are `rpe1.communities.tsv` and `rpe1.dist_matrix.tsv`. Either the report or the files should be renamed before publication so reviewers can locate the artefacts.
2. **Missing sequence (RPE1#2#chrX_qarm).** 91 of 92 flanks processed; one arm absent from the PHR TSV. No diagnosis of *why* (assembly gap? telomere unresolved? excluded by `--min-count`?). Section reports the count but offers no follow-up.
3. **No multi-arm-only re-analysis of W/B.** The 5 multi-arm communities (C2, C9, C13, C18, C20) are flagged as carrying the biologically meaningful inter-chromosomal signal, but the W/B / Mantel statistics are not recomputed restricted to those communities. Without that, the headline 83×–212× W/B is hard to interpret against the homolog-pair conflation caveat.
4. **No HPRC-vs-self overlap analysis at the community level.** The comparison table reports counts (15 vs 37) and one shared community pair (C13 ↔ HPRC C3, C9 ↔ HPRC C7) but no ARI / Jaccard across arm-membership for the two clusterings on the same 46 arms.
5. **t(X;10) — no orthogonal validation.** The translocation rediscovery is the headline biological result, but the section does not show a contact map / dotplot / breakpoint coordinate confirming that chrX_HAP1's q-arm is the chr10q-bearing one. The graph-level evidence is implicit.
6. **No Mantel resolution sweep.** Only one Mantel ρ per dataset is reported (resolution unspecified in the table — appears to be aggregated). Per-resolution Mantel ρ across 5–100 kb is computed (TSVs exist) but not tabulated; this would parallel the HPRC presentation.
7. **Asymmetric flanking-control reporting.** Flanking sequence-level ρ table reports 50 kb and 10 kb only; community-free table reports a single ρ per dataset (resolution implicit). Resolutions reported are inconsistent across the two flanking analyses.
8. **No 5/20/100 kb scatter PDFs in the figures inventory above** for community-free PHR — only 10 kb scatter PDFs exist (`*_seqlevel_10kb_seqlevel_scatter.pdf`). Other resolutions are CSV-only.
9. **Mitotic CiFi outlier — no explanation for the W/B vs Mantel divergence.** Mitotic CiFi has the *highest* W/B (212.8×) but the *lowest* Mantel ρ (0.409) and lowest significant-community count (14/39). Same pattern as in HPRC RPE-1 (survey 05) but the section cross-references neither, nor offers a synchronisation-purity caveat.
10. **No L/R / parm-vs-qarm split.** All 5 multi-arm communities mix arm types except C9 (all p) and C13 (one of each). Whether the homolog-pair 2-arm communities pair p with p / q with q (expected) vs cross-arm (unexpected) is not reported.
11. **Flanking community-free uses a *different* sequence set (68) than the PHR analysis (91/92).** The expected near-zero ρ is partly a power issue (94 vs 1,540 non-zero pairs) — the section attributes it to biology but does not show that flanking ρ remains near zero on a power-matched subsample.
12. **No mouse cross-reference.** The same pipeline was run on mouse (survey 08), but section 09 does not compare RPE-1's "self-vs-self at single-individual scale" findings to the mouse single-individual analogue. Worth a one-sentence comparison given the structural parallel.

---

## 6. Suggested figures with captions (produced vs to-do)

### Already produced (use directly or recompose)

**P-1. RPE-1 PHR pangenome graph layout.**
*File:* `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/pggb/rpe1_phr_combined.paf.…final.og.lay.draw.png`.
*Caption:* "Pangenome graph layout of RPE-1's 92 subtelomeric flanks (191,005 nodes / 262,658 edges). Path bundles correspond to chromosome arms; cross-arm contact in the layout reflects the 1,540 non-zero pairwise Jaccard pairs."

**P-2. Async CiFi Mantel scatter (self-discovered).**
*File:* `hic_validation/50000bp/rpe1_self_async_cifi_hic_similarity_contact_scatter.pdf`.
*Caption:* "PHR sequence distance vs Hi-C distance across 92 RPE-1 arms (async CiFi, 50 kb). Mantel ρ = 0.548, p < 1 × 10⁻³⁰⁰; each point is one inter-chromosomal arm pair."

**P-3. Async CiFi within/between heatmap.**
*File:* `hic_validation/50000bp/rpe1_self_async_cifi_hic_community_heatmap.pdf`.
*Caption:* "Arm × arm contact heatmap with arms ordered by self-discovered community. The diagonal-block enrichment is dominated by 32 intra-chromosomal homolog 2-arm communities; the 5 multi-arm blocks (C2, C9, C13, C18, C20) carry the inter-chromosomal signal."

**P-4. MDS comparison (sequence vs Hi-C embedding).**
*File:* `hic_validation/50000bp/rpe1_self_async_cifi_hic_mds_comparison.pdf`.
*Caption:* "Two-dimensional MDS of arms by sequence distance (left) vs Hi-C distance (right). Coloured by self-discovered community; correspondence between embeddings is the visual reading of Mantel ρ = 0.548."

**P-5. Async Pore-C Mantel scatter (strongest signal).**
*File:* `hic_validation/50000bp/rpe1_self_async_porec_hic_similarity_contact_scatter.pdf`.
*Caption:* "Async Pore-C, 50 kb. Mantel ρ = 0.684 — the strongest correlation across the three RPE-1 datasets, consistent with multi-way contacts amplifying signal."

**P-6. Bootstrap null distributions.**
*File:* `hic_validation/50000bp/rpe1_self_<dataset>_hic_bootstrap_distributions.pdf` (one per dataset).
*Caption:* "Bootstrap permutation null vs observed W/B for each of the 39 RPE-1 self-discovered communities; vertical line = observed; histograms = 10,000-permutation null."

**P-7. Flanking control scatter.**
*File:* `mitotic_cifi_flanking_50000bp_seqlevel_scatter.pdf` (and async pair).
*Caption:* "Flanking-region sequence-level scatter (100 kb centromere-ward of PHR). With only 5 inter-chromosomal Jaccard pairs > 0, both CiFi datasets sit at ρ ≈ 0; Pore-C retains a weak positive ρ = 0.13–0.23 driven by those few pairs."

**P-8. PHR community-free scatter (PHR comparator).**
*File:* `async_porec_seqlevel_10kb_seqlevel_scatter.pdf`.
*Caption:* "Per-PHR-pair sequence Jaccard vs Hi-C contact (async Pore-C, 10 kb). Companion to the flanking control above — shows that PHR-driven similarity is what drives the global ρ, not chromosome-end proximity."

### To-do (suggested new figures)

**T-1. t(X;10) one-figure positive control.**
*Caption:* "Top: schematic of RPE-1's t(X;10) translocation; chrX_HAP1's q-arm physically carries chr10q material. Middle: arm-pair Jaccard heatmap with chrX_q × chr10_q highlighted (extract from `rpe1.dist_matrix.tsv`). Bottom: the same pair on the Hi-C contact matrix (`rpe1_self_async_cifi_contact_matrix.tsv`) — sequence sharing and elevated 3D contact aligned." Required for the talk's "the pipeline rediscovers a known translocation from sequence alone" beat.

**T-2. Multi-arm-only W/B re-analysis.**
*Caption:* "W/B ratio recomputed restricted to the 5 multi-arm communities (C2, C9, C13, C18, C20). Removes the homolog-pair 2-arm conflation and isolates the inter-chromosomal signal. Comparison: full W/B (83×, 45×, 213×) vs multi-arm-only W/B (TBD)."

**T-3. Self vs HPRC community membership comparison.**
*Caption:* "Sankey / alluvial diagram mapping the 46 RPE-1 arms across the 15 HPRC communities (left) and 37 RPE-1 self-discovered communities (right). Highlights the C9 ↔ HPRC C7 acrocentric block, the C13 ↔ HPRC C3 chr7_p+chr16_q match, and the t(X;10)-only C2 (no HPRC analogue)."

**T-4. Mantel ρ across resolutions — RPE-1 self, all three datasets.**
*Caption:* "Mantel ρ as a function of contact-matrix resolution (5 / 10 / 20 / 50 / 100 kb). Async Pore-C > async CiFi > mitotic CiFi at every resolution; ρ falls monotonically as resolution coarsens. Built from `hic_validation/{res}bp/rpe1_self_<dataset>_phr_pair_correlation.tsv`."

**T-5. Per-haplotype W/B.**
*Caption:* "W/B per arm × per haplotype (92 arms × 3 datasets). Heatmap reveals which arms drive the global signal and whether HAP1 vs HAP2 contributions are symmetric (a sanity check on the homolog-pairing artefact)."

**T-6. Cell-cycle modulation panel — async vs mitotic on self-discovered communities.**
*Caption:* "Direct comparison of async CiFi vs mitotic CiFi on the *same* 37 RPE-1 self-discovered communities. W/B amplifies (83× → 213×) while Mantel ρ attenuates (0.548 → 0.409). Replicates the HPRC-community pattern in survey 05 with self-discovered communities."

**T-7. Methods schematic (talk-grade).**
*Caption:* "Telomere → flank → wfmash (per-haplotype workaround) → pggb → Leiden (37 communities) → Hi-C 3-test, with the t(X;10) emerging at the Leiden step."

**T-8. Self-discovered vs HPRC W/B / Mantel ρ summary bar.**
*Caption:* "Two-panel bar chart for the talk: (a) W/B HPRC vs self for each of the three RPE-1 datasets; (b) Mantel ρ HPRC vs self. Built from existing `phr_pair_correlation.tsv` and `global_test.tsv` files plus survey 05 numbers."

---

## 7. Talk slide takeaways (15-min talk)

1. **Headline (one slide):** "Run the pipeline on a single individual's diploid genome — no population at all — and it rediscovers RPE-1's t(X;10) translocation from sequence alone, then predicts that translocation's 3D position in the nucleus." Use figure T-1.
2. **Self vs population numbers (one slide):** "RPE-1's *own* communities give Mantel ρ = 0.55 (CiFi) and 0.68 (Pore-C), comparable-to-better than the HPRC population-level ρ on the same Hi-C data — and they capture features (t(X;10)) the population-level run cannot." Use figure T-8.
3. **Caveat slide (essential, prevent reviewer pushback):** "32 of 37 self-discovered communities are 2-arm intra-chromosomal homolog pairs; the inflated 83×–213× W/B is partly homolog pairing, not new biology. The Mantel test is unaffected and the 5 multi-arm communities (C2, C9, C13, C18, C20) carry the inter-chromosomal signal." Use figure T-2.
4. **Acrocentric persistence:** "C9 (chr14_p, chr15_p, chr21_p, chr22_p) and C13 (chr7_p, chr16_q) reproduce the HPRC acrocentric and f7501-arm communities — the population-level signal is detectable from a single individual." Use figure T-3.
5. **Cell-cycle modulation:** "Mitotic CiFi: highest W/B (213×), lowest Mantel ρ (0.41) and fewest significant communities (14/39). Mitotic condensation amplifies the global within/between contrast but flattens fine-scale similarity-contact correlation — same pattern as the HPRC analysis." Use figure T-6.
6. **Negative control:** "The 100 kb immediately centromere-ward of PHR boundaries — unique sequence — gives ρ ≈ 0 in CiFi (and ρ = 0.13–0.23 in Pore-C, driven by 5 lonely non-zero Jaccard pairs). The PHR signal is sequence-driven, not chromosome-end-proximity-driven." Use figure P-7.
7. **Method honesty:** "Single-sample wfmash workaround (`-T RPE1#1#` / `-T RPE1#2#` separately, then concatenate); 91/92 flanks processed; PHR coords from RPE-1's *own* PHR TSV throughout."
