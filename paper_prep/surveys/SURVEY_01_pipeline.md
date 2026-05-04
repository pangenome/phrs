# SURVEY 01 — Pipeline section

Source: `end-to-end-report/report/01_pipeline.md` (381 lines).
Scope: contig classification → flank extraction → all-vs-all alignment → inter-chromosomal region detection → pangenome graph + Jaccard similarity → arm-level + sequence-level community detection → cross-arm affinity characterization.

---

## 1. Key findings (with concrete metrics)

1. **Contig classification.** 465 PAF mapping files yielded **12,649 classified contigs** (9,557 pq-spanning, 1,598 q-only, 1,494 p-only) under MIN_LEN ≥ 1 Mb and a strict telomere-presence filter; **12,635 pass full validation**, 14 retained with mixed-strand caveats.
2. **Subtelomeric flank panel.** A 500 kb subtelomeric window (telomere-trimmed) generated **18,827 flanks** across all 48 chromosome arms (24 chrs × 2 arms). Per-arm counts range from **76 (chr13_p) to 458 (chr8_p)**; acrocentric p-arms and sex chromosomes are systematically under-represented because of assembly difficulty.
3. **All-vs-all alignment.** wfmash v0.23.0-41-gb5f0ff1c at `-p 95 -t 48` produced **18,827 PAF.gz files** (each flank as target vs all flanks as queries); the 95% identity floor was chosen to capture the high-identity peak of the bimodal duplicon identity distribution reported by Ambrosini et al. 2007 (peaks at 91% and 98%).
4. **Inter-chromosomal regions.** Sliding-window scan (5 kb window/step, 500 kb max distance, ≥0.95 identity, ≥2 chromosomes, ≥5 alignments/chr, 4 consecutive failing windows, ≥3 kb min region) yielded **18,826 rows**. **15,668 sequences (83.2%) carry an inter-chromosomal match** spanning **41 of 48 arms**; **3,158 (16.8%) are signal-free**. Region lengths: **median 105 kb, mean 144 kb, range 5–500 kb**.
5. **Chimeric contig flagged and removed.** NA18982#1 chr18_q (JBKABS010000018.1, 84.4 Mb) fuses chr18 with **966 kb of chrX PAR1 across a 100 bp NNN scaffold join** at query ~83.37–83.38 Mb; mapq 60 in both wfmash and minimap2 v2.30. A 2,826 bp terminal TTAGGG tract (~471 repeats) precedes the gap; sequence dropped → **15,668 PHRs retained**.
6. **One-copy region validation (Ambrosini 2007).** 4 of 6 historical "one-copy" arms confirmed by zero or PAR-restricted signal: **chr8_q, chr11_q, chr18_q (zero signal); chrX_p/chrY_p with 416/419 sequences (99.3%) matching exclusively chrX/chrY**. The remaining pair (**chr7_q, chr12_q**) is **redefined as a private pair**: **449/449 (100%) of chr12_q and 424/446 (95.1%) of chr7_q** match exclusively chr7/chr12 (chr7_q shared regions median 40 kb, chr12_q median 25 kb). Forms community **C4** at the arm level.
7. **Pangenome graph + similarity matrix.** `pggb -p 95` on 15,668 PHRs + `odgi similarity --all -P` produced a **15,668 × 15,668 Jaccard matrix (~12 GB compressed, ~245 M entries both directions)**, providing the distance measure used downstream.
8. **Arm-level communities.** Leiden (primary): **15 communities, optimal resolution 1.16, silhouette 0.347**. UPGMA (comparison): **14 communities, silhouette 0.342**. **Methods agree exactly on 12 of 15** Leiden communities; the 3 disagreements concentrate around the f7501-carrying cluster (C3) and chr15_q / chr20_p / chr2_q boundaries. Communities recapitulate known biology: D4Z4 (C1: chr4_q/chr10_q), recurrent chr10_p/chr18_p pair (C2), f7501 cluster (C3), chr7_q/chr12_q private pair (C4), RPL23A/WASH/DDX11L module (C5), q-arms incl. acrocentric (C6), acrocentric p-arms (C7), OR4F21 cluster (C11), PAR1 (C15), PAR2 (C14).
9. **f7501 (L78442) population validation reproduces and extends Mefford & Trask 2002.** Per-arm distribution across **465 haplotypes** (233 samples × 2; minus 1 for CHM13). Fixed sites confirmed: **chr3_q 91.8%, chr19_p 90.5%, chr15_q 85.6%**. AFR enrichment confirmed: **chr16_q OR=17.4, p=6.6e-27**; chr7_p OR=3.7, p=8.2e-03. Three new AFR-enriched arms: **chr8_p (90% AFR, p=8.5e-05), chr16_p (76% AFR, p=6.7e-07), chr9_q (43% AFR, p=1.9e-05)**. **chr2_q is SAS-enriched (77% SAS, p=6.8e-11, OR=23.0)**; **chr6_p is AMR-enriched (62% AMR, p=7.0e-04, OR=7.4)**. **chr15_q EUR enrichment p=2.5e-04** (98.5% EUR vs 64.9% AFR carriers). Three FISH-detected sites (chr4_q, chr19_q, plus partial chr1_p / chr20_p / chrX_q at very low frequency) are below sequence-alignment detection threshold or appear as novel low-frequency loci.
10. **Sequence-level communities.** Leiden on the full 15,668 × 15,668 graph at **k-NN = 75, resolution = 0.8 → 50 communities, modularity 0.97**. **18 pure (1 arm), 23 near-pure (>90%), 9 mixed**. Most polymorphic arms by # communities: **chr6_q (8), chr19_p (7), chr3_q (6), chr7_q/chr11_p/chr5_q/chr16_q/chr20_p (5 each)**. Major mixed communities include **C4 acrocentric-p (770 seqs)**, **C3 D4Z4 (712 seqs, chr4_q ≈ chr10_q ≈ 50/50)**, **C32 PAR2 (432)**, **C33 PAR1 (416)**, **C40 f7501 multi-arm (352)**, **C13 chr17_p/chr11_p (508)**, **C26 chr7_p/chr11_p (490)**, **C27 chr6_q/chr15_q (298)**.
11. **Arm × sequence nesting.** Most monolithic arm communities (≥97% in one seq community): **C1/D4Z4 99.4%, C7/acro-p 99.6%, C13/chr4_p 99.8%, C14/PAR2 99.8%, C15/PAR1 99.3%, C10/chr17_p 97.1%**. Most fragmented: **C3/f7501 → 16 seq communities (largest 18.0%)**, **C11/OR4F21 → 15 (24.7%)**, **C6/q-arms → 10 (16.7%)**, **C5/RPL23A-WASH → 8 (26.3%)**. Partition agreement Arm vs Seq: **ARI = 0.35, NMI = 0.76**. Silhouette: **arm 0.347 → sequence 0.602**.
12. **Cross-arm affinity (sequence absorption into foreign communities).** **1,740 / 15,668 sequences (11.1%) are cross-arm**. **7 arms 100% absorbed**: chrY_p → chrX_p, chrY_q → chrX_q (PAR1 / PAR2), chr15_p / chr21_p / chr13_p / chr22_p → chr14_p (rDNA-adjacent acrocentric homogenization), chr10_q → chr4_q (D4Z4). Most polymorphic autosomal arm: **chr11_p (62.7% cross-arm; 281/448, distributed across chr9_q, chr17_p, chr7_p, chr5_q communities)**. **8 arms with 0% cross-arm** (truly arm-private content): chr10_p, chr12_p, chr17_q, chr18_p, chr1_q, chr20_q, chr21_q, chr9_q.

---

## 2. Existing figures/plots referenced

The pipeline section text **does not embed figure paths** (no `.png` / `.pdf` / `.svg` references). Plot-producing scripts referenced in the section:

| Script (path as listed in source) | Likely figure outputs |
|---|---|
| `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R` | Distance heatmap, MDS, UMAP, UPGMA dendrogram on the 41×41 arm matrix |
| `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-seq-community-structure.R` | Sequence-level (50-community) UMAP / community heatmaps |
| `/moosefs/guarracino/HPRCv2/scripts/community/detect_communities.R` | Leiden modularity / silhouette scans (arm-level k-scan, sequence-level k×resolution heatmap) |

Comparative external figures invoked (not files in this repo):
- **Mefford & Trask 2002, Fig. 3** — f7501 FISH per-arm distribution (the basis of the 13-arm population-enrichment table at lines 165–179).
- **Linardopoulou et al. 2005, Fig. 5** — recurrent chr10_p/chr18_p transfer pair (C2) and OR4F21 block (C11, C5 boundary).
- **Ambrosini et al. 2007** — bimodal duplicon identity distribution (peaks 91% / 98%) and the 6 one-copy regions list.

> Action item for survey: file paths to actual rendered PDFs/PNGs are NOT present in `01_pipeline.md`. They likely exist under `/moosefs/guarracino/HPRCv2/.../similarity/` or a `figures/` sibling — these need to be discovered by checking the script outputs (out-of-scope for this survey but flagged in §5).

---

## 3. Existing data tables / CSVs referenced (paths verbatim)

**Input data**
- `/moosefs/pangenomes/HPRCv2/*.fa.gz` — 465 HPRCv2 assemblies.
- `/moosefs/guarracino/HPRCv2/PHR_III/pq-classification/chm13.centromeres.approximate.bed` — CHM13 centromere coordinates.

**Intermediate**
- `/moosefs/guarracino/HPRCv2/PHR_III/pq-classification/contig_classifications.tsv` — 12,649 classified contigs.
- `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz` — 18,827 subtelomeric flanks.
- `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95/*.paf.gz` — 18,827 all-vs-all PAFs (~88 GB).
- `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv` — inter-chromosomal region calls (one row per flank).
- `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz` — 15,668 PHR sequences.
- `/moosefs/guarracino/HPRCv2/PHR_III/pggb/.../similarity.tsv.gz` — Jaccard pairwise similarity (~10.8 GB compressed).

**Output**
- `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.dist_matrix.rds` — 15,668 × 15,668 distance matrix (R RDS).
- `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv` — arm-level (15 community) assignments.
- `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv` — sequence-level (50 community) assignments.
- `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv` — 41 × 41 arm distance matrix.

**Scripts** (full list reproduced in source lines 342–355)
- `classify_contigs.py`, `trim-telomeres.sh`, `phr_wfmash_array.sh`, `phr_post_wfmash.sh`, `find-multichr-regions-incremental.py`, `extract_flanking_sequences.py`, `plot-similarity-subtelo.R`, `community/detect_communities.R`, `community/extract-seq-assignments.R`, `community/community-utils.R`, `community/compare-community-levels.R`, `similarity/plot-seq-community-structure.R`.

---

## 4. Statistical / computational methods

- **Contig classification**: deterministic rule-based classifier from CHM13 PAF alignments (chromosome-alias-matched only); telomere-presence filter and per-arm telomere-count constraints.
- **Flank extraction**: deterministic 500 kb window inward from telomere; telomere repeat trimming.
- **All-vs-all alignment**: `wfmash v0.23.0-41-gb5f0ff1c -p 95 -t 48 --quiet`; one PAF per flank (each flank as target).
- **Inter-chromosomal region scan**: incremental sliding window via impg index; thresholds (5 kb / 5 kb / 500 kb / id ≥ 0.95 / ≥2 chrs / ≥5 alns/chr / 4 fail windows / ≥3 kb output).
- **Chimera flagging**: cross-aligner agreement (wfmash + minimap2 v2.30 -x asm20, mapq 60) plus Flagger annotation of the NNN gap.
- **Pangenome graph**: `pggb -p 95 -D /scratch` on 15,668 PHRs.
- **Pairwise similarity**: `odgi similarity --all -P` → Jaccard over shared graph nodes.
- **Arm-level Leiden**: distance = 1 − Jaccard; arm-pair distance = mean over sequence pairs; weighted graph with `w_ij = exp(-d_ij / median(d))`; resolution scan 0.1–3.0 step 0.01; selected by max mean silhouette.
- **UPGMA comparison**: `hclust(method = "average")` on the 41 × 41 arm distance matrix; cut at k ∈ [2, 20] maximizing silhouette.
- **Sequence-level Leiden**: k-NN graph (k ∈ {10, 25, 50, 75, 100, 125}) × resolution (0.1–3.0 step 0.1); modularity-maximizing pair within target community count 5–50 selected (chosen point: k=75, res=0.8, modularity 0.97).
- **Partition comparison**: silhouette score (−1 to 1; arm 0.347, seq 0.602); Adjusted Rand Index (ARI = 0.35); Normalized Mutual Information (NMI = 0.76).
- **f7501 carrier detection**: `minimap2 -x asm20`, one-vs-one alignment of L78442.1 (36.3 kb) against each flank; threshold ≥ 30 kb matching bases per haplotype.
- **Population enrichment**: per-arm, per-superpopulation **one-sided Fisher's exact (greater)** with sizes AFR=134, AMR=88, EAS=104, EUR=65, SAS=74; reported odds ratio + p-value for the most-enriched superpopulation.
- **Cross-arm affinity classification**: assignment of each sequence's seq-community dominant arm vs its own arm at the k50 partition.

No multiple-testing correction is reported in the source for the f7501 enrichment table — flagged below.

---

## 5. Open gaps for a Nature-style manuscript

1. **Multiple-testing correction.** The f7501 per-arm × per-superpopulation Fisher's exact tests (~13 arms × 5 populations ≈ 65 tests in the displayed subset) lack BH-FDR or Bonferroni correction. For Nature, all p-values must be q-values or carry an explicit "uncorrected" label and a corrected companion column.
2. **Resolution-stability evidence.** The Leiden choices (arm: r=1.16, k15; seq: k=75, r=0.8) are reported as silhouette- / modularity-maximizing single points. A stability assessment is needed: (a) silhouette/modularity surface across the (k, r) grid; (b) bootstrap or subsampling stability of community assignments; (c) sensitivity of `1,740` cross-arm count to (k, r).
3. **UPGMA vs Leiden disagreement at boundaries.** The 3-community disagreement around C3/C8/C12 (f7501 cluster + chr2_q / chr15_q / chr20_p) is described qualitatively. A confusion matrix with per-arm probabilistic membership across both methods would harden the boundary claims.
4. **Background / null model.** No explicit null is presented for "x% cross-arm affinity" or for the 50-community modularity (0.97). Compare against permuted graphs and shuffled arm labels.
5. **Sensitivity to identity threshold.** The 95% identity choice excludes the 91% Ambrosini peak. A sensitivity analysis at p ∈ {90, 92.5, 95, 97.5} on a sub-panel (e.g., one chromosome family) would show how communities depend on the threshold.
6. **Chimera screen completeness.** Only one chimera (NA18982#1 chr18_q) is reported. A systematic NNN-gap screen across all 15,669 retained sequences should be documented even if it returns zero further hits.
7. **Acrocentric / sex-chromosome under-representation.** Per-arm counts show chr13_p (76), chrY_p (92), chr21_p (116), chr15_p (119), chr22_p (223), chrY_q (101), chr14_p (230), chrX_p (327), chrX_q (330) substantially below the autosome modal value (~440–460). A rendered figure showing assembly completeness alongside community participation is needed.
8. **Effect-size confidence intervals.** OR values up to 23.8 are reported as point estimates without CIs; reviewers will want bootstrap or exact-method 95% CIs.
9. **Figure assets do not exist as paths in the source.** Plot scripts are listed but the rendered figures (heatmaps, UMAPs, network diagrams) are not pointed to. The manuscript will need a figure-asset inventory.
10. **Chimera detection metadata.** The Flagger labels (NNN, Hap) are referenced but not tabulated; for the manuscript a small supplementary file with the chimera coordinates and labels is needed.
11. **Method version pinning beyond wfmash.** wfmash and minimap2 are pinned; pggb, odgi, Leiden (igraph?), R packages are not. A `versions.yml` or methods table is needed.
12. **Reproducibility / code release plan.** Scripts are referenced by absolute filesystem path on `/moosefs/...`. A repo snapshot with relative paths is needed.

---

## 6. Suggested figures (main + extended/SI)

Status legend: **EXISTS** = a script in §3 plausibly already emits this; **NEW** = needs generation.

### Main figures

| ID | Status | One-line caption idea |
|----|--------|---|
| F1 | NEW | Pipeline overview schematic — 465 assemblies → 18,827 flanks → 15,668 PHRs → 15 arm + 50 sequence communities (single panel, with counts on each arrow). |
| F2 | EXISTS | 41 × 41 arm-level distance heatmap with Leiden community blocks annotated (C1–C15) and UPGMA dendrogram on top — driven by `arm_dist_matrix.tsv`. |
| F3 | EXISTS | Sequence-level UMAP / force-directed layout colored by 50-community partition; key communities (C3 D4Z4, C4 acrocentric-p, C32 PAR2, C33 PAR1, C40 f7501) labelled — from `plot-seq-community-structure.R`. |
| F4 | NEW | f7501 carrier-frequency × superpopulation matrix — heatmap of % carriers per (arm × pop), starred cells = significant after BH-FDR; Mefford 2002 Fig. 3 reference row alongside. |
| F5 | NEW | Cross-arm affinity diagram — circular plot of 41 arms with edges weighted by # sequences absorbed across arm boundaries; emphasizes the chr11_p hub and the 7 fully-absorbed arms. |

### Extended / Supplementary

| ID | Status | One-line caption idea |
|----|--------|---|
| EF1 | NEW | Per-arm flank counts (bar chart, 48 arms) with assembly QC overlay; highlights acrocentric / sex chromosome under-representation. |
| EF2 | NEW | Inter-chromosomal region length distribution (histogram + ECDF) on 15,668 PHRs; median 105 kb / mean 144 kb annotated. |
| EF3 | NEW | Chimera evidence panel for NA18982#1 chr18_q — wfmash + minimap2 alignment dotplot, NNN gap, Flagger track. |
| EF4 | NEW | Leiden resolution scan: silhouette vs resolution for the arm-level run (0.1–3.0); peak at r=1.16, silhouette 0.347. |
| EF5 | NEW | Sequence-level (k × resolution) modularity heatmap; selected point (k=75, r=0.8, mod=0.97) starred. |
| EF6 | EXISTS | UPGMA dendrogram with cut at k=14; communities colored to compare with Leiden 15-community labels. |
| EF7 | NEW | Confusion matrix Arm-Leiden vs Sequence-Leiden membership (15 × 50) with ARI=0.35, NMI=0.76 in figure title. |
| EF8 | NEW | f7501 hits per haplotype across superpopulations — ridge plot or violin per arm, annotated with OR and p (corrected). |
| EF9 | NEW | One-copy region validation panel — 6 historical regions × 4 status outcomes (zero-signal / PAR / private-pair / undetected). |
| EF10 | NEW | Pipeline parameter sensitivity: identity ∈ {90, 92.5, 95, 97.5} → community count and ARI vs the p=95 reference. |
| EF11 | NEW | Pangenome graph stats: node count, edge count, average node degree, component sizes. |

---

## 7. Suggested talk slides (15-min talk)

| # | Slide title | One-line takeaway |
|---|---|---|
| 1 | Why subtelomeres? | Subtelomeres are recombination/exchange hotspots whose population-scale architecture has been invisible until pangenome assemblies. |
| 2 | The dataset | 233 HPRCv2 samples → 465 haplotypes → **18,827 telomere-anchored 500 kb flanks** across all 48 arms. |
| 3 | Pipeline at a glance | Classify → flank → all-vs-all (95% id) → inter-chr region → graph → Jaccard → Leiden — five stages, all reproducible. |
| 4 | What is a PHR? | **15,668 / 18,827 (83.2%)** flanks share sequence with another chromosome; 16.8% are arm-private. |
| 5 | One-copy regions revisited | Confirms 4 of Ambrosini's 6 historical one-copy regions; redefines **chr7_q/chr12_q as a private pair**. |
| 6 | Arm-level communities | **15 Leiden communities** (silhouette 0.347) recapitulate D4Z4, PAR1/2, acrocentric-p, f7501 cluster — biology-faithful. |
| 7 | Method robustness | UPGMA vs Leiden agree on **12 of 15** arm-level communities; differences localize to f7501 boundaries. |
| 8 | f7501 lives — but reshaped | We reproduce Mefford & Trask 2002 at 465 haplotypes and add **3 new AFR-enriched arms** + chr2_q SAS + chr6_p AMR. |
| 9 | f7501 chr15_q EUR shift | **98.5% EUR vs 64.9% AFR (p=2.5e-04)** — ongoing African-population loss or incomplete lineage sorting. |
| 10 | Sequence-level zoom | At k=75 resolution=0.8 we obtain **50 communities (modularity 0.97)** — chr6_q is the most polymorphic arm (8 communities). |
| 11 | Cross-arm absorption | **11.1% of sequences cluster with a foreign arm**; chrY → chrX absorbed entirely; chr10_q → chr4_q (D4Z4); 4 acrocentric p-arms → chr14_p. |
| 12 | chr11_p — the polymorphic outlier | **62.7% cross-arm**; sequences distribute across chr9_q, chr17_p, chr7_p, chr5_q communities. |
| 13 | Truly arm-private content | **8 arms with 0% cross-arm** define the lower bound of subtelomeric exchange. |
| 14 | What this lets us validate | Output feeds 3D-genome (Hi-C / Dip-C) validation, gene/repeat enrichment, mouse / RPE1 comparators. |
| 15 | Summary + next steps | Pipeline reproducible end-to-end; manuscript-ready figures need FDR correction + stability bootstraps + chimera screen. |
