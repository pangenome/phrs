# Provenance register

One entry per artefact: what produced it, from what inputs, with what parameters, and what may and may
not be cited from it. Every note in `response/notes/` opens with a pointer here. Tool versions and pins:
`response/notes/tool_versions_and_pins.md`.

Conventions: paths under `/moosefs/...` are read-only externals; `data/...` are vendored in-repo.
"Verdict-safe" means the artefact may be quoted in the manuscript; "diagnostic only" means it may not.

## P1 — duplicon monophyly (R1-M1, topology)
- script: `scripts/cladistics/duplicon_monophyly.py`; outputs `duplicon_monophyly.{tsv,json}`; note `M1c_duplicon_monophyly.md`
- inputs: `data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz`, `data/chm13-annotations.bed`,
  `data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`, HPRC 100 kb telomeric flank windows
  (`/moosefs/.../flanking/hprcv2.1Mb.flanking_100kb.fa.gz` + `.fai`)
- params: FastGA `-f512` and `-f8192`, min alignment 500 bp; NJ with ≥100 block-bootstrap replicates
- verdict-safe: per-family chromosome-label monophyly and the flank control. **Caveat:** trees resolve
  tandem copy units, not individuals (haplotype-label monophyly 0/20).
- supersedes: `M1_concerted_evolution.md` (distance-based; its negative conclusion is withdrawn for D4Z4,
  OR4F, SEPTIN14 and the acrocentrics — arrays saturate, topology does not)

## P2 — copy enumeration (R1-M1 corroboration)
- script: `scripts/cladistics/impg_copy_enumeration.py`; output `impg_copy_enumeration.tsv`; note `impg_transitivity_copy_enumeration.md`
- inputs: CHM13 annotation/GFF3, communities table, FastGA alignments of 5–20 kb cores
- params: FastGA `-f8192 -c85 -s500 -T32` plus a ≥95% post-filter; own small IMPG index
- **not used:** the shipped 37 GB index `/moosefs/.../hprc25272.bidirectional.impg` — unusable (every
  subcommand demands `--alignment-files`, and the surviving PAF sets are newer than the index, so BGZF
  virtual offsets are stale). Recorded in the note as F8.
- verdict-safe: placements per core and their arm distribution. Caveat: small cores (15 kb) can miss their
  own arm (SEPTIN14 F7); the 95% identity floor is **not** the binding constraint (6/10 cores unchanged at 85%).

## P3 — specificity control (R1-M5 / R3-M3)
- script: `scripts/hic/matched_sd_control.py`; outputs `matched_sd_control.{tsv,json}`; note `M5_matched_sd_control.md`
- inputs: `data/all-vs-all.1Mb.p95.id95.len.tsv`, `data/all-vs-all.1Mb.p95.id95.len.per_window.tsv`,
  Hi-C/Pore-C arm-level contacts (see `end-to-end-report/report/05_hic_validation.md`)
- params: Mantel with average-rank ties; permutation p and bootstrap CI; window classes T (telomeric) vs I (≥1 Mb)
- verdict-safe: T ≈ 0.50 vs I ≈ 0.13, Δ +0.39, p 0.006. **Scope:** this is a WINDOW-CLASS control, not the
  referee's literal matched-segdup control (no SD annotation exists on disk).
- abandoned and not used: wfmash-based SD derivation (out of scope, killed)

## P4 — mouse stage ordering (R1-M7, and corrects M6)
- script: `scripts/hic/mouse_stage_profile_disambiguation.py`; output `mouse_stage_profile_disambiguation.json`; note `M_mouse_stage_ordering.md`
- inputs: `data/mouse_<stage>_phr_{20000,50000}bp_seqlevel.tsv`, `data/mouse_meiosis_sweep/seqlevel/…`,
  `/moosefs/.../mouse_T2T/similarity_1Mb/mouse.dist_matrix.tsv`
- params: dyadic two-way clustered SEs, MRQAP double semi-partialling, arm-block bootstrap 10,000, one seed 20260518
- verdict-safe: per-stage ρ with CIs; zygotene peak; paired zygotene−pachytene contrast −0.0381 [−0.067,−0.009] p 0.017
- **one axis only:** the published odgi matrix equals the arm-pair mean of the per-pair Jaccard column (1.11e-16)

## P5 — mouse satellite fractions (R1-M6), and the M6 retraction
- script: `scripts/hic/mouse_satellite_sensitivity.py`; output `mouse_satellite_sensitivity.json`; note `M6_mouse_satellite.md`
- **RETRACTED:** the masked/stratified axis is void — `mouse_satellite_sensitivity.py:706` stores
  `(QLEN, QSTART)` instead of `(qstart, qend)`, so 4,228/4,228 records are reversed and every derived
  Jaccard is 0. Its stratified numbers may not be cited.
- verdict-safe from this same artefact: per-arm satellite fractions (p-end 0.2597–0.9928, mean 0.8848;
  q-end 0.0021–0.0077, mean 0.0044) and the annotation-vs-proxy status (major `GSAT_MM` real; minor proxied)
- the confound itself is established on the real axis: Corr(satellite, Jaccard) = +0.790 (P4)

## P6 — flank censoring (R1-M8)
- script: `scripts/cladistics/flank_censoring_sensitivity.py`; output `flank_censoring_sensitivity.tsv`; note `M8_flank_censoring.md`
- inputs: `data/all-vs-all.1Mb.p95.id95.len.tsv` (+ `.per_window.tsv`, read-only from /moosefs),
  `data/chm13.phrs.bed`, `data/chm13.chrom.sizes`
- verdict-safe: 1,087 PHRs (6.94%) censored at 500 kb; 1,110 within one window; 43 resumable; caller rule
  replayed exactly on 18,826/18,826 flanks. **Not performed:** 2–5 Mb re-extraction (no longer-flank human
  similarity set exists; source PAFs unrecoverable) — so censoring is bounded, not measured.

## P7 — gene enrichment family dependence (R1-M9)
- script: `scripts/cladistics/enrichment_family_dependence.py`; output `enrichment_family_dependence.tsv`; note `M9_enrichment_dependence.md`
- inputs: `submission/supporting_material/ED_Fig2_arm_functional_components.tsv`,
  `data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz`, `data/chm13.phrs.bed`
- gates: exact hypergeometric tail vs brute force; frozen p-values reproduced with max |Δ| = 0.000 for BH and Holm
- verdict-safe: K = 187; 187 copies from 39 families; 8 families hold 142/187; 124/209 terms single-family;
  spermatid term = 9 loci on 10 chromosomes. **Wording:** the 215 are `LOC*`-without-ontology-source, not "unannotated"

## P8 — pedigree completions (R1-M2)
- scripts: `scripts/pedigree/identity_margin_distribution.py`, `permuted_parent_child_fp.py`;
  outputs `identity_margin_distribution.tsv`, `permuted_parent_child_fp.json`; note `M2_pedigree_completions.md`
- inputs: `data/fig5_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz` + `.fai`,
  assembly `.fai` under `/moosefs/.../washu_v131/assemblies_v11/`
- verdict-safe: 1.02% of windows have margin > 0; 558 windows exceed the retained tract; contig check clean
- note the margin metric is IMPG dice/coverage identity, not base-level % identity

## P9 — PHR-call threshold sensitivity (R3-M6) and the 2,830 pairs (R5-L193)
- script: `scripts/cladistics/phr_call_threshold_sensitivity.py`; output `phr_call_threshold_sensitivity.tsv`; note `R5_pairs_filter_and_R3_thresholds.md`
- inputs: `data/all-vs-all.1Mb.p95.id95.len{,.per_window}.tsv`, `data/chm13.phrs{,.no_acro}.bed`
- gate: replays the published caller exactly (15,668 PHRs, 41 arms, all 18,826 region sets identical)
- verdict-safe: 2,830 = C(77,2) − 96 same-chromosome pairs; arm count invariant to flank size and the
  consecutive-fail rule; moves 41→37 only under the chromosomes-per-window rule

## P10 — superpopulation descriptives (R1-m11)
- script: `scripts/cladistics/phr_by_superpopulation.py`; output `phr_by_superpopulation.tsv`; note `D1_superpopulation_descriptives.md`
- inputs: `data/chm13.phrs.bed`, superpopulation mapping in
  `slides/v2-review-zoom/_revision_assets/v9/labels_superpop_stats_polish/nearest_same_superpop_mds_distances.tsv`
- reproducibility gate: rerun md5-identical, `diff` identical
- verdict-safe: descriptive only; differentiation lives in P3/P11

## P11 — matched F_ST baseline (R1-M6 in the population sense; M12 CIs)
- script: `scripts/popgen/matched_fst_from_blocks.py`; output `scripts/popgen/matched_fst_d_m6.tsv`;
  notes `06_answers_stats_and_bib.md`, `03_change_log.md`
- inputs: `scripts/ci/fst_per_arm_per_pair.tsv` + published baselines (Bhatia 2013; Patterson 2012)
- verdict-safe: ratios 0.72–1.44, all equivalent; block heterogeneity −0.08 to +0.35; r = +0.07 vs sequence sharing
- upstream variant needing moosefs: `matched_fst_d_m6.py` (per-arm contributions)

## P12 — module architecture (R3-M7, ED3)
- scripts: `arm_module_architecture.py`, `arm_module_order.py`, `arm_module_orientation.py`, `peerq3_arm_type_asymmetry.py`;
  outputs `arm_module_{tiling,blocks,inversions,order,orientation}.tsv`, `peerq3_*.tsv`; notes `05_answers_arm_asymmetry_and_architecture.md`
- inputs: HPRCv2 all-vs-all flank alignments (`/moosefs/.../all-vs-all.1Mb.p95/*.paf.gz`),
  `data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
- verdict-safe: 98.6% same telomere-relative orientation (764 pairs, 9.6 Gbp); repertoire Jaccard 0.83 vs 0.57;
  order not conserved; C6 all-Q-Q, C7 all-P-P

## P13 — data integrity (cross-cutting)
- note `L5_data_reconciliation.md`; committed evidence `submission/fig/MainFigures/arm_inclusion_audit.tsv`
- finding: the vendored `data/chm13.phrs.bed` is a 38-arm/15,088-region call, NOT the published 41-arm/
  15,668 set; **`chr13_p` is absent entirely** (a 500 kb published PHR). 500 kb saturation is 8/38 (15.8%),
  not 6/37 (5.4%). Gene copies 402 → 404; ED2's 209 terms unchanged (max |Δp| ≤ 3.7e-7).
- recommendation: regenerate the BED, then re-run `build_inputs` → burden nulls → ED2 as one chain

## P14 — external resources (read-only, no artefact shipped)
- PCLAI local ancestry, S3 `release2_pclai_beds`, three coordinate releases (use `asm_coord`/`chm13_coord`;
  the `grch38_coord` release overstates subtelomeric resolvability). Only HG06807 of the pedigree is present.
- RepeatMasker per haplotype, index `/moosefs/.../repeat_masker_bed_hprc_r2_v1.0.index.csv` (462 haplotypes)
- shipped IMPG index `hprc25272.bidirectional.impg` — **unusable** (see P2)
