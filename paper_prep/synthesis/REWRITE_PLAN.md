---
title: "REWRITE PLAN — canonical materials → Nature companion + BoG-this-week"
author: "agent-747 (audit-canonical-materials)"
date: 2026-05-05
anchor: paper_prep/synthesis/ABSTRACT.md
companion_to: AUDIT_REPORT.md
---

# REWRITE PLAN — canonical materials → Nature companion + BoG-this-week

## Constraints

- **Target venue.** *Nature*, as a **companion paper to the HPRC v2 main paper**. Pre-arranged for parallel review alongside HPRC v2. The companion-paper framing is non-negotiable for the introduction (claim **C1**).
- **Submission horizon.** ~2–3 weeks from today (2026-05-05).
- **Near-term event.** Lead author Erik Garrison gives a talk at **Biology of Genomes (BoG) THIS WEEK** based on this work. Figure improvements and synthesis material that strengthen the talk are **HIGH PRIORITY** and live in **Lane A**, separate from the longer manuscript lane.

### Nature-format constraints to respect
(per current Nature *Article* guidelines — confirm against `nature.com/nature/for-authors/formatting-guide` at submission time):
- **Abstract:** ≤ 200 words. Current `ABSTRACT.md` is ~360 words → needs trimming for submission, but the canonical claims must remain intact. Trimming is part of Lane B.
- **Main text:** ~3,000 words (Article format). Current `MANUSCRIPT_DRAFT.md` is ~3,500–4,000 words and built around an off-spec wider thesis → needs full rewrite.
- **Display items:** typically 4–6 main figures + tables; up to 10 Extended Data figures.
- **Methods:** structured, at the end (not interleaved).
- **References:** typically ≤ 50 for Articles.
- **Reflows / formats:** display items are formatted Nature-style (panel labels, captions ≤ 200 words). The **rewrite plan does not need to reformat existing material in this task** — Lane B tasks flag where current materials violate these constraints so future tasks address them.

### How tasks are organised
- **Lane A (BoG-this-week):** highest priority. Deliverables within 1–4 days. Figure improvements and synthesis content that materially strengthen the talk. 5–10 tasks.
- **Lane B (Manuscript-for-Nature):** priority, 2–3 week horizon. Methods writing, results sections (one per major claim), intro with HPRC v2 companion framing, discussion with concerted-evolution thesis, references, Nature-format assembly, render-with-figures-verified.

### Total task count
**Total: 33 tasks** (Lane A: 8 + Lane B: 25). Within the 25–60 cap mandated by the task brief.

### Task entry template
Every task uses this template:

```
### TASK-NN: <imperative title>
Lane: A or B
Inputs: <specific files / data sources>
Output: <specific filename + format>
Acceptance: <one checkable condition; for any render task, MUST include
            'figures visible inline at expected positions, verified by image-object count >= N
             or page screenshots'>
Depends on: <prior TASK-NN ids, or 'none'>
```

---

## Lane A: BoG-this-week (highest priority, deliverable within days)

### TASK-01: Build NJ tree from arm-level Jaccard distance and annotate the abstract's named clades
Lane: A
Inputs: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv` (41 × 41 arm Jaccard distance); R `ape::nj()`; `ABSTRACT.md` clade list.
Output: `paper_prep/figures/fig1/figure_fig1_nj_panel.pdf` + `figure_fig1_nj_panel.png` + `nj_tree.newick` + `nj_clade_annotations.tsv`.
Acceptance: NJ tree renders 41 arm tips; the named clades from ABSTRACT.md (Xp/Yp via PAR1, Xq/Yq via PAR2, acrocentric short arms, 10p–18p, {22q,21q,19q,1q,13q,17q}, 4q–10q DUX4) are each visibly annotated and labelled on the rendered PNG. Bootstrap support (≥ 100 replicates) reported on each named-clade node.
Depends on: none.

### TASK-02: Add explicit clade labels to existing Fig 1c arm-distance heatmap for talk
Lane: A
Inputs: `paper_prep/figures/fig1/figure_fig1.{R,pdf,png}` (current Leiden+UPGMA arm heatmap); TASK-01 NJ assignments.
Output: `paper_prep/figures/fig1/figure_fig1_v2.pdf` + `figure_fig1_v2.png` (annotated overlay with abstract's clade names).
Acceptance: Annotated PNG visibly labels each abstract-named clade (10p–18p block; {22q,21q,19q,1q,13q,17q} block; 4q–10q DUX4 block; PAR1 and PAR2 blocks; acrocentric block) directly on the arm-distance heatmap. PNG file ≥ 200 KB and visibly carries all annotations when opened.
Depends on: TASK-01.

### TASK-03: Produce PAR2-comparable subtelomere identity calibration panel for headline 1a
Lane: A
Inputs: `chm13.phrs.bed`, `CHM13-HG002.sub-telo-phrs.bed`, `hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz`, `plot-impg-coverage.R` (top-level repo).
Output: `paper_prep/figures/fig1/figure_fig1_par2_inset.pdf` + `.png` + `par2_vs_subtelo_lengths.tsv`.
Acceptance: One panel showing PAR2 length distribution (~334 kb anchor) overlaid against per-arm subtelomeric homology block length distribution (median 105 kb, mean 144 kb from `all-vs-all.p95.id95.len.tsv`). Text annotation: "subtelomeric homology comparable in scale to PAR2." PNG visibly carries both distributions and the annotation.
Depends on: none.

### TASK-04: Re-cut TALK_OUTLINE for BoG (15-min slot, anchored on canonical 8 claims C1-C8)
Lane: A
Inputs: `paper_prep/synthesis/TALK_OUTLINE_15MIN.md`, `paper_prep/synthesis/ABSTRACT.md`, AUDIT_REPORT.md §1 figures audit.
Output: `paper_prep/synthesis/TALK_OUTLINE_BOG.md` (replaces or supplements existing TALK_OUTLINE_15MIN.md).
Acceptance: Slide list explicitly maps slides → C1..C8 (cell per slide). Pedigree, RPE-1, mouse, sperm, Dip-C content reduced to ≤ 1 slide total (single "validation" slide if retained at all). HPRC v2 companion framing on slide 1 or 2.
Depends on: TASK-01, TASK-02, TASK-03 (so the talk references the new figures).

### TASK-05: BoG dry-run rendered-slides PDF
Lane: A
Inputs: TASK-04 outline; figures from `paper_prep/figures/fig1/` and `ed5/` and `ed8/`; `figure_fig1_nj_panel.pdf` from TASK-01; `figure_fig1_par2_inset.pdf` from TASK-03.
Output: `paper_prep/synthesis/talk_BoG_2026.pdf` (slide deck).
Acceptance: PDF renders ≥ 14 slides; **every slide that references a figure has the figure visible inline (verified by `pdfimages -list talk_BoG_2026.pdf | wc -l` ≥ N where N = number of figure-bearing slides)**; per-slide PNG snapshots in `talk_BoG_2026_pages/` directory.
Depends on: TASK-04.

### TASK-06: Annotate ED2c chord plot with abstract clade names for talk
Lane: A
Inputs: `paper_prep/figures/ed2/figure_ed2.{R,pdf,png}` (current 50-community chord); TASK-01 NJ clades.
Output: `paper_prep/figures/ed2/figure_ed2c_clade_annotated.pdf` + `.png`.
Acceptance: Chord plot has the named edges (PAR1, PAR2, acrocentric→14p, 10q→4q D4Z4, 10p–18p) explicitly text-labelled on the rendered PNG.
Depends on: TASK-01.

### TASK-07: Produce a single one-pager handout summarising C1-C8 for in-person BoG follow-up conversations
Lane: A
Inputs: `ABSTRACT.md`, AUDIT_REPORT.md, fig1/2/3 captions.
Output: `paper_prep/synthesis/BoG_handout.pdf` + `BoG_handout.md`.
Acceptance: 1-page handout (US Letter or A4); each of C1-C8 named with one-sentence finding + one figure reference + one supporting number. PDF embeds the headline genome-wide identity heatmap as a thumbnail (image-object count ≥ 1, verified).
Depends on: TASK-03 (PAR2 calibration available for the C4 row of the handout).

### TASK-08: Pre-talk regression check — run all Lane A figure scripts cleanly under pinned environment
Lane: A
Inputs: `paper_prep/synthesis/VERSIONS.md`, all Lane A figure scripts.
Output: `paper_prep/synthesis/lane_A_render_log.md` (per-task: command, exit code, output paths, timing).
Acceptance: All Lane A figure scripts (TASK-01, -02, -03, -06) exit 0 in a clean env reproducible from `VERSIONS.md`. Log records every input file's size + checksum.
Depends on: TASK-01, TASK-02, TASK-03, TASK-06.

---

## Lane B: Manuscript-for-Nature (priority, 2–3 week horizon)

### TASK-09: Establish HPRC v2 companion-paper framing in introduction (C1)
Lane: B
Inputs: `ABSTRACT.md`, `framing_synthesis.md` (top-level repo) §1 "three-paper arc", HPRC v2 main-paper draft if available, `MANUSCRIPT_DRAFT.md` §1 introduction.
Output: `paper_prep/synthesis/sections/01_introduction.md`.
Acceptance: Introduction is ≤ 600 words; **first paragraph explicitly names this paper as the subtelomeric-companion-to-HPRC-v2 publication**; Guarracino 2023 (acrocentric) and de Lima/Guarracino 2025 (Robertsonian) are cited as predecessors; HPRC v2 main paper is cited and the companion-paper relationship is stated explicitly. References section auto-builds the new HPRC v2 entry.
Depends on: none.

### TASK-10: Write implicit-pangenome-graph methods section (C2)
Lane: B
Inputs: `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95/` (PAF directory), `/moosefs/guarracino/HPRCv2/scripts/phr_wfmash_array.sh`, `find-multichr-regions-incremental.py`, `wfmash.sh`, `extract_telo_flanks.mouse.sh`; `ABSTRACT.md` paragraph 1.
Output: `paper_prep/synthesis/sections/06_methods_implicit_pangenome.md` (sub-section of Methods).
Acceptance: Section explains (a) reference-free all-to-all design (no chromosomal partitioning); (b) **the ~12 % pairwise sampling figure — what fraction of possible pairs were evaluated, computed from the actual PAF set, and why** (this number must be derived from the on-disk PAFs, not asserted); (c) wfmash parameters (asm20, p95, id95, len ≥ 30 kb); (d) why this is "implicit pangenome graph" rather than a built pggb graph (i.e., the analysis is on the PAF set itself rather than on a constructed GFA). One paragraph each (≤ 200 words each).
Depends on: none.

### TASK-11: Compute and tabulate the ~12 % pairwise sampling fraction from existing PAFs
Lane: B
Inputs: `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv` and PAF list at `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.paf.list`.
Output: `paper_prep/synthesis/pairwise_sampling.tsv` + a one-paragraph note `paper_prep/synthesis/pairwise_sampling_note.md`.
Acceptance: TSV has columns (`total_possible_pairs`, `evaluated_pairs`, `fraction`, `derivation`). Fraction is computed as `evaluated_pairs / total_possible_pairs`. The note explains how the value relates to the abstract's "~12 %".
Depends on: none.

### TASK-12: Resolve 465 vs 466 haplotype headcount and produce dataset table (C3)
Lane: B
Inputs: `/moosefs/guarracino/HPRCv2/assemblies/`, `/moosefs/guarracino/HPRCv2/PHR_III/HiC/sample_info.tsv`, HPRC v2 main-paper sample table.
Output: `paper_prep/synthesis/dataset_table.tsv` + `paper_prep/synthesis/sections/06_methods_dataset.md`.
Acceptance: TSV has one row per sample × haplotype with included/excluded flag + reason; the total = 466 and the 466 number is reconciled with HPRC v2's published count. Methods sub-section cites the canonical 466.
Depends on: none.

### TASK-13: Run NJ on sequence-level distance matrix + bootstrap; export .newick
Lane: B
Inputs: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.dist_matrix.tsv` (15,668 × 15,668 sequence Jaccard) and `hprcv2.1Mb.subtelo.arm_dist_matrix.tsv` (41 × 41 arm Jaccard).
Output: `paper_prep/synthesis/nj_outputs/{arm,seq}.nj.newick`, `{arm,seq}.nj.bootstrap.tsv`.
Acceptance: Both .newick files load with `ape::read.tree()` without error; bootstrap support (≥ 100 replicates) computed for every internal node; all named clades from `ABSTRACT.md` resolve as monophyletic with bootstrap support reported.
Depends on: none (Lane A TASK-01 is a quick first cut; this is the manuscript-quality version with the full sequence-level tree as well).

### TASK-14: Write Results sub-section on the genome-wide identity survey (C4)
Lane: B
Inputs: `paper_prep/figures/fig1/` (panels 1a/1b), TASK-03 PAR2-calibration TSV, `chm13.phrs.bed`.
Output: `paper_prep/synthesis/sections/02_results_identity_survey.md`.
Acceptance: ≤ 500 words; cites the 18,827 telomere-anchored 500 kb flanks (or, if the abstract is rewritten to drop "500 kb", the new value); cites the median 105 kb / mean 144 kb PHR length; explicitly states the "comparable to PAR2" claim with the PAR2 calibration number from TASK-03; references Fig 1a/1b inline.
Depends on: TASK-03.

### TASK-15: Write Results sub-section on the NJ-tree cladistic analysis (C5)
Lane: B
Inputs: TASK-13 NJ outputs, TASK-01 annotated NJ figure, `cross_arm_affinity_sequences.tsv`, `gene_copy_summary.csv` (DUX4 row from `paper_prep/_brainstorming/`).
Output: `paper_prep/synthesis/sections/03_results_cladistics.md`.
Acceptance: ≤ 700 words; explicitly names every clade from the abstract (Xp/Yp via PAR1, Xq/Yq via PAR2, acrocentric short arms, 10p–18p, the {22q,21q,19q,1q,13q,17q} clade, 4q–10q DUX4 with copy-number diversity, large moderate-similarity clade) with bootstrap support per clade; references the NJ-tree figure inline; cites DUX4 copy-number diversity from the salvaged copy-number table.
Depends on: TASK-01, TASK-13.

### TASK-16: Run PCA on similarity matrix and write Results sub-section on PCA + community detection (C6)
Lane: B
Inputs: `hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`, `hprcv2.1Mb.subtelo.dist_matrix.tsv`, `cross_arm_superpop_enrichment.tsv`, `fst_superpop_matrix.tsv`, sample × superpopulation labels.
Output: `paper_prep/synthesis/pca_outputs/pca_arm.tsv` + `pca_seq.tsv` + `pca_seq.png`; `paper_prep/synthesis/sections/04_results_pca_communities.md`.
Acceptance: PCA computed on a feature matrix (haplotype × PHR-presence-binary or haplotype × arm-membership-vector); the first 3 PCs are reported with variance-explained; superpopulation hulls overlaid; the section text references the PCA figure inline and the community-detection partition (Leiden / NJ-derived) inline.
Depends on: TASK-13.

### TASK-17: Write Results sub-section on Hi-C 3D nuclear-envelope-proximity hypothesis (C7)
Lane: B
Inputs: `paper_prep/figures/fig3/` panel 3a (HG002 Pore-C contact matrix), `paper_prep/figures/ed5/` (multi-resolution + exclusion controls), `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/` and `analysis/human/exclusion_controls/`; LAD or Lamin B1 ChIP-seq overlay if available.
Output: `paper_prep/synthesis/sections/05_results_hic_envelope.md`.
Acceptance: ≤ 600 words; cites HG002 within-vs-between contact ratio and per-resolution stability; **explicitly addresses the "nuclear envelope" wording from the abstract** — either by citing the LAD overlay (preferred) or by acknowledging the gap and proposing the test as a follow-up; references ED5b (exclusion controls) as evidence the signal is not driven by nucleolar / PAR contacts alone.
Depends on: TASK-22 (LAD overlay).

### TASK-18: Trim canonical ABSTRACT.md to ≤ 200 words for Nature submission while preserving C1-C8
Lane: B
Inputs: `ABSTRACT.md`.
Output: `paper_prep/synthesis/abstract_submission.md` (do **not** overwrite `ABSTRACT.md`; treat the original as the long-form anchor).
Acceptance: Word count ≤ 200; every C1-C8 claim still present (verified by checklist in the file's frontmatter); HPRC v2 companion framing in the first sentence.
Depends on: none.

### TASK-19: Write Discussion section anchored on the concerted-evolution-and-unorthodox-recombination thesis (C8)
Lane: B
Inputs: `paper_prep/figures/ed8/` (panels a, b: feedback loop + D4Z4-CTCF-lamin), TASK-15 cladistics text, TASK-17 Hi-C text.
Output: `paper_prep/synthesis/sections/07_discussion.md`.
Acceptance: ≤ 500 words; explicitly closes the abstract's title-thesis ("ongoing recombination shapes subtelomeres → concerted evolution and unorthodox recombination"); cites the named NJ-tree clades; cites the within-cluster identity peak; cites the Hi-C envelope-proximity test (or honestly flags the missing LAD overlay); references ED8a/b inline.
Depends on: TASK-15, TASK-17.

### TASK-20: Build canonical references.bib for the rewritten manuscript
Lane: B
Inputs: existing `paper_prep/synthesis/REFERENCES.bib`, plus HPRC v2 main paper, NJ-tree methods (Saitou & Nei 1987 or modern equivalent), and any new C2/C5 citations.
Output: `paper_prep/synthesis/REFERENCES_v2.bib` (do not overwrite REFERENCES.bib until rewrite freezes).
Acceptance: ≤ 50 entries (Nature Article cap); HPRC v2 main paper cited; Saitou & Nei 1987 cited for NJ method; every citation key used in any `paper_prep/synthesis/sections/*.md` file resolves to an entry in this .bib.
Depends on: TASK-09, TASK-10, TASK-12, TASK-13, TASK-14, TASK-15, TASK-16, TASK-17, TASK-19.

### TASK-21: Update fig1 multi-panel composite for manuscript (C4 + C5 + C6)
Lane: B
Inputs: existing `paper_prep/figures/fig1/figure_fig1.{R,pdf,png}`, TASK-01 NJ panel, TASK-03 PAR2 calibration, TASK-16 PCA outputs.
Output: `paper_prep/figures/fig1/figure_fig1_v3.{R,pdf,png}` + updated `caption.md` + updated `sources.tsv`.
Acceptance: 4-panel composite: (a) genome-wide identity heatmap with PAR2 calibration overlay; (b) NJ tree with annotated clades; (c) arm-level Jaccard heatmap with NJ-derived ordering and named-clade boxes; (d) PCA scatter with superpopulation hulls. Caption ≤ 200 words. Composite PNG ≥ 500 KB and visibly carries all 4 panels when opened.
Depends on: TASK-01, TASK-03, TASK-16.

### TASK-22: Compute and overlay LAD / Lamin B1 ChIP track on subtelomeric Hi-C contacts (C7 envelope test)
Lane: B
Inputs: HG002 Hi-C mcool (`/moosefs/guarracino/HPRCv2/PHR_III/HiC/HG002/hicpro_output/cool/HG002.mcool`); public Lamin B1 ChIP-seq or LAD calls for a matched cell type (lookup external resource or flag as missing); `chm13.phrs.bed` for subtelomeric coordinates.
Output: `paper_prep/figures/fig3/figure_fig3_envelope.{R,pdf,png}` + `lad_overlap_stats.tsv`.
Acceptance: One-panel figure showing subtelomeric Hi-C contacts overlaid with LAD coverage at chromosome ends; per-arm enrichment statistic. **If no matching LAD dataset is available, the task fails honestly and the abstract is amended to soften the "nuclear envelope" claim — that decision is logged in the task output.**
Depends on: none.

### TASK-23: Update fig3 to focus on Hi-C only (drop Dip-C/sperm/mouse/CiFi from main fig)
Lane: B
Inputs: existing `paper_prep/figures/fig3/figure_fig3.{R,pdf,png}` and `caption.md`; TASK-22 envelope panel.
Output: `paper_prep/figures/fig3/figure_fig3_v2.{R,pdf,png}` + updated caption.
Acceptance: 4-panel main fig: (a) HG002 Hi-C contact matrix ordered by NJ-tree partition; (b) within-vs-between summary across HPRC Hi-C samples (Hi-C only — drop Pore-C, CiFi, Dip-C, sperm, mouse from main); (c) ED5b-style exclusion-control summary; (d) LAD-overlap envelope panel from TASK-22 (or honest gap flag). Composite PNG visibly carries all 4 panels.
Depends on: TASK-22.

### TASK-24: Demote off-spec figures (Fig 4, ED3, ED4) to SI or scrap
Lane: B
Inputs: `paper_prep/figures/fig4/`, `paper_prep/figures/ed3/`, `paper_prep/figures/ed4/`; AUDIT_REPORT.md §1 verdicts.
Output: `paper_prep/figures/_si/` directory containing the demoted figures + a `_si/README.md` explaining placement; updated `MANUSCRIPT_SKELETON_v2.md`.
Acceptance: Files moved with `git mv`; the new `_si/README.md` lists each demoted figure and the section-level reason from the audit; the new skeleton no longer cites these as main / ED slots.
Depends on: none.

### TASK-25: Restructure Methods to be end-of-paper (Nature format) and structured into named sub-sections
Lane: B
Inputs: TASK-10, TASK-12, TASK-13, TASK-16, TASK-17, TASK-22 outputs.
Output: `paper_prep/synthesis/sections/06_methods.md` (assembled from the per-task method stubs).
Acceptance: Methods is its own .md; sub-sections in this order: (i) Dataset (TASK-12); (ii) Implicit pangenome graph (TASK-10); (iii) Pairwise sampling (TASK-11); (iv) NJ tree construction (TASK-13); (v) PCA + community detection (TASK-16); (vi) Hi-C 3D analysis (TASK-17); (vii) Statistics. Each sub-section ≤ 250 words. References resolve to `REFERENCES_v2.bib`.
Depends on: TASK-10, TASK-11, TASK-12, TASK-13, TASK-16, TASK-17.

### TASK-26: Update CAPTIONS.md to match the rewritten figure set
Lane: B
Inputs: TASK-21 (new fig1), TASK-23 (new fig3), TASK-24 (demotions); existing `CAPTIONS.md`.
Output: `paper_prep/synthesis/CAPTIONS_v2.md` (do not overwrite until rewrite freezes).
Acceptance: One caption per surviving display item (≤ 200 words each); statistical conventions carried over from `CAPTIONS.md`; every demoted figure's caption is moved to the SI section of the file.
Depends on: TASK-21, TASK-23, TASK-24.

### TASK-27: Re-run STATS_AUDIT.md against the new figure / claim set
Lane: B
Inputs: existing `STATS_AUDIT.md` framework, TASK-15 / TASK-16 / TASK-17 / TASK-19 sections.
Output: `paper_prep/synthesis/STATS_AUDIT_v2.md` + per-test corrections under `stats_audit_v2/`.
Acceptance: Every p-value, OR, ρ in the rewritten sections has either a BH-FDR q-value within its declared family or an explicit `(uncorrected; single combined test)` annotation; every OR has an exact / conditional-MLE 95 % CI; family definitions stated.
Depends on: TASK-15, TASK-16, TASK-17, TASK-19.

### TASK-28: Assemble manuscript draft v2 (single .md) from all sections
Lane: B
Inputs: TASK-09, TASK-14, TASK-15, TASK-16, TASK-17, TASK-19, TASK-25, TASK-18 (abstract); TASK-20 (refs).
Output: `paper_prep/synthesis/MANUSCRIPT_DRAFT_v2.md`.
Acceptance: Single .md file in this order — front matter / abstract (≤ 200 words from TASK-18) / introduction / results 1-N / discussion / methods / references. Word count of main text (excluding methods, captions, refs) is 2,500–3,500. Every figure / table referenced ≥ once. Every citation resolves.
Depends on: TASK-09, TASK-14, TASK-15, TASK-16, TASK-17, TASK-18, TASK-19, TASK-20, TASK-25, TASK-26.

### TASK-29: Render manuscript draft v2 PDF with mandatory figure-presence verification
Lane: B
Inputs: `MANUSCRIPT_DRAFT_v2.md`, all figure PDFs / PNGs in `paper_prep/figures/`, `REFERENCES_v2.bib`.
Output: `paper_prep/synthesis/MANUSCRIPT_DRAFT_v2.pdf` + per-page screenshots in `paper_prep/synthesis/MANUSCRIPT_DRAFT_v2_pages/`.
Acceptance: **Figures visible inline at expected positions, verified by `pdfimages -list MANUSCRIPT_DRAFT_v2.pdf | wc -l` ≥ N where N = (number of main figures + number of ED figures + number of in-text image objects)**, and per-page PNG screenshots show inline figures at the expected positions for at least the headline figure (Fig 1) and Hi-C figure (Fig 3). Failure of either check fails the task. PDF page count ≥ expected for a Nature Article.
Depends on: TASK-21, TASK-23, TASK-28.

### TASK-30: Pre-render figure-file existence pre-check task
Lane: B
Inputs: `MANUSCRIPT_DRAFT_v2.md` (or its preceding draft).
Output: `paper_prep/synthesis/figure_existence_check.tsv` + go/no-go decision in the task log.
Acceptance: For every `paper_prep/figures/<id>/figure_<id>.{pdf,png}` referenced in the manuscript .md, the file exists on disk. If any is missing, the task fails (TASK-29 cannot run).
Depends on: TASK-28.

### TASK-31: Build extended-data figures index for v2 (now without ED3/ED4/ED6/ED7)
Lane: B
Inputs: surviving ED figures (`ed1`, `ed2`, `ed5`, `ed8`); TASK-24 demotions.
Output: `paper_prep/synthesis/extended_data_index.md` + per-ED `caption.md` consolidated into the rewritten CAPTIONS_v2.md.
Acceptance: ≤ 10 ED figures (Nature Article cap); each ED figure has a corresponding source script and caption; the index lists ED1, ED2, ED5, ED8 + any new ED slot the rewrite produces.
Depends on: TASK-24.

### TASK-32: Trim REFERENCES_v2.bib to ≤ 50 entries (Nature Article cap)
Lane: B
Inputs: `REFERENCES_v2.bib`, every citation key actually used in `MANUSCRIPT_DRAFT_v2.md`.
Output: `paper_prep/synthesis/REFERENCES_final.bib`.
Acceptance: Entry count ≤ 50; every entry is cited at least once in `MANUSCRIPT_DRAFT_v2.md`; HPRC v2 main-paper companion citation present.
Depends on: TASK-28.

### TASK-33: Final pre-submission checklist + Nature-format compliance scan
Lane: B
Inputs: `MANUSCRIPT_DRAFT_v2.{md,pdf}`, `REFERENCES_final.bib`, `CAPTIONS_v2.md`, `STATS_AUDIT_v2.md`, `extended_data_index.md`.
Output: `paper_prep/synthesis/SUBMISSION_CHECKLIST.md`.
Acceptance: Every Nature constraint (abstract ≤ 200 words; main text 2,500–3,500 words; ≤ 6 main figures; ≤ 10 ED figures; ≤ 50 references; structured Methods at end; every figure caption ≤ 200 words) is verified with a checkbox + the actual measured value. Process-failure note (figure-presence verification) is checked.
Depends on: TASK-29, TASK-31, TASK-32.

---

## Total task count

- **Lane A (BoG-this-week):** 8 tasks (TASK-01 through TASK-08)
- **Lane B (Manuscript-for-Nature):** 25 tasks (TASK-09 through TASK-33)
- **Total:** 33 tasks. **Within the 25–60 cap.**

## Specific tasks that MUST appear (cross-check against the task brief)

The brief mandates these specific tasks. Each is present:

- ✓ A task that establishes HPRC v2 companion-paper framing in the intro → **TASK-09**.
- ✓ A task (or tasks) that write the implicit-pangenome-graph methods section (~12 % pairwise sampling explanation) → **TASK-10** (methods section) + **TASK-11** (pairwise-sampling computation).
- ✓ One task per abstract claim C4–C7 producing the corresponding results subsection:
  - C4 → **TASK-14** (genome-wide identity survey).
  - C5 → **TASK-15** (NJ-tree cladistic analysis).
  - C6 → **TASK-16** (PCA + community detection).
  - C7 → **TASK-17** (Hi-C 3D + envelope hypothesis).
- ✓ A figure-render task that includes inline-figure verification → **TASK-29** (and the pre-check at TASK-30).

## Lane order-of-operations (informal)

- **Day 1–2 (Lane A):** TASK-01 (NJ tree first cut), TASK-03 (PAR2 calibration), TASK-04 (talk outline) → quick wins for BoG.
- **Day 2–3 (Lane A):** TASK-02 (annotated fig1), TASK-06 (annotated chord), TASK-05 (talk slides), TASK-07 (handout), TASK-08 (regression check).
- **Week 1 (Lane B):** TASK-09–TASK-13 (intro, methods, dataset, NJ at scale).
- **Week 2 (Lane B):** TASK-14–TASK-19 (results sections, abstract trim), TASK-21–TASK-24 (figure rebuild + demotion).
- **Week 2–3 (Lane B):** TASK-22 (LAD overlay — risky / external-resource dependent), TASK-25–TASK-27 (methods assembly, captions, stats), TASK-28–TASK-33 (assemble + render-with-verification + format compliance).

---

*End of REWRITE_PLAN.md.*
