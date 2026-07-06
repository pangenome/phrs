---
title: "AUDIT REPORT — canonical materials vs ABSTRACT.md"
author: "agent-747 (audit-canonical-materials)"
date: 2026-05-05
anchor: paper_prep/synthesis/ABSTRACT.md
target_venue: Nature (companion to HPRC v2 main paper)
submission_horizon: ~2-3 weeks
near_term_event: Biology of Genomes (BoG) talk this week — Erik Garrison lead-presenting
---

# AUDIT REPORT — canonical materials vs ABSTRACT.md

**Anchor.** `paper_prep/synthesis/ABSTRACT.md` — "Concerted evolution and unorthodox recombination of human subtelomeres" (Guarracino & Garrison). Eight load-bearing claims, labelled C1–C8 in the task brief:

- **C1.** HPRC v2 companion-paper framing (Nature, alongside HPRC v2 main).
- **C2.** Implicit pangenome graph: reference-free, all-to-all, ~12 % pairwise sampling, no chromosomal partitioning.
- **C3.** 465 near-complete assemblies from HPRC v2.
- **C4.** Genome-wide identity survey: extended (10s–100s kb) interchromosomal homology at nearly all subtelomeres, comparable to PAR2.
- **C5.** NJ-tree cladistic analysis: Xp/Yp + Xq/Yq via PARs, acrocentrics; 10p–18p; the {22q,21q,19q,1q,13q,17q} clade; 4q–10q DUX4 with copy-number diversity; large moderate-similarity clade.
- **C6.** PCA + community detection on similarity matrix → subtelomere clustering across human populations.
- **C7.** Hi-C 3D maps testing nuclear-envelope-proximity recombination hypothesis.
- **C8.** Synthesizing thesis: ongoing recombination shapes subtelomeres → "concerted evolution and unorthodox recombination".

This audit follows the brief: figures-first (Erik needs them for BoG this week), then synthesis docs, then data/code, then brainstorming, then a process-failure note.

---

## 1. Figures audit (PRIORITY — Erik needs this for BoG this week)

The figure inventory under `paper_prep/figures/` contains **10 directories**: `fig1` … `fig4` (main) and `ed1, ed2, ed3, ed4, ed5, ed8` (Extended Data). **`ed6` and `ed7` are absent on disk despite being specified in `MANUSCRIPT_SKELETON.md`** — flagged in section 5.

Each directory has been inspected: rendered PDF + PNG, R or Python source script, `caption.md`, `sources.tsv` (ed2 also has `extract_within_community_jaccard.sh`; fig1 also has `architecture_per_arm.tsv`). PNG file sizes range 151 KB–735 KB, all visibly carry the panels described in their captions.

Mapping uses ABSTRACT.md claim labels: **C1** (HPRC-v2 companion framing — narrative only, not figure-mappable), **C2** (implicit pangenome graph), **C3** (465 near-complete assemblies — narrative), **C4** (genome-wide identity survey), **C5** (NJ-tree cladistic clades), **C6** (PCA + community detection), **C7** (Hi-C 3D), **C8** (synthesis — discussion display item).

| Figure | What it actually shows | Maps to abstract claim | Aligned? | Action | Reasoning | BoG-talk value |
|---|---|---|---|---|---|---|
| **fig1** | (a) genome-wide stacked identity heatmap, 465 near-complete assemblies × 24 chr, 100 kb windows (with chr18q inset); (b) genome-wide n-chromosomes-sharing heatmap with PHR BED overlay; (c) 41 × 41 arm-level Jaccard distance heatmap with **Leiden k=15** community blocks and **UPGMA k=14 dendrogram on top** (14/15 agreement); (d) per-arm architecture-category bar (homogeneous 4 / polymorphic 28 / fully interchangeable 9). | C4 (1a, 1b directly), C5 (1c is the Leiden/UPGMA cladistic visual — proxy for the abstract's NJ tree but **not the same algorithm**), C6 (1c is the partitioned similarity matrix; PCA itself not shown) | **partial** | **MINOR-REVISE** | Panels 1a/1b are the canonical genome-wide identity survey for C4 — keep as-is. Panel 1c is the closest existing analogue to C5 (cladistic structure across arms) but uses Leiden + UPGMA rather than the abstract's NJ tree, and the named clades from the abstract (10p–18p; {22q,21q,19q,1q,13q,17q}; 4q–10q DUX4) are **not annotated on the dendrogram** — needs annotation overlay or a re-run as NJ. Panel 1d is correct architecture taxonomy but is built around the "communities" framing of the off-target draft, not the abstract's cladistic framing. | **high** — 1a/1b are the headline visual; 1c needs clade-name annotation before the talk |
| **fig2** | (a) allele-vs-paralog Wilcoxon paired test across 9 multi-arm communities (C7 = acrocentric reversed); (b) two-domain model: per-arm Spearman ρ + piecewise breakpoint (39/48 arms); (c) cross-arm × superpopulation Fisher heatmap + Hudson Fst matrix; (d) UPGMA "out-of-Africa" tree from cross-arm Fst. | C5 partial (1d is a cross-arm Fst tree, not the within-similarity NJ tree the abstract names); C6 (cross-arm × population enrichment is a population-detection product on the similarity matrix — partially aligned). **Allele-vs-paralog and two-domain are NOT in the abstract.** | **partial** (allele/paralog and two-domain are off-spec; population panels are partially on-spec) | **REDO** (a, b panels) / **MINOR-REVISE** (c, d) | The abstract has no allele-vs-paralog or Flint–Mefford two-domain framing. Panels 2a and 2b carry the prior-session "communities" thesis. Panels 2c/2d are population-structure work that is closer to C6 but is built on the cross-arm Fst representation, not on PCA. **For the BoG talk, the population panels (2c/2d) carry value; 2a/2b should be cut from the talk and reworked or moved to ED for the manuscript.** | **medium** — keep 2c/2d for talk; drop 2a/2b for talk |
| **fig3** | (a) HG002 Pore-C inter-arm contact matrix, 50 kb, ordered by sequence community (B/W = 0.056, p = 3.9e-85); (b) convergent-evidence forest plot, 14 tests across 6 technologies (Hi-C, Pore-C, CiFi, Dip-C, sperm scHi-C, mouse meiotic Hi-C); (c) S_all negative-control box plots; (d) flanking paradox (PHR vs flanking 100 kb centromere-ward). | C7 (Hi-C 3D maps — yes, but the abstract specifies *Hi-C*, not the wider single-cell/sperm/mouse-meiotic suite shown in 3b) | **partial** | **MINOR-REVISE for talk; REDO panel-set selection for manuscript** | Panel 3a is on-spec for C7 — the canonical headline 3D visual. Panels 3b/3c bring in Dip-C, sperm scHi-C, mouse meiotic Hi-C, CiFi, Pore-C — these are *not* in the canonical abstract, which specifies only Hi-C 3D maps. For the Nature companion, the 3D figure should narrow to bulk Hi-C across HG002 + a small set of HPRC samples, with single-cell / mouse / sperm material moved to ED or cut. Panel 3d (flanking paradox) is a good methodological control regardless. | **high** for 3a/3d; medium for 3b (replace with Hi-C-only forest); cut 3c for talk |
| **fig4** | (a) WashU 3-gen T2T pedigree untangle ribbons; (b) CEPH1463 cross-assembler 11-feature parent matrix; (c) RPE-1 t(X;10) rediscovery; (d) mouse zygotene per-PHR-pair Jaccard vs Hi-C contact (ρ = 0.715, p = 4.4e-55, n = 344). | **None of C1–C8.** Pedigrees, RPE-1 self-validation, and mouse cross-species are entirely outside the canonical abstract. | **N** | **SCRAP** (for the canonical Nature companion) | The abstract makes no claim about pedigrees, RPE-1 single-individual validation, or mouse meiosis. Fig 4 is the most off-target main figure. Panels here may be retained as **supplementary** (e.g., a single SI figure showing pedigree as evidence of ongoing exchange — supports C8 thesis indirectly) but should not occupy a main-text figure slot. The mouse panel (4d) is methodologically nice but cross-species generalisation is not in the abstract. | **low** — only 4a (pedigree as direct evidence of ongoing recombination) has talk value; consider one slide |
| **ed1** | (a) pipeline schematic (465 → 18,827 flanks → 18,827 PAFs → 15,668 PHRs → 15/50 communities); (b) per-arm flank counts (48 arms); (c) PHR length distribution (median 105 kb, mean 144 kb); (d) chr18q chimera (NA18982#1) evidence panel. | C2 partial (pipeline schematic is the closest current artifact to "implicit pangenome graph methods" but does **not** explain the ~12 % pairwise-sampling design, all-to-all reference-free aspect, or the no-chromosomal-partitioning principle); C3 partial (cites 465, abstract says 465). | **partial** | **REDO** | The pipeline schematic must be re-cast as the **implicit pangenome graph schematic** with the ~12 % pairwise-sampling design, reference-free framing, and the all-to-all-without-chromosomal-partitioning principle made explicit. Current schematic reads as a wfmash/impg/pggb pipeline diagram — correct but missing the conceptual frame. Update assembly count to 465 near-complete assemblies (232 HPRC v2 individuals plus CHM13, per abstract). | **high** — Erik will be asked "what is the implicit pangenome graph?" — needs the explicit schematic |
| **ed2** | (a) UMAP of 15,668 PHRs by 50-community Leiden; (b) within-community Jaccard distance bimodality across 8 communities; (c) cross-arm affinity radial chord plot (top edges: chrY_q→chrX_q PAR2 329, chrY_p→chrX_p PAR1 287, chr15_p→chr14_p 280, chr10_q→chr4_q D4Z4 220); (d) confusion matrix arm-Leiden × seq-Leiden (ARI 0.35, NMI 0.76). | C5 partial (the chord plot in 2c surfaces the PAR1, PAR2, acrocentric, and 4q–10q DUX4 edges named in the abstract); C6 partial (UMAP is a non-linear projection, not PCA proper). | **partial** | **MINOR-REVISE** | Panel 2c chord plot is a strong asset for C5 — visually shows the cladistic edges the abstract names (PAR1/2, acrocentric, D4Z4). Should be promoted toward a main figure or made one of the headline ED panels. Panel 2a UMAP and 2d confusion matrix are off-spec methodology surrogates for the abstract's PCA + community detection. | **medium** — 2c is a great talk asset; 2a/2b/2d less so |
| **ed3** | (a) TAR1 prevalence per arm (94.6 % overall); (b) (TTAGGG)n island length distribution and motif composition; (c) terminal telomere length by community; (d) per-arm TAR1 distance-from-telomere. | **None of C1–C8.** | **N** | **SCRAP** for the canonical paper | Annotation panels (TAR1, ITS, telomere length) are interesting but not part of the canonical abstract. They belong in a separate annotation-focused paper or in supplementary if retained at all. | **none** for Nature companion; medium for the talk if Erik wants to mention annotation context |
| **ed4** | (a) GO:BP top terms (PHR-only, 23 protein-coding genes); (b) copy-weighted vs deduplicated GO enrichment; (c) high-copy gene families (top 15 — DUX4, OR4F, IL9R, FRG2, etc.); (d) OR4F pseudogenisation gradient by arm (11.1 % chr7p → 99.8 % chr15q). | **None of C1–C8** directly; panel (c) names DUX4 and OR4F that surface in the C5 4q–10q DUX4 clade context, but enrichment-/pseudogene-framing is not in the abstract. | **N** | **SCRAP** for canonical | Pure gene-enrichment / pseudogene material is the off-target ORA / OR4F / DUX4 framing that previous-session synthesis hill-climbed on (now in `paper_prep/_brainstorming/`). Not in canonical abstract. The 4q–10q DUX4 copy-number diversity in the abstract is supported by raw copy-number tables, not by these enrichment panels. | **none** |
| **ed5** | (a) W/B contact across 5 mcool resolutions × 8 datasets; (b) Mantel ρ before/after acrocentric+sex exclusion (50 kb); (c) O/E-normalised within-vs-between contact; (d) per-community reproducibility heatmap (15 communities × 11 datasets). | C7 (Hi-C robustness panels — directly support the Hi-C 3D claim with multi-resolution + confound-exclusion stability). | **Y** | **KEEP** (with one revision: drop the per-community columns reliant on 50-community sequence-Leiden if the manuscript switches off Leiden in favour of NJ + PCA partition) | This is the strongest existing ED for the canonical claim C7. The exclusion controls (panel b) are the right rigorous-defense panel for nuclear-envelope-proximity claims because they prove the signal is not driven by nucleolar (acrocentric) or PAR contacts alone. | **high** — Erik should show panel (b) as the "we tested for confounds" slide |
| **ed8** | (a) causal feedback loop schematic (PHR self-reinforcement); (b) D4Z4-CTCF-lamin tethering for C1 (chr4q ↔ chr10q); (c) recombination rate vs cross-arm affinity null (ρ = -0.35 → ρ = -0.01 after callability filter); (d) compartment identity at tips (HG002, A vs B compartment). | C7 partial (panel b mentions lamin tethering — directly aligned with the nuclear-envelope-proximity hypothesis); C8 partial (panel a is the synthesizing-feedback-loop schematic). | **partial** | **MINOR-REVISE** | Panels (a) and (b) are the most direct visual support for C8 (synthesizing thesis) and the nuclear-envelope mechanism for C7. They should be retained and possibly moved into the main Discussion as the closing display item. Panel (c) is good honest-null work but is not in the canonical abstract; consider keeping in ED as a negative-control panel. Panel (d) is supplementary. | **high** for (a) and (b); medium for (c) |

### Figures missing for the abstract's claims

The canonical ABSTRACT.md commits to specific visualisations and analyses that have **no current figure**:

1. **NJ-tree cladistic figure (C5).** The abstract names "neighbor-joining trees of subtelomeric similarity" with specific clades — Xp/Yp via PAR1, Xq/Yq via PAR2, acrocentric short arms, 10p–18p, the {22q,21q,19q,1q,13q,17q} tight clade, 4q–10q DUX4 with copy-number diversity, and a large moderate-similarity clade. **No NJ tree currently exists** in `paper_prep/figures/` or in `/moosefs/guarracino/HPRCv2/PHR_III/similarity/` (UPGMA k=14 and Leiden k=15 partitions exist; UMAP and MDS exist; no NJ output). This is the **single largest figure gap**.
2. **PCA on the similarity matrix (C6).** The abstract specifies "Principal component … analyses of the similarity matrix … further resolve subtelomeric clustering across human populations." The repository has **MDS** (`hprcv2.1Mb.subtelo.full_mds.rds`, MDS scatter PNGs by community / superpopulation) and **UMAP** but no explicit PCA. MDS on a Jaccard distance matrix is metric-equivalent to PCoA, not PCA, and the abstract specifically says PCA. Either re-run as PCA on a feature matrix, or replace abstract wording with "MDS / PCoA". (Reframing the abstract is **out-of-scope** for this audit — the rewrite plan flags a small task to either run PCA or re-word the abstract claim.)
3. **Identity heatmap "comparable to PAR2" comparison (C4).** The abstract claims subtelomeric homology blocks span "tens to hundreds of kilobases at nearly all subtelomeres — comparable in scale to canonical pseudohomologous systems such as PAR2 on the sex chromosomes." Fig 1a shows the genome-wide identity heatmap but does **not include a PAR2 reference panel** for visual side-by-side calibration. A small comparative panel (e.g., PAR2 identity track overlaid on a representative subtelomere of similar size) would directly support the headline analogy.
4. **Hi-C nuclear-envelope-proximity test (C7).** The abstract phrases the Hi-C role as "evaluat[ing] the hypothesis [that PHR exchange is] facilitated by the physical proximity of subtelomeres at the nuclear envelope." Existing Fig 3 + ED5 + ED8b mostly show within-community-vs-between contact enrichment. A **direct nuclear-envelope-association panel** — e.g., subtelomeric contact-with-LADs (lamin-associated domains) or radial-position from a Lamin B1 ChIP / DamID overlay — is **not present**. ED8b cites lamin-tethering for C1 but does not provide the genome-wide envelope-proximity test the abstract names.
5. **465-assembly headcount panel (C3).** Current ED1a/b cite 465. Either the abstract is wrong by one (CHM13 reference may or may not count toward the 465) or the dataset description needs to be updated. A 1-line headcount panel and an explicit Methods sentence would resolve this.

### Figures that exist but are off-spec for the abstract (candidates to demote / cut)

- **fig4** (pedigree + RPE-1 + mouse) — none in abstract; demote to SI or cut.
- **ed3** (TAR1 / ITS / telomere length) — none in abstract; demote to SI or cut.
- **ed4** (gene enrichment / pseudogene gradient) — off-target; this is the ORA / OR4F material that was correctly flagged in the prior `park-off-target` task, but the gene-enrichment ED4 *figure itself* is still in the canonical figure tree. Demote to SI or cut.
- **fig2 panels (a) and (b)** (allele-vs-paralog, two-domain Flint–Mefford) — off-spec methodology surrogates for the cladistic and population claims.

---

## 2. Synthesis docs audit

For each `.md` and supporting file under `paper_prep/synthesis/` (excluding the anchor `ABSTRACT.md`):

| Filename | Topic actually covered | Aligned with abstract? | Salvageable content (specific section refs) | Reason for verdict |
|---|---|---|---|---|
| `MANUSCRIPT_SKELETON.md` | Nature-format skeleton (~3000 words) built around "Population-scale subtelomeric communities mirror 3D nuclear organisation across human and mouse" thesis. 4 main figs + 8 ED. Headline numbers: 232 individuals, 465 near-complete assemblies, 18,827 flanks, 15,668 PHRs, 41 arms → 15 / 50 communities. Includes Pore-C, Dip-C, sperm, mouse meiotic, pedigree integration. | **partial** | Headline numbers (232 / 465 / 18,827 / 15,668), pipeline-stage list (§Headline numbers and §Methods), Fig 1a/1b/1c entries, Fig 3a (HG002 Pore-C contact matrix), and ED1, ED5, ED8 entries. | Skeleton is anchored on a wider thesis than the canonical abstract: it integrates Dip-C / sperm scHi-C / mouse meiotic / pedigree material that is **not** in ABSTRACT.md. Sections covering implicit-pangenome-graph methods (C2), NJ tree (C5), PCA (C6) are **absent** in the skeleton; sections covering pedigree / mouse / sperm scHi-C are **present but off-spec**. Needs a major prune + insertion of C2/C5/C6 sections. |
| `MANUSCRIPT_DRAFT.md` (~409 lines) | Compiled rough draft of the Nature manuscript per the wider thesis. Abstract block, 7 main sections (Intro, Communities, Heterogeneity, 3D, Pedigree+Cross-species, Discussion, Methods). Includes 24-entry bibliography. | **partial** | Intro paragraph #2 (the "three classes of question newly addressable") is reusable; Methods paragraph 1 ("Assemblies and flank extraction") and paragraph 4 ("3D analysis"); Discussion paragraph 1 ("we have shown…"). | Built around "15 communities" framework, not the canonical NJ-tree + PCA + community detection + nuclear-envelope-proximity-Hi-C framework. Pedigree / mouse / sperm sections are off-spec. The abstract block in this draft does not match `ABSTRACT.md` — needs replacement with the canonical anchor abstract verbatim. |
| `MANUSCRIPT_DRAFT.typ` | Typst version of the draft, mirrors `.md`. | **partial** | Same as `.md`. | Same misalignment as `.md`. The render output is the off-target figureless PDF that triggered the process-failure note (§5). |
| `CAPTIONS.md` | Deduplicated captions for Fig 1–4 + ED1–8 with statistical conventions. | **partial** | Fig 1 caption (matches the figure on disk and is on-spec for C4/C5), Fig 3a caption (matches, on-spec for C7), ED5 caption (on-spec for C7), ED8a/b captions (on-spec for C7/C8). | Fig 4 / ED3 / ED4 / ED6 / ED7 captions are off-spec (Fig 4 = pedigree+mouse; ED3 = TAR1/ITS; ED4 = enrichment; ED6 = Dip-C+sperm+RPE-1; ED7 = mouse). Some captions reference figures (ED6, ED7) that **do not exist on disk**. |
| `STATS_AUDIT.md` (268 lines) | BH-FDR + 95 % CI audit on every p-value in the headline-numbers block. Family definitions, q-values, conditional-MLE OR for f7501 chr16q. | **partial** | Family-definition framework, BH-FDR conventions, 95 % CI for any odds ratio cited in the canonical paper; the specific Mantel multi-resolution table (40 tests) — directly applicable to C7 robustness claims. | Most p-values audited belong to off-spec analyses (allele-vs-paralog, sperm W/B, mouse zygotene). The framework and conventions are universal and salvageable; the per-test corrections need to be re-scoped after the rewrite picks the canonical p-value family. |
| `SCRIPT_INVENTORY.md` (301 lines) | Catalogue of every script cited in any survey. Lives at `/moosefs/guarracino/HPRCv2/scripts/`. | **Y** (universal) | Entire file. Especially §1 (top-level pipeline), §2 (community detection), §3 (similarity), §4 (Hi-C). | Path-only inventory; venue-agnostic. Required for Methods regardless of which thesis the paper picks. |
| `VERSIONS.md` (96 lines) | Pinned tool versions (wfmash, pggb, odgi, igraph, R, Python). | **Y** | Entire file. | Universal. Required for Nature Methods. |
| `WORK_DECOMPOSITION.md` (191 lines) | Phase-3 figure tasks + Phase-4 validation tasks for the prior-session manuscript build. | **N** | Decomposition pattern (fan-out-merge with file-scope rule) is a useful template. Specific task IDs are obsolete. | Built around the off-spec figure set (includes ed6/ed7 that were never produced). Replaced by the new REWRITE_PLAN.md. |
| `ARCHITECT_TASK_BRIEF.md` (122 lines) | Phase-2 architect brief for the prior-session synthesis. Describes the pipeline that produced the off-spec MANUSCRIPT_SKELETON.md. | **N** | Process notes on file-scope rule and integrator pattern are reusable. | Anchored on a thesis (15 communities, 3D + pedigree + mouse) that is wider than the canonical abstract. Superseded by the audit-canonical-materials brief. |
| `ACCEPTANCE_CHECKLIST.md` (105 lines) | "12 anchoring findings" cross-reference for the off-spec skeleton; 27-entry novel-contributions ledger. | **partial** | Findings 1, 3, 8, 9, 12 (community structure, cross-arm exchange rate, 3D mirroring sequence communities, flanking paradox, two-domain) overlap with C4 / C5 / C7. | The "12 anchoring findings" framing is the off-spec thesis. Replaced by the canonical 8-claim (C1–C8) framing in this audit. |
| `TALK_OUTLINE_15MIN.md` (146 lines) | 15-minute conference talk outline. Same off-spec wider thesis. | **partial** | Slides 1 (cytogenetic-era continuity) and 4 (41 → 15 communities, named clades) are reusable for the BoG talk. Slide 11 (pedigrees) and slide 14 (acks) need replacement. | The talk outline is currently calibrated for a venue different from BoG — too much pedigree/mouse content. Needs a focused BoG re-cut. |
| `NOVEL_CONTRIBUTIONS.tsv` (28 lines, 27 + header) | 27-row novel-contributions ledger. | **partial** | Rows that map to C4 / C5 / C7 are reusable. | Built around 27 contributions of the wider thesis; many (sperm, mouse, pedigree, RPE-1) are off-spec. |
| `LIMITATIONS_X_FINDINGS.tsv` (13 lines) | 12 anchoring findings × applicable limitations × bound-on-interpretation. | **partial** | Limitation entries (95 % identity threshold, 500 kb flank, multi-mapping at PHRs) are universal. | Findings axis is the off-spec set; limitation axis is universal. |
| `REFERENCES.bib` (415 lines, 24 entries) | Bibliography for the prior-session draft. | **Y** (mostly) | All entries: MeffordTrask2002, Flint1997, Linardopoulou2005, Ambrosini2007, Riethman2004, Stong2014, Bailey2002, Trask1998, vanDeutekom1996, Stout1999, Rouyer1986, Lemmers2010, Masny2004, Ottaviani2009, plus modern entries. | Will need to add HPRC v2 main paper (companion citation for C1), the implicit-pangenome-graph methodology citation (C2), and any NJ tree methods citation (C5). |
| `pandoc_convert.log`, `render.log` | Render diagnostics for the off-spec MANUSCRIPT_DRAFT.{md,typ}. | n/a | n/a | Build artifacts; no content audit needed. |
| `stats_audit/` (subdir) | Per-test TSVs (f7501 fisher, mantel multires) backing STATS_AUDIT.md. | **Y** (universal data tables) | All TSVs. | Directly reusable for any p-value report regardless of thesis. |

### Synthesis-doc summary

- **Reusable as-is (Y):** SCRIPT_INVENTORY.md, VERSIONS.md, REFERENCES.bib (with additions), stats_audit/ TSVs.
- **Partially reusable (partial):** MANUSCRIPT_SKELETON.md, MANUSCRIPT_DRAFT.md/typ, CAPTIONS.md, STATS_AUDIT.md, ACCEPTANCE_CHECKLIST.md, TALK_OUTLINE_15MIN.md, NOVEL_CONTRIBUTIONS.tsv, LIMITATIONS_X_FINDINGS.tsv.
- **Off-target (N):** WORK_DECOMPOSITION.md, ARCHITECT_TASK_BRIEF.md.

The dominant failure pattern is that the prior synthesis session anchored on a *different thesis* than `ABSTRACT.md` — it integrated single-cell 3D, mouse meiosis, pedigrees, and out-of-Africa, none of which are in the canonical abstract.

---

## 3. Data & code assets audit

For each abstract claim C1–C8, what's on disk that supports it.

### C1. HPRC v2 companion-paper framing
- **Code:** none directly applicable (this is a framing claim).
- **Data:** none directly applicable.
- **Reference:** `framing_synthesis.md` at repo root contains a long discussion of "the third in a series" (Guarracino 2023 acrocentric → de Lima/Guarracino 2025 Robertsonian → this paper) — directly supports the companion-paper arc but does not name HPRC v2 explicitly enough.
- **MISSING — needs to be produced:** explicit companion-paper anchoring (HPRC v2 citation in REFERENCES.bib; intro paragraph framing this paper as the subtelomeric companion to HPRC v2 main).

### C2. Implicit pangenome graph (~12 % pairwise sampling, no chromosomal partitioning)
- **Pipeline scripts:** `/moosefs/guarracino/HPRCv2/scripts/phr_wfmash_array.sh`, `phr_post_wfmash.sh`, `find-multichr-regions-incremental.py`, `wfmash.sh`, `extract_telo_flanks.mouse.sh`, `partitionPGGB.sh`.
- **Data:** `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv` (the all-vs-all PAF aggregate), `all-vs-all.1Mb.p95.id95.len.tsv`. PAF files at `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95/` (18,827 PAF.gz).
- **MISSING — needs to be produced:** the explicit "~12 % pairwise sampling" calculation. The all-vs-all PAFs exist, but the *sampling fraction* — what fraction of the 18,827 × 18,827 / 2 possible pairs were actually evaluated and why — needs to be (a) computed from the existing PAFs and (b) explained as a methodological choice in a Methods sub-section. The "no chromosomal partitioning" principle exists by construction (`partitionPGGB.sh` is *not* used for this paper's main partition) but needs to be stated explicitly. **No schematic figure of the implicit pangenome graph methodology currently exists** (cf. ED1 audit row above).

### C3. 465 near-complete assemblies from HPRC v2
- **Data:** `/moosefs/guarracino/HPRCv2/assemblies/` (HPRC v2 v1.1 assemblies: two haplotype assemblies for each of the 232 individuals plus the haploid CHM13 anchor = 465 near-complete assemblies). Sample list: `/moosefs/guarracino/HPRCv2/PHR_III/HiC/sample_info.tsv`.
- **Code:** `classify_contigs.py`, `trim-telomeres.sh` (telomere detection / trimming).
- **MISSING — needs to be produced:** an authoritative dataset table (sample × haplotype × included/excluded) that resolves the 465-vs-466 discrepancy and is citation-quality for Nature Methods. Current `MANUSCRIPT_DRAFT.md` Methods says "233 individuals × 2 haplotypes = 465 (CHM13 added → 466 in some counts)" — this needs to be a single canonical number with a clear cell.

### C4. Genome-wide identity survey (extended interchromosomal homology comparable to PAR2)
- **Figure assets (already in tree):** `p_genome_wide_identity_heatmap.{pdf,png}`, `p_genome_wide_identity_heatmap_no_inset.{pdf,png}`, `p_genome_wide_numchrom_heatmap.{pdf,png}`, all 24 `identity_heatmap_chr*.pdf` per-chromosome zooms (top-level repo).
- **Data:** `hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz` (top-level, ~genome-wide identity windowed); `chm13.phrs.bed`, `chm13.phrs.no_acro.bed`, `CHM13-HG002.sub-telo-phrs.bed` (PHR call BEDs).
- **Code:** `plot-impg-coverage.R` (top-level repo) — produces the genome-wide identity heatmap.
- **MISSING — needs to be produced:** the explicit "comparable to PAR2" calibration panel (see §1 missing-figures item 3). PHR length distribution from `all-vs-all.p95.id95.len.tsv` exists and is plotted in ED1c — the comparison to PAR2 length is a single number that needs to be added to the figure or caption.

### C5. NJ-tree cladistic analysis (specific clades)
- **Data — distance matrix (input to NJ):** `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv` (41 × 41 arm-level Jaccard distance, exists). `hprcv2.1Mb.subtelo.dist_matrix.tsv` (15,668 × 15,668 sequence-level, exists).
- **Existing tree-like artifacts:** `hprcv2.1Mb.subtelo.arm-upgma-k14.assignments.tsv` (UPGMA k=14 partition, exists). UPGMA dendrogram on top of Fig 1c. Leiden k=15 partition assignments. **No NJ tree files** (`*.nj.*`, `*.newick`, `*.nexus`, `*.tree`) exist in the similarity/ directory or elsewhere in the repo.
- **Code for tree construction:** `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R` builds the UPGMA dendrogram. **No NJ-tree script** is present.
- **MISSING — needs to be produced:** (a) NJ-tree construction from `hprcv2.1Mb.subtelo.arm_dist_matrix.tsv` (R `ape::nj()` or equivalent), (b) annotation of the named clades from the abstract — Xp/Yp via PAR1, Xq/Yq via PAR2, acrocentric short arms, 10p–18p, the {22q,21q,19q,1q,13q,17q} tight clade, 4q–10q DUX4 with copy-number diversity, large moderate-similarity clade — with bootstrap support if feasible. (c) DUX4 copy-number table per haplotype for the 4q–10q clade.

### C6. PCA + community detection on similarity matrix
- **Existing PCA-like artifacts:** MDS / PCoA outputs (`hprcv2.1Mb.subtelo.full_mds.rds`, `hprcv2.1Mb.subtelo.mds.arm-leiden-k15.communities.{pdf,png}`, `hprcv2.1Mb.subtelo.mds.arm-leiden-k15.superpop-hulls.{pdf,png}`); UMAP outputs (`hprcv2.1Mb.subtelo.umap.color-by-{arm,chromosome,superpopulation}.{pdf,png}`). MDS is metric-equivalent to PCoA on the Jaccard distance matrix.
- **Existing community detection:** Leiden k=15 (arm-level), Leiden k=50 (sequence-level), UPGMA k=14. Files in `/moosefs/guarracino/HPRCv2/PHR_III/similarity/`.
- **Population-structure data:** `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_affinity_sequences.tsv`, `cross_arm_superpop_enrichment.tsv`, `fst_superpop_matrix.tsv`. Population labels at `/moosefs/guarracino/HPRCv2/PHR_III/HiC/sample_info.tsv`.
- **MISSING — needs to be produced:** explicit PCA on a feature matrix (rather than PCoA on a distance matrix). Either (a) run PCA on a haplotype × PHR-presence-binary or PHR-identity-fingerprint feature matrix, (b) clarify in the abstract / methods that the existing MDS is the PCA-equivalent (PCoA on Jaccard distance), or (c) compute PCA on the binary haplotype-by-arm community-membership matrix.

### C7. Hi-C 3D maps testing nuclear-envelope-proximity hypothesis
- **Hi-C data (mcool):** `/moosefs/guarracino/HPRCv2/PHR_III/HiC/HG002/{cifi_output,hicpro_output,porec_output}/*.mcool`; `HG00658/hicpro_output/cool/HG00658.mcool`; `NA19036`, `HG02148`, `HG02559`, `RPE1`, `CHM13` similar.
- **Analysis outputs:** `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/{5,10,20,50,100}000bp/<sample>_global_test.tsv`; `no_acrocentric/`, `no_sex/`, `no_acro_p_sex/`, `no_acro_pq_sex/`, `no_strong/` exclusion-control directories under `analysis/human/exclusion_controls/`.
- **Code:** `analyze_hic_communities.py`, `sequence_hic_correlation.py`, `parm_qarm_3d_enrichment.py` (under `scripts/community/`).
- **MISSING — needs to be produced:** (a) explicit nuclear-envelope-proximity test — overlay subtelomeric Hi-C contacts with lamin-associated domain (LAD) calls or a Lamin B1 ChIP track. The abstract specifically frames Hi-C as testing the nuclear-envelope-proximity hypothesis; current Fig 3 / ED5 / ED8 test community-vs-non-community contact enrichment, which is *consistent with* envelope tethering but not a *direct* test of it. (b) ED5d already has within-community reproducibility across samples — strong support for the spatial-clustering claim. The envelope-specific test is the missing piece.

### C8. Synthesizing thesis: ongoing recombination shapes subtelomeres
- **Conceptual support:** ED8a (causal feedback loop schematic — direct visual support); ED8b (D4Z4-CTCF-lamin tethering for C1 — direct mechanistic support). Discussion section of MANUSCRIPT_DRAFT.md §6 is partially on-spec.
- **Direct evidence of ongoing exchange (recombination):** Cross-arm-affinity sequences (`cross_arm_affinity_sequences.tsv`) — 2,484 sequences with cross_arm_affinity > 1; allele-vs-paralog Wilcoxon (`allele_vs_paralog_distance.tsv`) — 5,946 paired distances.
- **MISSING — needs to be produced:** a Discussion paragraph that closes the abstract's title-thesis ("concerted evolution and unorthodox recombination") with a direct citation to (a) the named NJ-tree clades, (b) the within-cluster identity peak, and (c) the Hi-C envelope-proximity test as the integrative model. Currently the Discussion §6 of `MANUSCRIPT_DRAFT.md` is built around feedback-loop / D4Z4 / honest-null-recombination — three of those are reusable, the fourth is off-spec.

---

## 4. Brainstorming inventory (one-pager)

`paper_prep/_brainstorming/` contains 158 files and a categorised README (`paper_prep/_brainstorming/README.md`). Categories per the README:

- **A.** Copy-number-weighted ORA / hypergeometric / phyper methodology — theory, R/Python implementation, validation, benchmarks, edge-case suites, reports (~80 files).
- **B.** Gene-family-specific deep research (OR4F, DUX4 / FRG2, miRNA, OR biology) — 8 .md files.
- **C.** TUBB8 paralog deep research — 5 .md files.
- **D.** g:Profiler / GO / KEGG enrichment runs and outputs — ~20 files.
- **E.** Off-target synthesis docs and supporting tables (executive summaries, decision frameworks, fact-checks, gene-copy tables) — ~10 files.
- **F.** GSEA figures and Excel summary outputs — 4 files (`Figure1_GSEA_BP_vertical.pdf`, `Figure_GSEA_MF_vertical.pdf`, two .xlsx).
- **G.** Wrapper scripts (R, functions) — 2 files.

**Default expectation (per the task brief): nothing in `_brainstorming/` becomes canonical.** The audit confirms this. The narrow exceptions are:

- The **DUX4 copy-number tables** (`gene_copy_summary.csv`, `all_gene_copies_by_arm.csv`, `genome_wide_gene_copies.csv`, `enriched_genes_per_arm.md`) — relevant to the abstract's "4q–10q DUX4 with copy-number diversity" sub-claim of C5. These could be cited as a footnote-level supplementary table to support the DUX4 copy-number diversity statement; but the abstract does not require enrichment / GO context.
- The **deep-research synthesis** (`deep_research_dux4_frg2.md`, `deep_research_synthesis.md`) — provides background on D4Z4 / DUX4 biology that may inform a single sentence in the Discussion, no more.
- **`andrea_phr_reconciliation.md`** — may contain reconciled methodology notes from the co-author; worth a one-pass read for any C2 (implicit pangenome graph) or C5 (cladistic) argument that should be canonical.

Everything else (ORA / phyper / TUBB8 / OR4F / GSEA / GO enrichment) is firmly off-spec for the canonical abstract. The README's note "the canonical paper does not depend on per-family minutiae beyond the abstract's mention of the 4q–10q DUX4 clade" is correct.

---

## 5. Process-failure note

### What happened

The prior task `render-manuscript-draft-4` (commit `3f859cb`) produced `paper_prep/synthesis/MANUSCRIPT_DRAFT.pdf` (now removed by `chore: paper-prep gitignore + remaining brief and survey`, commit `3db7023`) which **passed acceptance** ("≥5 pages, file=PDF") **but contained NO visible figures** when opened by a human reader. The render task's acceptance criterion checked file metadata (page count, file type) but not figure-presence in the rendered PDF.

Additionally, the figures specified in `MANUSCRIPT_SKELETON.md` for ED6 (Dip-C + sperm + RPE-1) and ED7 (mouse meiotic) were **never produced as files** — they exist only as caption stubs in `CAPTIONS.md`. The render task did not detect this either: the absence of `paper_prep/figures/ed6/` and `paper_prep/figures/ed7/` did not surface in any acceptance check.

### Recommended fix for future render tasks

Every future render task in `REWRITE_PLAN.md` (Lane B) **MUST** include the following acceptance criteria in its `## Validation` section, beyond the trivial "file is a PDF" check:

1. **Embedded-image-object count.** Run a programmatic check (e.g., `pdfimages -list <PDF> | wc -l` or `qpdf --json <PDF> | jq` for embedded /Image XObjects) and assert the count is **≥ N**, where N is the expected number of figures in the document. Lane B render tasks must specify N explicitly.
2. **Per-page screenshot.** Render every page to PNG (e.g., `pdftoppm -r 100 <PDF> page`) and verify the screenshots show inline figures at the expected positions. The screenshots should be uploaded as artefacts of the render task so a human / evaluator can spot-check.
3. **Figure-file existence pre-check.** Before the render runs, check that every figure file referenced in the manuscript source actually exists on disk (e.g., `grep -oE 'paper_prep/figures/[a-z0-9_]+/figure_[a-z0-9_]+\.(pdf|png)' <SRC> | xargs -I{} test -f {}`). If any referenced figure is missing, the task **fails** before render rather than producing a silent figure-less PDF.
4. **Caption-figure cross-check.** For each figure file present, verify that a corresponding `caption.md` exists in the same directory. Conversely, every caption referenced in the manuscript main text must resolve to a figure file.

These four checks are mandatory in every Lane B render task in the REWRITE_PLAN.

---

*End of AUDIT_REPORT.md.*
