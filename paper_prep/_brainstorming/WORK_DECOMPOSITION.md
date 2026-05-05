# WORK DECOMPOSITION — Phase-3 (figures) and Phase-4 (validation)

**Inputs.** `MANUSCRIPT_SKELETON.md` + `TALK_OUTLINE_15MIN.md` + 13 `paper_prep/surveys/SURVEY_*.md`.
**Constraint.** Under 15 dispatched tasks total (architect brief). Reached at **13 worker tasks + 1 integrator = 14 dispatched tasks**.
**Pattern.** Fan-out-merge — workers in parallel; `compile-manuscript-draft` integrator depends `--after` all workers.
**File-scope rule.** Each task owns a unique pair of (figure-file path, caption-stub file path). Each task description includes the relevant excerpt of `MANUSCRIPT_SKELETON.md` so workers can render captions in context.

---

## Phase-3: figure generation (9 tasks)

> Each task produces one figure deliverable (main figure or ED figure) plus a caption stub committed under `paper_prep/figures/<id>/`. Figure deliverable is the renderable PDF/PNG (preferred) plus the rendering script. Tag: `paper-prep,figure`.

| # | Task ID | Figure(s) | Status | File scope |
|---|---|---|---|---|
| 1 | `figure-1-landscape-communities` | Fig 1 (4 panels) | 2 READY + 2 GENERATE | `paper_prep/figures/fig1/` |
| 2 | `figure-2-heterogeneity-population` | Fig 2 (4 panels) | 1 READY + 3 GENERATE | `paper_prep/figures/fig2/` |
| 3 | `figure-3-3d-convergence` | Fig 3 (4 panels) | 1 READY + 3 GENERATE | `paper_prep/figures/fig3/` |
| 4 | `figure-4-pedigree-mouse` | Fig 4 (4 panels) | 2 READY + 2 GENERATE | `paper_prep/figures/fig4/` |
| 5 | `figure-ed1-ed2-pipeline-seqlevel` | ED1 + ED2 | 1 READY + 7 GENERATE | `paper_prep/figures/ed1/`, `paper_prep/figures/ed2/` |
| 6 | `figure-ed3-ed4-annotation-genes` | ED3 + ED4 | 1 READY + 7 GENERATE | `paper_prep/figures/ed3/`, `paper_prep/figures/ed4/` |
| 7 | `figure-ed5-hic-robustness` | ED5 (4 panels) | 0 READY + 4 GENERATE | `paper_prep/figures/ed5/` |
| 8 | `figure-ed6-ed7-singlecell-mouse` | ED6 + ED7 | 4 READY + 4 GENERATE | `paper_prep/figures/ed6/`, `paper_prep/figures/ed7/` |
| 9 | `figure-ed8-discussion-models` | ED8 (4 panels; 2 schematics) | 0 READY + 4 GENERATE | `paper_prep/figures/ed8/` |

## Phase-4: validation, captions, methods (4 tasks)

| # | Task ID | Deliverable | Tag |
|---|---|---|---|
| 10 | `validate-statistics-fdr-cis` | Statistics audit: BH-FDR coverage + bootstrap CIs + effect-size table for all main-text claims | `paper-prep,validation` |
| 11 | `validate-reproducibility-paths` | Reproducibility audit: vendor `/moosefs/.../scripts/` to repo (or document explicit paths); pin tool versions; verify every figure-source path is reachable | `paper-prep,validation` |
| 12 | `validate-captions-references` | Caption + reference cross-check: every figure caption compiles, every cited number traceable to a TSV row, bibliography cleaned | `paper-prep,validation` |
| 13 | `validate-acceptance-checklist` | Final acceptance walk-through: §1.4–§1.12 anchoring metrics each have a figure or table; SI tables present; novel-contributions ledger filled | `paper-prep,validation` |

## Integrator (1 task)

| # | Task ID | Deliverable | Depends on |
|---|---|---|---|
| 14 | `compile-manuscript-draft` | `paper_prep/synthesis/MANUSCRIPT_DRAFT.md` (rough draft, ~3,000 words) integrating all 9 figure tasks + 4 validation tasks + skeleton + survey content | `--after` all 13 above |

---

## Per-task acceptance criteria (excerpts that go into each `wg add -d`)

### Phase-3 figure tasks

Each figure task must produce, under `paper_prep/figures/<id>/`:

1. `figure_<id>.pdf` — final composed figure (or links to existing READY PDFs in a manifest).
2. `figure_<id>.{R,py}` — rendering script (vendored or with explicit `/moosefs/...` source path).
3. `caption.md` — preliminary caption (≤ 200 words) using metrics drawn from named TSVs in the corresponding survey.
4. `sources.tsv` — table of (panel, status, source-PDF-or-TSV-path) with status ∈ {READY/GENERATE}.
5. `## Validation` block in the task description checked off:
   - [ ] All panels in `MANUSCRIPT_SKELETON.md` Fig X (or EDX) entry are addressed (one panel per row in `sources.tsv`)
   - [ ] If GAP (no extant artefact), task adds a flag in this `WORK_DECOMPOSITION.md` `## Gaps` section (don't run new wet-lab or large new analyses — see Guardrails)
   - [ ] Caption cites at least 2 metrics with source TSV paths
   - [ ] PDF + PNG both produced (PNG for talk re-use)

### Phase-4 validation tasks

#### 10. `validate-statistics-fdr-cis`
- [ ] Every p-value in `MANUSCRIPT_SKELETON.md` headline-numbers section has either a corrected-q-value column or an explicit "uncorrected" annotation
- [ ] Every odds-ratio (e.g. f7501 chr16_q OR = 17.4) has a 95 % CI (Fisher exact mid-p or bootstrap)
- [ ] Multi-resolution Mantel ρ aggregated to a single SI table (5/10/20/50/100 kb × 8 datasets)
- [ ] Output: `paper_prep/synthesis/STATS_AUDIT.md` summarising additions; per-test corrections written into the source TSVs alongside

#### 11. `validate-reproducibility-paths`
- [ ] Every script path cited in surveys (`/moosefs/guarracino/HPRCv2/scripts/...`) either copied into `scripts/` of this repo, OR listed explicitly in a `paper_prep/synthesis/SCRIPT_INVENTORY.md` with absolute paths
- [ ] Tool versions pinned (wfmash, pggb, odgi, impg, igraph/Leiden, R packages) → `paper_prep/synthesis/VERSIONS.md`
- [ ] Every figure-source path in `MANUSCRIPT_SKELETON.md` re-verified to exist on disk; missing files added to `## Gaps` section below
- [ ] `rpe1_subtelo.communities.tsv` vs `rpe1.communities.tsv` filename inconsistency resolved (`SURVEY_09 §5 #1`)

#### 12. `validate-captions-references`
- [ ] Every caption in `MANUSCRIPT_SKELETON.md` has p-value, n, and effect-size if applicable
- [ ] Bibliography assembled: all named-author citations from `SURVEY_FRAMING Part 4` plus Riethman 2004, Stout 1999, van Deutekom 1996, Cechova et al. 2025, Francis et al. 2025, Porubsky et al. 2025, Xu et al. 2025, Zuo et al. 2021, Tan et al. 2018, Lemmers 2010, Masny 2004, Ottaviani 2009, Patel 2019, Lalli 2025, Gershman 2022 (Stergachis Fiber-seq paper to be tracked down)
- [ ] Output: `paper_prep/synthesis/REFERENCES.bib` and `paper_prep/synthesis/CAPTIONS.md`

#### 13. `validate-acceptance-checklist`
- [ ] Each of the 12 anchoring findings in `SURVEY_10/11/12 §1.2` is supported by a Main / ED panel — table of which figure
- [ ] 27 novel contributions ledger written as `paper_prep/synthesis/NOVEL_CONTRIBUTIONS.tsv` (5 columns: id / claim / metric / figure / TSV)
- [ ] Limitations × findings cross-reference (`SURVEY_10/11/12 §6 T-5`) written as `paper_prep/synthesis/LIMITATIONS_X_FINDINGS.tsv`
- [ ] Output: `paper_prep/synthesis/ACCEPTANCE_CHECKLIST.md` with all rows ✓ or with a flag pointing to the gap entry

### Integrator task

#### 14. `compile-manuscript-draft`
- Reads all `paper_prep/figures/*/caption.md` + `paper_prep/synthesis/{MANUSCRIPT_SKELETON.md,STATS_AUDIT.md,REFERENCES.bib,ACCEPTANCE_CHECKLIST.md}` + the 13 surveys.
- Produces `paper_prep/synthesis/MANUSCRIPT_DRAFT.md`:
  - Abstract (~200 words)
  - 6 main-text sections (Intro, Communities, Heterogeneity, 3D, Pedigree+Cross-species, Discussion) following the section outline in `MANUSCRIPT_SKELETON.md`
  - Inline figure callouts `[Fig 1]..[Fig 4]` and `[ED1]..[ED8]` and `[Table 1]`
  - Methods stubs (~500 words; pointers to script inventory)
  - Bibliography appended from `REFERENCES.bib`
- Acceptance:
  - [ ] Word count main text 2,500–3,500
  - [ ] Every figure / ED / table referenced ≥ once
  - [ ] No claim without a TSV / figure citation
  - [ ] Honest-null section (recombination ρ = 0.00 after confound removal) explicitly present in Discussion

---

## Gaps and out-of-scope items (flagged per architect guardrails)

The following are **not** new tasks (would violate "no new wet-lab or large new computational analyses" guardrail) but are surfaced here so reviewers and downstream agents can see them:

1. **Human meiotic Hi-C** — flagged in `SURVEY_07 §5 #1` as the "single most informative missing experiment". Belongs in Discussion as named future work.
2. **CTCF-density × Hi-C contact correlation** — datasets exist (Gershman 2022 T2T-CHM13 ENCODE; Stergachis Fiber-seq) but the cross-correlation has not been computed (`SURVEY_07 §5 #4`). Could be a follow-up paper; noted in ED8 as a testable prediction.
3. **F1-hybrid (B6 × CAST) phased mouse Hi-C** — `analyze_hic_communities.py` advertises support but the source uses non-F1 Zuo 2021 (`SURVEY_08 §5 #6`).
4. **De-novo vs inherited split** of pedigree patches — data are present in `patches.tsv`, but a per-patch parental-presence call is not in current outputs (`SURVEY_14 §5 #7`).
5. **Source-stratified (LCL vs blood) cross-arm validation** — `hprc-sequence-production.tsv` lacks DNA-source column (`SURVEY_10/11/12 §4.3 T5`).
6. **Stability bootstraps for Leiden community assignments** — `SURVEY_01 §5 #2`. Could be added inside `figure-ed5-hic-robustness` if budget allows; otherwise belongs in SI Methods.
7. **Repeat-class composition of mouse PHRs** (RepeatMasker / Ensembl GFF intersected with PHR) — `SURVEY_08 §5 #1`. Annotations downloaded but not intersected; flagged for the Methods footnote, not for a new analysis task.
8. **PHR-only re-run of the gene-enrichment analysis** (Angela's GSEA used 1 Mb window; PHR-only re-run pending per `TODO.md`) — flagged for the gene-enrichment figure caption; if the PHR-only re-run materialises before manuscript freeze, swap the panels.
9. **BH-FDR correction on f7501 per-arm × per-superpopulation enrichment** (~65 tests) — `SURVEY_01 §5 #1`. Captured by `validate-statistics-fdr-cis` (task 10) — counts as validation, not new analysis.
10. **Filename mismatch `rpe1.communities.tsv` vs `rpe1_subtelo.communities.tsv`** in `SURVEY_09 §5 #1` — **RESOLVED** by `validate-reproducibility-paths` (agent-685). Both files exist on disk (`/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/`); they are independent Leiden runs (37 communities `rpe1.*` 2026-03-24 — current; 39 communities `rpe1_subtelo.*` 2026-03-20 — earlier sweep). **Manuscript-canonical name is the `rpe1.*` prefix.** See `paper_prep/synthesis/SCRIPT_INVENTORY.md §13` for the full reasoning. Documentation TODO: add a `README.md` inside the `RPE1_subtelo/` directory disambiguating the two prefixes.

### Gaps surfaced by the reproducibility audit (2026-05-05)

The following figure-source paths in `MANUSCRIPT_SKELETON.md` exist conceptually but require either a path-canonicalisation, a new TSV emission, or a schematic / composite that has not yet been authored. They are NOT new analyses (per architect guardrails) and are tracked here so figure-task workers can attend to them or surface them as task-internal sub-blockers.

11. **Fig 1c — `arm-leiden-k15.assignments.tsv` filename canonicalisation.** The canonical file on disk is `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv` (15 communities; per `detect_communities.R` output naming). The skeleton's `arm-leiden-k15.assignments.tsv` name is descriptive, not literal. Figure-1 worker should cite the actual filename and verify the `k=15` content matches the headline numbers.

12. **ED2b — `similarity.tsv.gz` per community subsets.** SURVEY_01 §1 cites a ~10.8 GB compressed similarity TSV under `/moosefs/guarracino/HPRCv2/PHR_III/pggb/.../similarity.tsv.gz` but does not give the exact subdirectory. ED2b worker must locate the specific pggb-output subdirectory (Andrea/Pjotr's pggb run) before the within-community Jaccard distance histogram can be generated.

13. **ED3c — `.telo.tsv` canonical path.** SURVEY_02 §5 #1 explicitly flags this gap. The terminal-telomere TSV used for the Kruskal-Wallis test exists upstream of `hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz` (likely RUKKI / `seqtk-telo` output) but its absolute path is not documented. ED3 worker must locate or re-derive this TSV.

14. **ED4d — per-arm pseudogene fraction TSV.** OR4F pseudogenisation gradient (62.1 % pseudogene; 11.1 % chr7_p → 99.8 % chr15_q) is reported in SURVEY_10/11/12 but the per-arm pseudogene-fraction TSV is not yet emitted. ED4 worker may need to derive it from the Liftoff GFF3 (biotype filter) — small computation but pure derivation, not a new analysis.

15. **ED6c — per-arm radial-position TSV.** Table embedded in SURVEY_06 §1.4 (46 arms, 3 radial categories). No per-arm-radial TSV is emitted by current scripts; ED6 worker should export the table from `community_3d_enrichment.py` or copy directly from §1.4.

16. **ED7d — mm39 → hg38 syntenic net for mouse private pairs.** SURVEY_08 §1.6 lists the mouse private pairs but the explicit syntenic-net TSV mapping mm39 mouse-arm pairs to hg38 syntenic intervals is not produced. ED7 worker may need to derive from a UCSC syntenic-net BED (small derivation).

17. **ED8d — HG002 100 kb compartment-eigenvector TSV.** SURVEY_07 §6 T-6 references a compartment-identity-at-tips diagnostic (e1 distribution, 68 % A; mean +0.007) but the eigenvector TSV's path is not in the surveys. ED8 worker must locate the cooltools-`call_compartments` output for HG002 at 100 kb.

18. **ED1a, ED1d, ED8a, ED8b — schematics / composites.** Pipeline schematic, NA18982#1 chr18_q chimera evidence, causal-feedback-loop schematic, and D4Z4-CTCF-lamin tethering schematic are all "new authoring" panels (Inkscape / Affinity / composite). Workers should produce them as the last step of their respective figure tasks; not blockers.

19. **`SURVEY_01 §5 #11` version-pinning TODO.** pggb / odgi / igraph / R-package versions must be captured from a fresh `sessionInfo()` and Guix profile manifest before manuscript freeze. Captured by `validate-reproducibility-paths` in `paper_prep/synthesis/VERSIONS.md` §7 as a Methods TODO; downstream `compile-manuscript-draft` should flag in Methods.

20. **`SURVEY_02 §5 #2` Liftoff version pinning.** HPRC annotation-index Liftoff version (used for 462 of the 465 haplotypes) must be added alongside the JHU v0.6 already documented for HG002. Methods TODO; tracked in `VERSIONS.md`.

21. **`SURVEY_01 §1.5` minimap2 version pinning.** Cross-aligner agreement was run with minimap2 v2.30; the current audit host runs 2.24-r1122. Methods must explicitly cite v2.30 for the chimera-flagging analysis to remain reproducible. Tracked in `VERSIONS.md`.

22. **wfmash `do_not_overfilter` feature branch (RPE-1).** SURVEY_09 §1.1 used `wfmash` from the `do_not_overfilter` feature branch (single-prefix workaround). For publication, this should be replaced by either the merged-canonical wfmash equivalent or an explicit commit-hash pin in the Methods. Tracked in `VERSIONS.md`.

---

## Why the architect grouped some figures into shared tasks

Grouping rationale (to keep the count under 15):

- **ED1 + ED2** share the pipeline-and-similarity data tree (`/moosefs/.../similarity/`); same scripts, same input TSVs.
- **ED3 + ED4** share the annotation/repeat data tree (`/moosefs/.../enrichment/`, `/moosefs/.../ttaggg_analysis/`, repo-root gene-enrichment artefacts).
- **ED6 + ED7** share the single-cell + mouse PHR-coordinate analysis pattern (Mantel + radial + per-arm-pair) — same plotting conventions.

A worker handling a grouped task should produce the panels as separate PDFs but commit them under the two ED subdirectories cleanly. This avoids the same agent re-reading the same source TSVs twice.

---

## Why no figures are tagged "split into more sub-tasks"

Each Phase-3 task is sized to a single rendering session: 4 panels per task, all panels share the same source data tree and rendering toolchain. Sub-task fanout would only multiply coordination overhead. The integrator (task 14) is the only join point.

If a worker hits a real blocker on one panel (e.g., source TSV missing), they should `wg add` a focused fix-task rather than re-decomposing the figure. The fix-task should list the missing artefact and unblock the parent figure.

---

## Dispatch summary

```
Phase 3 (9 tasks):
  figure-1-landscape-communities          (F1, 4 panels)
  figure-2-heterogeneity-population       (F2, 4 panels)
  figure-3-3d-convergence                 (F3, 4 panels)
  figure-4-pedigree-mouse                 (F4, 4 panels)
  figure-ed1-ed2-pipeline-seqlevel        (ED1+ED2, 8 panels)
  figure-ed3-ed4-annotation-genes         (ED3+ED4, 8 panels)
  figure-ed5-hic-robustness               (ED5, 4 panels)
  figure-ed6-ed7-singlecell-mouse         (ED6+ED7, 8 panels)
  figure-ed8-discussion-models            (ED8, 4 panels)

Phase 4 (4 tasks):
  validate-statistics-fdr-cis
  validate-reproducibility-paths
  validate-captions-references
  validate-acceptance-checklist

Integrator (1 task):
  compile-manuscript-draft (--after all 13 above)

Total dispatched: 14 (architect cap: 15) ✓
```
