# synthesize-paper-architect — task brief

This file is the description for the Phase-2 fan-in task. It is referenced by the `wg add` command that creates that task; the agent that picks up the task should read this file in full.

## Context

We are preparing **two coupled deliverables** from the existing end-to-end report at `end-to-end-report/` and accompanying repo artifacts:

1. A **Nature manuscript** on subtelomeric sequence sharing across HPRCv2 (233 samples, 465 haplotypes, 15 arm-level + 50 sequence-level communities) validated by 3D genome data across multiple technologies, cell types, and species.
2. A **15-minute conference presentation** keyed to the manuscript's main figures.

The work is organized as a divide-and-conquer pipeline:

- **Phase 1 (done by the time you read this):** ~14 parallel survey tasks that each digested one section of the report (`SURVEY_*.md` in `paper_prep/surveys/`) plus a figure inventory, a data-table inventory, and a cytogenetic/FISH literature framing survey.
- **Phase 2 (you):** read every Phase-1 survey output and produce the manuscript skeleton + talk outline + work decomposition.
- **Phase 3 — figures (you create these tasks):** implementation agents that mostly *generate figures* from existing data — no large new analyses unless flagged. Each agent gets the rough draft + its own slot. Output: figure files (PNG/PDF/SVG) under `paper_prep/figures/` plus per-figure caption .md.
- **Phase 4 — section drafting (you create these tasks):** one agent per manuscript section file in `paper_prep/synthesis/sections/`. Each agent receives (a) the rough manuscript skeleton (their slot highlighted), (b) the relevant Phase-3 figures + captions, (c) the relevant Phase-1 SURVEY outputs and the source `end-to-end-report/report/*.md` files. They produce tight publication-quality Markdown prose for that section with DOI/URL citations and figure references.
- **Phase 5 — validation (you create these tasks):** fine-tooth-comb agents, one per drafted section. They check every numeric/statistical/literature claim against (a) the data we actually have on disk, (b) the cited literature. They produce `paper_prep/synthesis/validation/VALIDATE_<section>.md` with a per-claim PASS/FAIL/REVISE list. Validation runs after the section draft AND its figures both exist.
- **Phase 6 — synthesize + harmonize + render (final integrator):** this is **NOT** a pure file concatenation. The integrator agent reads every section draft AND every validation report, then *actively edits* the section files to:
  - Resolve every FAIL/REVISE item from Phase 5 validation reports.
  - Harmonize voice, tense, and terminology across sections (e.g. "PHR" vs "subtelomeric homologous region" — pick one and apply globally).
  - De-duplicate claims that appear in multiple sections; route each to the correct home (intro vs results vs discussion).
  - Sand transitions between sections and tighten the abstract so it accurately summarizes the final body.
  - Verify every figure is referenced in prose at least once and every prose figure-reference resolves to a real figure file.
  - Verify every cited DOI/URL is non-broken (curl -I or wget --spider checks).
  - **Only after harmonization** concatenate the edited sections into `DRAFT_v1.md` and render `DRAFT_v1.pdf` via pandoc (xelatex → wkhtmltopdf → weasyprint fallback).
  - Produce a `paper_prep/synthesis/INTEGRATOR_NOTES.md` listing every change made during harmonization, so the user sees what was reconciled.

## Inputs you must read

All files in `paper_prep/surveys/SURVEY_*.md`. List them with `ls paper_prep/surveys/`. Each file follows a 7-section template:

1. Key findings with metrics
2. Existing figures referenced (paths)
3. Existing CSVs/data files (paths)
4. Methods used
5. Open gaps for figure (re)generation
6. Suggested figures (main + extended/SI) with caption ideas, marked produced-vs-todo
7. Suggested talk slide takeaways

The framing survey (`SURVEY_FRAMING_cytogenetic_fish.md`) has a different structure — read it for how to position this work against prior FISH / flow-cytometry / cytogenetic literature on subtelomeric exchange (Riethman, Mefford, Linardopoulou, Stong, Trask, Bailey, etc.).

## Deliverables (produce all four)

### A. `paper_prep/synthesis/MANUSCRIPT_SKELETON.md`

**Authoring format (mandatory).** The paper is written as a *set* of Markdown files — one per top-level manuscript section — that are later concatenated into a single PDF via pandoc for review. Typst/LaTeX migration is deferred (Nature accepts a PDF for submission, so the source format does not matter to them). Plan the section breakdown so it maps cleanly onto separate .md files: e.g. `paper_prep/synthesis/sections/00_abstract.md`, `01_introduction.md`, `02_results_communities.md`, `03_results_3d.md`, `04_results_pedigree.md`, `05_discussion.md`, `06_methods.md`, `07_references.md`. The integrator at the end of the workflow stitches them.

Nature-format skeleton:

- **Title** (1-2 candidate titles)
- **Abstract** — ~150 words, complete prose draft
- **Main text** — ≤3000 words target. Section headings + bullet outlines for each:
  - Introduction (with the cytogenetic/FISH framing — what prior work could see, what we now see at scale)
  - Results (4-7 results subsections, each keyed to one main figure)
  - Discussion (the meiotic-bouquet / D4Z4-CTCF-lamin / nucleolar-association integration; testable predictions; limitations)
- **Main figures** — exactly **4 figures** target (Nature norm: 4-6 main; aim low to keep the story tight). For each: figure number, working title, panel-by-panel description (A, B, C…), data source, status (already-produced / needs-regeneration / needs-creation), one-line caption.
- **Extended Data figures** — ~6-10. Same per-figure schema, lower bar.
- **Methods** — bullet outline only, citing the existing report sections that already contain method text.
- **References** — list of ~30-50 must-cite papers with author + year + topic. Include the cytogenetic/FISH framing citations.

### B. `paper_prep/synthesis/TALK_OUTLINE_15MIN.md`

15-min talk: ~12-14 slides. For each slide: number, title, one-line takeaway, primary figure (cite the main-fig number from MANUSCRIPT_SKELETON.md or Extended fig). Include opening framing (cytogenetic-era → flow+FISH → HPRC pangenome) and closing predictions/future-work.

### C. `paper_prep/synthesis/WORK_DECOMPOSITION.md`

Four ordered lists:

- **Phase 3 (figure tasks).** One row per main + extended figure. Columns: figure-id, source-data path(s), tool-or-script suggestion, dependencies, acceptance criteria. Mark `STATUS=READY|REGENERATE|GENERATE|GAP` per the figure-status discipline below.
- **Phase 4 (section-draft tasks).** One row per file under `paper_prep/synthesis/sections/`. Columns: section-id, target file path, manuscript-skeleton excerpt that defines the slot, input figures (Phase-3 ids), input surveys, target word-count.
- **Phase 5 (validation tasks).** One row per Phase-4 section. Columns: section-id, claims-to-verify (numbered, each pointing at a specific data file or literature ref), output-file path `paper_prep/synthesis/validation/VALIDATE_<section>.md`.
- **Phase 6 (integrator).** One row — the single harmonize+render task with the full --after list.

### D. Dispatch Phases 3, 4, 5, and 6 via `wg add`

The downstream pipeline has four fan-outs and one fan-in. Wire them up like this:

- **Phase 3 (figures)** — `--after synthesize-paper-rough`. Tag `paper-prep,figure`.
- **Phase 4 (section drafts)** — `--after synthesize-paper-rough` AND `--after` the relevant Phase-3 figure task(s) for that section. Tag `paper-prep,section-draft`. One task per file in `paper_prep/synthesis/sections/` (abstract, intro, each results section, discussion, methods).
- **Phase 5 (validation)** — `--after` the matching Phase-4 section-draft task AND `--after` the matching Phase-3 figure task(s). Tag `paper-prep,validation`. One task per drafted section. Output `paper_prep/synthesis/validation/VALIDATE_<section>.md`.
- **Phase 6 (integrator: synthesize + harmonize + render PDF)** — `--after` *every* Phase-4 and Phase-5 task. Tag `paper-prep,integrator`. Single task. Outputs: edited section files (in-place), `paper_prep/synthesis/DRAFT_v1.md`, `paper_prep/synthesis/DRAFT_v1.pdf`, `paper_prep/synthesis/INTEGRATOR_NOTES.md`.

Every task description must include: (1) a one-line summary, (2) the relevant excerpt of MANUSCRIPT_SKELETON.md so the agent knows the slot it is filling, (3) explicit input file paths (Phase-1 surveys, source report, sibling figures, etc.), (4) explicit output file path, (5) acceptance criteria.
The Phase-6 integrator (`wg add "Synthesize, harmonize, and render manuscript v1" --after <all phase-4 and phase-5 ids> --tag paper-prep,integrator`) does the following work, in order:
  1. Read every section draft AND every validation report.
  2. Resolve every FAIL/REVISE item from validation by editing the section files in place.
  3. Harmonize voice, tense, and terminology across sections.
  4. De-duplicate claims that appear in multiple sections; route each to its correct section.
  5. Sand transitions between sections; tighten the abstract so it matches the final body.
  6. Verify every figure is referenced in prose at least once, and every figure-reference resolves to a real figure file.
  7. Verify every cited DOI/URL resolves (curl -I or wget --spider).
  8. **Only then** concatenate the edited sections (numeric-prefix order) into `paper_prep/synthesis/DRAFT_v1.md`.
  9. Render `paper_prep/synthesis/DRAFT_v1.pdf` via `pandoc --pdf-engine=xelatex --citeproc` (fall back to `wkhtmltopdf`, then `weasyprint` — log which succeeded).
 10. Write `paper_prep/synthesis/INTEGRATOR_NOTES.md` listing every change made during harmonization (so the user can review what was reconciled).
 11. Record DRAFT_v1.md, DRAFT_v1.pdf, and INTEGRATOR_NOTES.md as `wg artifact`s on the integrator task.

**Reference style (mandatory).** Every section .md and the merged DRAFT_v1.md must use Markdown citations of the form `[Author Year](https://doi.org/...)` or `[Author Year](https://www.ncbi.nlm.nih.gov/pubmed/...)`. Every claim drawn from prior literature must have such a link. The MANUSCRIPT_SKELETON.md reference list (your deliverable A) should already collect these DOI/URL pairs so downstream Phase-3/4 agents can reuse them. (Format-portability note: this DOI-link style survives pandoc → typst → LaTeX conversion later when Nature wants a different source.)

**Figure status discipline.** WORK_DECOMPOSITION.md must clearly separate:
  - **READY (already produced)** — exists at a concrete path; only needs polish/relabel/inclusion. Cite the path.
  - **REGENERATE** — exists but needs different parameters, layout, or styling for publication.
  - **GENERATE** — does not exist; needs to be made from existing data files (cite which CSV/RData/etc).
  - **GAP** — would require new analysis; flag clearly and explain why existing data is insufficient.

**Respect the decomposition guardrails.** `max_child_tasks_per_agent` defaults to 10 — if you need more than ~10 sub-tasks, group several figures into a single task (e.g. "Figure 3 + Extended 5 + Extended 6: HPRC community Hi-C panels"). Aim for ~8-15 total Phase-3/4 tasks. The integrator task makes ~16. If you need to exceed the cap, raise it via `wg config --max-child-tasks 20` first.

## What NOT to do

- Do **not** propose new wet-lab experiments or new sequencing.
- Do **not** propose any large new computational analysis without flagging it as a `GAP` in WORK_DECOMPOSITION.md and explaining why the existing data is insufficient. Default assumption: every claim in the report is already supported by data we have.
- Do **not** write the full prose of the manuscript yourself — that is the integrator's job after Phase 3+4 land.
- Do **not** invoke `Task`/`TaskCreate` built-in tools. Use `wg add` only. (See CLAUDE.md.)

## Acceptance criteria

- All four deliverables (A-D) exist at the paths above.
- MANUSCRIPT_SKELETON.md has exactly 4 main figures and 6-10 Extended Data figures listed with full panel descriptions.
- TALK_OUTLINE_15MIN.md has 12-14 slides each keyed to a manuscript figure.
- WORK_DECOMPOSITION.md has every figure tagged READY-vs-GENERATE.
- All Phase-3/Phase-4 tasks have been dispatched via `wg add` and are visible in `wg list`.
- A final `compile-manuscript-draft` task exists with the right `--after` chain.
