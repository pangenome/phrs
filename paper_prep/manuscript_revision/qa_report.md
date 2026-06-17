# Manuscript Revision QA Report

Date: 2026-06-17
Task: `manuscript-revision-qa`
Scope: lightweight manuscript-facing QA after the guarded manuscript patch. No heavy analyses were run.

## Summary

The repo-standard manuscript compile was attempted from `submission/` and failed before manuscript processing because the local TeX environment is missing `geometry.sty`. No `submission/paper.pdf` was produced in this QA run.

The active manuscript source audit found no absolute `/moosefs` paths, no internal survey/report citations, and no repo script-path construction pointers in `submission/paper.tex` or `submission/bibliography.bib` under the searched patterns. The only active-source construction markers found are TODO comments in `submission/paper.tex` for author-confirmed front/back matter.

The required revision artifact directory `paper_prep/manuscript_revision/` exists and contains the expected revision package, patch report, audit records, continuum artifacts, and mouse-shape artifacts.

## Commands Run

From the repository root unless noted:

```bash
wg msg read manuscript-revision-qa --agent "$WG_AGENT_ID"
wg show manuscript-revision-qa
git status --short
git log --oneline main..HEAD
ls -la paper_prep/manuscript_revision submission
find paper_prep/manuscript_revision -maxdepth 2 -type f | sort
```

Compile command, from `submission/`:

```bash
make
```

Manuscript-facing path and construction-artifact audits:

```bash
rg -n --hidden --glob 'submission/**' --glob '!submission/paper.log' --glob '!submission/*.aux' --glob '!submission/*.blg' --glob '!submission/*.bbl' '/moosefs|TODO|FIXME|PLACEHOLDER|placeholder|speaker[- ]?note|baked in|construction|TBD|XXX|DRAFT'

rg -n '/moosefs|TODO|FIXME|PLACEHOLDER|placeholder|speaker[- ]?note|baked in|construction|TBD|XXX|DRAFT' submission/paper.tex submission/bibliography.bib submission/README.md submission/BUILD_LOG.md paper_prep/manuscript_revision/*.md

rg -n '/moosefs|TODO|FIXME|PLACEHOLDER|placeholder|speaker[- ]?note|baked in|TBD|XXX|NATURE_DRAFT|SURVEY_|end-to-end-report|\.md:[0-9]|\.md line|scripts/(cladistics|hic|hic-realign|mouse|pedigree)' submission/paper.tex submission/bibliography.bib
```

Artifact inventory commands:

```bash
find paper_prep/manuscript_revision -maxdepth 1 -type f -printf '%f\t%s bytes\n' | sort
find paper_prep/manuscript_revision -mindepth 1 -maxdepth 1 -type d -printf '%f/\n' | sort
find paper_prep/manuscript_revision -type f | wc -l
find paper_prep/manuscript_revision/C0_continuum -maxdepth 1 -type f | wc -l
find paper_prep/manuscript_revision/F3_mouse_shape -maxdepth 1 -type f | wc -l
```

## Build Status

Status: blocked by local TeX dependency.

`make` in `submission/` invokes:

```bash
PAPER="paper" LATEX="pdflatex" BIBTEX="bibtex" LATEX_FLAGS="-interaction=nonstopmode -halt-on-error" bash compile.sh
```

The first `pdflatex` pass fails while loading `jnl.cls`:

```text
! LaTeX Error: File `geometry.sty' not found.
...
!  ==> Fatal error occurred, no output PDF file produced!
make: *** [Makefile:21: paper.pdf] Error 1
```

Output PDF path if available: `submission/paper.pdf`

Result in this QA run: unavailable; `submission/paper.pdf` was not produced because the compile stopped before the first pass completed.

This matches the upstream patch task's environment blocker and does not indicate a newly observed manuscript syntax error.

## Manuscript-Facing Audit

### Active manuscript source

Files audited directly as active manuscript source:

- `submission/paper.tex`
- `submission/bibliography.bib`

Targeted search result for forbidden construction artifacts and private paths:

- No hits for `/moosefs`.
- No hits for `NATURE_DRAFT`, `SURVEY_`, `end-to-end-report`, Markdown line pointers, or the searched local script-path citation patterns.
- No hits for speaker-note or baked-in-placeholder markers.
- Three TODO comments remain in `submission/paper.tex`:
  - `submission/paper.tex:273`: Author Contributions CRediT roles require author confirmation.
  - `submission/paper.tex:279`: Competing Interests require author-confirmed text.
  - `submission/paper.tex:287`: Tables placeholder comment remains.

These TODOs are comments, so they would not print in the PDF, but they remain in the uploadable source.

### Other submission-facing support files

The broader `submission/**` search also found expected internal provenance in support documentation:

- `submission/README.md` references the previous `paper_prep/synthesis/NATURE_DRAFT_v6.md` and notes TODO author-stub status.
- `submission/BUILD_LOG.md` references `NATURE_DRAFT_v6.md` as historical build provenance.
- `submission/notes/20260521_pptx-slide-image-provenance.md` contains multiple `/moosefs` paths and is an internal slide/figure provenance note.
- `submission/scripts/figures/README.md` notes default `/moosefs` paths for figure-generation helpers.
- `submission/bibliography.bib` contains ordinary article-title/note uses of "construction" for pangenome graph methods; these are bibliographic content, not construction artifacts.

Recommendation: do not bundle `submission/notes/` or `submission/scripts/figures/` as public manuscript source-data notes without converting private paths to public sources/accessions. The active `paper.tex` and bibliography do not currently leak these paths.

## Revision Artifact Inventory

`paper_prep/manuscript_revision/` exists and contains 62 files total, including this QA report.

Top-level revision documents:

- `00_inventory.md`
- `01_fanout_graph.md`
- `02_operating_rules.md`
- `A1_sample_counts.md`
- `A2_bibliography_audit.md`
- `A3_artifacts_audit.md`
- `A4_guarracino_doi.md`
- `A5_A7_mechanical_audit.md`
- `A_mechanical_fixes.md`
- `B0_3d_inventory.md`
- `B1_B3_3d_decision_record.md`
- `B4_pvalue_mantel_audit.md`
- `B5_3d_apparatus_essentiality.md`
- `BF_3d_contact_synthesis.md`
- `C0c_D1_resolution_sampling.md`
- `C1_tree_essentiality.md`
- `C2_bootstrap_audit.md`
- `C3_qarm_language.md`
- `CD_continuum_community_synthesis.md`
- `E_pedigree_audit.md`
- `F1_F2_orphan_audit.md`
- `G_abstract_intro_title_package.md`
- `final_revision_package.md`
- `manuscript_patch_report.md`
- `qa_report.md`

Nested artifact directories:

- `C0_continuum/`: 20 files, including continuum reports, assignment/scan tables, diagnostic PDF/PNG, and helper scripts.
- `F3_mouse_shape/`: 5 files, including the mouse shape analysis script, report, input inventory, and output TSV summaries.

Required post-patch artifacts confirmed present:

- `paper_prep/manuscript_revision/final_revision_package.md`
- `paper_prep/manuscript_revision/manuscript_patch_report.md`
- `paper_prep/manuscript_revision/A_mechanical_fixes.md`
- `paper_prep/manuscript_revision/A3_artifacts_audit.md`
- `paper_prep/manuscript_revision/C0_continuum/`
- `paper_prep/manuscript_revision/F3_mouse_shape/`

## Remaining Blockers

1. Local compile environment is incomplete: `geometry.sty` is missing, so `submission/paper.pdf` cannot be regenerated here with the repo-standard command.
2. `submission/paper.tex` still contains TODO comments for author-confirmed Author Contributions and Competing Interests. Agents should not invent these statements.
3. `submission/paper.tex` still contains a TODO comment for optional tables. This is source-only, not printed, but should be resolved or removed before source submission.
4. Public source-data/accession details remain an author/source-data follow-up where the revision package says exact accessions were not found in inspected artifacts.

## QA Conclusion

The guarded manuscript patch passes the lightweight active-source artifact audit for private `/moosefs` paths and internal construction pointers in `paper.tex` and `bibliography.bib`. The final compile remains unvalidated because the local TeX installation lacks `geometry.sty`, and the active source still contains non-printing TODO comments requiring author decisions.
