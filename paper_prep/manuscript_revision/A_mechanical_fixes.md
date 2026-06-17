# A Mechanical Fixes Fan-in

Date: 2026-06-17
Task: `manuscript-revision-a-fanin`

## Scope

Integrated the completed Cluster A audits:

- `A1_sample_counts.md`
- `A2_bibliography_audit.md`
- `A3_artifacts_audit.md`
- `A4_guarracino_doi.md`
- `A5_A7_mechanical_audit.md`

Only mechanical corrections supported by those audits were applied. J-marked
framing questions and author-judgment items were not resolved here.

## Changed Files

- `submission/paper.tex`
- `submission/bibliography.bib`
- `submission/compile.sh`
- `submission/README.md`
- `submission/BUILD_LOG.md`
- `paper_prep/manuscript_revision/A_mechanical_fixes.md`

## Line-level Rationale

### `submission/paper.tex`

| Current line(s) | Change | Rationale |
| ---: | --- | --- |
| 16-18 | Removed `multibib` setup. | A2 found duplicate main/Methods references caused by separate `\cite` and `\citeMeth` lists. A single numbered reference list removes the duplicate printed entries mechanically. |
| 33 | Replaced the source comment pointing to `submission/scripts/figures/` with a public-repository wording. | A3 flagged construction paths in submission source comments as source-cleanliness issues. |
| 49 | Changed abstract denominator to `465 near-complete assemblies (464 HPRC v2 haplotypes from 232 individuals, together with CHM13)`. | A1 confirmed 232 HPRC individuals, 464 HPRC haplotypes, and 465 analyzed assemblies including CHM13. |
| 96 | Changed HPRC v2 population denominator from 233 to 232 individuals. | A1 count reconciliation. |
| 107-108 | Changed the Results denominator from 233 individuals / 465 HPRC haplotypes to 232 individuals / 464 haplotype assemblies plus CHM13. | A1 count reconciliation. |
| 301-302 | Changed Fig. 1 legend from 465 HPRC v2 haplotypes plus CHM13 to 464 HPRC v2 haplotypes plus CHM13. | A1 count reconciliation. |
| 331-332 | Replaced `Leiden $k = 15$ community` with `15-community Leiden arm-level partition`. | A5/A7 symbol-collision audit: reserve `k_{\mathrm{NN}}` for nearest-neighbor degree and describe community counts in words. |
| 427-432 | Changed sample-selection Methods text to 232 HPRC individuals, 464 HPRC haplotype assemblies, 465 total assemblies with CHM13, and removed the incorrect superpopulation count breakdown. | A1 found the prior superpopulation counts summed to 233 and should not be retained without recomputation. |
| 437-443 | Removed internal `01_pipeline.md` provenance pointers from sample-exclusion prose while preserving the exclusion facts. | Same construction-artifact class as A3 line-pointer findings; submission text should state the method/fact directly rather than cite internal Markdown. |
| 456, 547, 579, 583-586, 595, 607, 616, 620 | Converted remaining Methods `\citeMeth{...}` commands to `\cite{...}`. | A2 single-reference-list cleanup. |
| 465 | Replaced the IMPG Methods citation target with `pangenome_graphs_impg_IMPG2023`. | A2 found the previous ODGI and Minigraph-Cactus keys were wrong targets for the IMPG transitive-closure sentence. |
| 476 | Removed the likely stray VG citation from PHR output-count reporting. | A2 found `Garrison2018` did not support the reported IMPG sliding-window output counts. |
| 482 | Added the ODGI citation key to the PGGB/odgi Jaccard subsection. | A2 identified ODGI as misplaced in the IMPG subsection and appropriate for the `odgi similarity/layout/draw` sentence. |
| 493-499 | Rewrote Leiden parameter text to use `k_{\mathrm{NN}}` for nearest-neighbor degree and describe 15 arm communities / 50 sequence communities in words. | A5/A7 symbol-collision audit. |
| 503-504 | Replaced bare `$k = 14$` with `tree cut into 14 clusters`. | A5/A7 symbol-collision audit. |
| 510-513 | Replaced local bootstrap script path and runtime note with direct bootstrap-method prose plus public-repository implementation statement. | A3 construction-artifact cleanup. |
| 540-543 | Removed the strict-MAPQ local script path and retained the method/control result. | A3 construction-artifact cleanup. |
| 566-569 | Replaced `scripts/hic-realign/` pointer and reproduction warning with public-repository workflow wording. | A3 construction-artifact cleanup. |
| 587-589 | Removed `05_hic_validation.md` line pointer and cited the public Dip-C source while retaining the statistic. | A3 construction-artifact cleanup. |
| 607-608 | Removed local mouse Mantel script path and stated that implementations are in the public repository. | A3 construction-artifact cleanup. |
| 614-617 | Removed local pedigree-classification script path and described classification basis in prose. | A3 construction-artifact cleanup. |
| 621-624 | Removed local Monte Carlo script path and stated that the permutation implementation is in the public repository. | A3 construction-artifact cleanup. |
| 634-643 | Replaced private `/moosefs` roots and repo-relative Hi-C realignment path in Data and Code Availability with public GitHub and public-release/accession wording. | A3 highest-priority path-leak cleanup. |
| 660 | Moved the single bibliography to the end of the main manuscript/Methods block. | A2 recommended replacing the early main bibliography plus Methods bibliography with one final reference list. |

### `submission/bibliography.bib`

| Current line(s) | Change | Rationale |
| ---: | --- | --- |
| 24-32 | Corrected `acrocentric_Guarracino2025ape` metadata while keeping the existing key. | A4 verified DOI `10.64898/2025.12.22.696095` is real but the first author/date metadata were wrong. Keeping the existing key avoids a cite-command rename in this mechanical pass. |

### `submission/compile.sh`

| Current line(s) | Change | Rationale |
| ---: | --- | --- |
| 9 | Clean stale `Meth.*` and `Supp.*` auxiliaries. | Prevents obsolete multibib artifacts from surviving after conversion to one reference list. |
| 14-21 | Removed `bibtex Meth` and `bibtex Supp` passes; build sequence is now `pdflatex -> bibtex paper -> pdflatex -> pdflatex`. | Required by the A2 single-reference-list conversion in `paper.tex`. |

### `submission/README.md`

| Current line(s) | Change | Rationale |
| ---: | --- | --- |
| 64-70 | Updated documented build sequence to one BibTeX pass. | Keeps submission build docs consistent with `compile.sh`. |
| 107, 111-115, 120-127 | Updated clean/output/layout descriptions to remove active Methods/Supp bibliography claims while noting stale auxiliaries are removed. | Keeps submission docs consistent with the single-reference-list conversion. |
| 169-176 | Rewrote bibliography section to describe a single numbered reference list using `\cite{...}` throughout. | A2 duplicate-reference cleanup. |
| 226-236 removed | Removed obsolete multibib troubleshooting and `\S`/multibib warning. | No active `multibib` setup remains in `paper.tex`. |

### `submission/BUILD_LOG.md`

| Current line(s) | Change | Rationale |
| ---: | --- | --- |
| 5-25 | Added current 2026-06-17 build attempt and local `geometry.sty` blocker. | The task required compiling from `submission` if edits were made; this records the command and environment-level failure. |

## Validation

Commands run from repository root unless noted:

```sh
rg -n "citeMeth|bibliographyMeth|bibliographystyleMeth|newcites|multibib|Methods References|bibtex +Meth|bibtex +Supp" submission/paper.tex submission/compile.sh submission/README.md
```

Result: only the README note that `make clean` removes stale `Meth.*` and
`Supp.*` auxiliaries remains.

```sh
rg -n "/moosefs|05_hic_validation|01_pipeline|scripts/(cladistics|hic|hic-realign|mouse|pedigree)|NATURE_DRAFT|SURVEY_[0-9]|\\.md.*L[0-9]" submission/paper.tex submission/bibliography.bib
```

Result: no matches.

```sh
rg -n "233 individuals|From the 233|465 HPRC|466 near|233 HPRC|233 samples|AFR 67|k = 75|k = 50|k = 14|Leiden \\$k|\\bXTR\\b" submission/paper.tex
```

Result: no matches.

Build command run from `submission/`:

```sh
make clean && make
```

Result: failed before manuscript content was processed because the local TeX
installation is missing `geometry.sty`:

```text
! LaTeX Error: File `geometry.sty' not found.
```

Available local TeX binaries were checked (`pdflatex`, `xelatex`, `lualatex`,
`latexmk`), and `kpsewhich geometry.sty` did not find the missing package.

## Unresolved Author Decisions

- `submission/paper.tex` still has TODO-only `Author Contributions` and
  `Competing Interests` sections. A5 supplied possible text, but the author
  contributions require author confirmation and the competing-interests
  declaration must be true for all authors. No author-judgment statement was
  invented in this mechanical pass.
- The active manuscript text has no `XTR` occurrence, so no acronym edit was
  needed. If Fig. 1 is regenerated with a visible `XTR` label, expand the first
  visible occurrence as `X-transposed region (XTR)` or use `X-transposed region`
  as the visible label.
- The public source-data/accession table referenced by Data and Code
  Availability is not created in this pass. This report only removes the private
  `/moosefs` roots from the active manuscript and points to public releases,
  accessions, citations, and source-data tables at the appropriate level.
