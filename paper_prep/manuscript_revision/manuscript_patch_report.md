# Manuscript Patch Report

Date: 2026-06-17

Task: `manuscript-revision-paper-patch`

## Summary

Applied the guarded final revision package to the active manuscript where the
upstream records supplied evidence-backed wording. No author-only title,
declaration, accession, or optional source-data decisions were invented.

## Changed Files

- `submission/paper.tex`
  - Replaced the abstract's "linked q-arm group" and zygotene-peak wording with
    the conservative two-tier continuum and broad contact-support framing.
  - Reframed the community Results around locally dense blocks embedded in a
    broader subtelomeric similarity continuum.
  - Replaced UPGMA/tree wording with average-linkage display/order language.
  - Softened q-arm language to "enriched q-arm neighborhood."
  - Reordered and softened the 3D Results: human pointwise PHR contacts are
    aggregate coordinate-level measurements, flanking unique-sequence analyses
    are the anti-artifact control, and mouse meiotic Hi-C supports a broad
    prophase-I association rather than a resolved zygotene-specific maximum.
  - Reframed the pedigree Results and Fig. 5 caption as candidate patch calls
    nested within three assayed WashU child-haplotype transmissions; community
    labels are assigned after patch calling.
  - Updated Fig. 2, Fig. 4, Fig. 5 and Extended Data Fig. 1 captions to remove
    load-bearing p-values, "caught in the act", phylogeny/clade framing, and
    zygotene-specific overclaiming.
  - Updated Methods language for Leiden resolution scans, UPGMA comparison,
    exploratory NJ/bootstrap, MAPQ0/random-primary handling, mouse stage
    interpretation, pedigree untangle, and cM/Mb limitation provenance.
- `submission/BUILD_LOG.md`
  - Added the 2026-06-17 post-patch compile attempt and the persistent
    `geometry.sty` environment blocker.
- `paper_prep/manuscript_revision/manuscript_patch_report.md`
  - This report.

## Files Considered But Not Changed

- `submission/bibliography.bib`: no revision-package item required a
  bibliography edit. Exact public accessions/source-data rows were not present
  in the inspected artifacts, so no accession metadata was invented.

## Commands Run

```
wg msg read manuscript-revision-paper-patch --agent $WG_AGENT_ID
wg show manuscript-revision-paper-patch
git status --short
git log --oneline main..HEAD --decorate
sed -n '1,260p' paper_prep/manuscript_revision/final_revision_package.md
sed -n '1,260p' paper_prep/manuscript_revision/G_abstract_intro_title_package.md
sed -n '1,560p' paper_prep/manuscript_revision/CD_continuum_community_synthesis.md
sed -n '1,620p' paper_prep/manuscript_revision/BF_3d_contact_synthesis.md
sed -n '1,620p' paper_prep/manuscript_revision/E_pedigree_audit.md
nl -ba submission/paper.tex | sed -n '1,720p'
rg -n 'discrete blocks|phylogeny|clade|clades|closed class|q-arm sextet|tight|linked q-arm|previously uncharacterized|peak at zygotene|peaks at zygotene|caught in the act|same signal increases|pointwise Spearman.*p|Mann-Whitney|strict-MAPQ|strict MAPQ|random-ligation|within-Leiden filter|/moosefs|end-to-end-report|SURVEY_|XTR|citeMeth|fig:ed[2-9]' submission/paper.tex submission/bibliography.bib
perl -0pi -e 's/\t//g; s/closed q-arm class/hard q-arm partition/g' submission/paper.tex
cd submission && bash compile.sh
git diff --name-only
git diff --stat
```

## Build Result

Command:

```
cd submission
bash compile.sh
```

Result: failed with exit code 1 during the first `pdflatex` pass before
manuscript content was processed:

```
! LaTeX Error: File `geometry.sty' not found.
```

Interpretation: this is the same local TeX environment/toolchain blocker
reported by the upstream mechanical-fixes fan-in. The attempted compile did not
reach manuscript syntax, citation, cross-reference, or figure processing.

## Targeted Search Result

The post-patch targeted search returned no hits for the requested stale or
unsupported manuscript patterns:

```
discrete blocks
phylogeny
clade / clades
closed class
q-arm sextet
tight
linked q-arm
previously uncharacterized
peak at zygotene / peaks at zygotene
caught in the act
same signal increases
pointwise Spearman ... p
Mann-Whitney
strict-MAPQ / strict MAPQ
random-ligation
within-Leiden filter
/moosefs
end-to-end-report
SURVEY_
XTR
citeMeth
fig:ed2..fig:ed9
```

## Unresolved Blockers and TODOs

- TODO(author decision): title remains unchanged. The final package treats title
  strength as author-sensitive; no explicit downstream author decision record
  was found.
- TODO(author decision): author contributions and competing interests remain
  comment placeholders in `submission/paper.tex`; these require author-confirmed
  truth statements.
- TODO(author/source data): exact public accessions/source-data rows for HG002
  Pore-C/Hi-C/CiFi, Dip-C, sperm scHi-C and mouse meiotic Hi-C were not present
  in the inspected artifacts. The manuscript keeps public-citation/source-data
  wording but does not invent accession identifiers.
- TODO(author decision): optional placement of F_ST and cM/Mb as background or
  honest-null material remains an author-level framing choice. The active
  manuscript only keeps the cM/Mb limitation with corrected provenance.
- TODO(compute): no completed formal sample-resampling stability analysis exists
  for Leiden assignments. The Methods now state that resolution/UPGMA scans are
  fixed-matrix checks, not sample-resampling stability.
- TODO(environment): local compile validation remains blocked by missing
  `geometry.sty`.
