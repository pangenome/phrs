# Reviewer-response change log

Manuscript text is **corrected**, never phrased as a reply to reviewers. This log
tracks what changed and why, so the response letter can be written from it.

Naming: `PEDIGREE_NAMING.md` (MOM=PAN010, DAD=PAN011, SELF=PAN027, KID=PAN028,
PARTNER=KID's father, unsampled).

## Decisions

- One ectopic event is retained: SELF.p.chr9q terminus derived from DAD.h2.chr3q.
- The chr1p/chr5q signal is ancestral homology, not an event.
- KID-side candidates are withdrawn: PARTNER is unsampled, and those scans had used
  the paternal chromosome copy anyway.
- Fig5 is simplified to the ectopic event: A pedigree, B identity-crossover zoom,
  C paternal scan, and a new **D: pairwise alignments of the sequences involved**
  (SELF.p.chr9q vs DAD.h2.chr3q / DAD.h2.chr9q / DAD.h1.chr9q). Panels D/E
  (maternal scan, grandchild scan) retired. 1p/5q relabelled ancestral; 21p/4p
  dropped.
- Population-frequency support for the standing arrangement: NOT pursued.

## Done

| # | change | where |
|---|---|---|
| 1 | Results: one documented exception, not "two putative autosomal exchanges"; chr1p/chr5q stated as ancestral homology | `paper.tex` §Ectopic subtelomeric recombination within a pedigree |
| 2 | Results: maternal-scan paragraph trimmed; granddaughter candidates removed | same section |
| 3 | Results: closing paragraph now says the family transmits an ectopic product rather than a de novo exchange with resolved breakpoints | same section |
| 4 | Fig5: panels D/E removed; legend 21p/4p dropped, 1p/5q relabelled ancestral; caption rewritten, title now "Pedigree view of an ectopic subtelomeric exchange" | `paper.tex` Fig5 block |
| 5 | Methods: "Inference rule" paragraph (best parental match = source; exceptions PAR1, acrocentric arms, this case) | `paper.tex` §Pedigree direct-alignment windows |
| 6 | Methods: "Haplotype labelling" paragraph (per-chromosome parent of origin; maternal copy resolved by identity; unsampled PARTNER ⇒ KID paternal scans not reported) | same |
| 7 | Methods: Fig5C--E references reduced to C | same |
| 8 | Evidence record with all numbers, retractions and revisions rev 1-12 | `REVIEWS.txt` |
| 9 | Canonical naming + per-chromosome maternal table | `PEDIGREE_NAMING.md` |
| 10 | New **Fig5D** pairwise-alignment panel built and added to the figure block + caption; asset, box coordinates and generator script committed | `submission/fig/MainFigures/Fig5D_pairwise_9q3q.{pdf,png}`, `Fig5D_highlights.tsv`, `submission/scripts/figures/make_fig5d_pairwise.sh` |

## Pending

1. **Rebuild Fig5B/Fig5C**: `Fig5B_zoom9q3q.pdf` and `Fig5C_paternal.png` need
   regenerating with the new labels (SELF.p, DAD.h1/h2); `Fig5A_pedigree.pdf`
   should print the nickname/copy scheme. Fig5D (new) is already built.
2. **Retire or supplement the D/E assets**: `Fig5D_maternal.png`,
   `Fig5E_pan028.png`, and the scripts that produce them
   (`submission/scripts/figures/make_fig5_panels.sh`, `make_fig5_zoom.R`,
   `make_fig5_pedigree_supporting_material.py`) — move to supporting material or
   delete, and drop their rows from
   `submission/supporting_material/Fig5_pedigree_supporting_material.tsv`.
3. **Audit other PAN028-dependent analyses** for the same label error
   (`scripts/pedigree/untangle_multimap_tracts.py`, `untangle_native_merge_tracts.py`,
   `fig5_followup_transmission_check.py`, `patch_tract_*`, report §14) — each needs
   its chromosome copy resolved before reuse.
4. **Acrocentric paragraph**: state the flagger breakdown for the acrocentric arms
   specifically (global: 66/538 patches overlap problematic regions, 14.2% of patch
   bp vs 6.7% flank background, 0/133 gene-conversion tracts) so the error-prone
   status is explicit.
5. **Margin table**: per-event identity margin + alternative-haplotype exclusion for
   every retained call, for Methods/supplement.
6. **Build verification**: LaTeX is not installed on this machine; `paper.tex` changes
   are unverified by compilation. Run `bash compile.sh` on a machine with
   pdflatex+bibtex and check `grep -c undefined paper.log`.
7. **Discussion/abstract audit** for "recurrent exchange" phrasing that implies
   multiple pedigree events; keep the population-level recurrent-exchange claim
   (supported by the pangenome survey) separate from the single germline event.
