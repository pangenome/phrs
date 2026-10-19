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

## Audit result: PAN028 label exposure (pending item 3, partially done)

`all_pedigrees_patches.tsv`, WashU dataset = 538 patches:

| dataset | patches | comment |
|---|---|---|
| `PAN028 maternal (hap1) vs PAN027` | **310 (58%)** | KID-side; hap1 is maternal only on chr2,6,8,10,11,13,16,17,19,20 |
| `PAN027 maternal (hap1) vs PAN010` | 167 | SELF side; labels consistent (`PAN027_mat`) |
| `PAN027 paternal (hap2) vs PAN011` | 61 | SELF side; labels consistent (`PAN027_pat`) |

Of the 310 KID patches, **228 (74%) sit on chromosomes where hap1 is the PATERNAL
copy**, i.e. they are calls about PARTNER's unsequenced chromosome and cannot support
a transmission claim; only 82 are label-correct. Wrong-label counts by chromosome:
chr15 94, chr21 38, chr14 25, chr22 18, chr5 15, chr9 10, chr1 10, chr3 7, chr7 7, chr4 4.

`crossover_like` patches: 18 total, 15 in KID. **7 of those 15 are wrong-label**
(chr1, chr3, chr4, chr5, chr22) — including the two the drafts quote explicitly:
chr3 q 366,221-368,310 → chr9q:2 (2,089 bp: the call already retracted in rev 10) and
chr3 q 262,953-290,922 → chr7p:2 (27,969 bp, the "largest" event). The 8 label-correct
crossover-like calls are all chr13p → chr15p/chr21p/chr22p, i.e. acrocentric p-arm sharing.

Consequences:
- The submitted manuscript (`submission/paper.tex`) contains **no** 494/538 or
  "13 of 16 crossover-like" claim — checked: clean. So the submitted pedigree text
  needed only the corrections already made (items 1-7).
- The claim still stands in working documents and must be fixed before reuse:
  `paper_prep/figures/fig4/caption.md` line 9 ("494/538 = 92 %"),
  `paper_prep/synthesis/NATURE_DRAFT_v4.md` line 48 and
  `NATURE_DRAFT_v5.md` line 53 ("Thirteen of sixteen crossover-like events are in
  PAN028, confirming that meiotic-resolution inter-chromosomal breakpoints transmit
  across generations"; "largest 27.97 kb on PAN028 maternal chr3q").
- The 92% within-community figure is a property of query/ref arm membership, so it is
  not automatically invalidated, but the 228 wrong-label calls must be re-run on the
  correct maternal copy before the count is reported as evidence of transmission.

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
3. **PAN028 label exposure (audited, see above)**: re-run the untangle for the 228
   wrong-label KID patches with the per-chromosome maternal copy; fix the 92% and
   "13 of 16 crossover-like" wording in `paper_prep/figures/fig4/caption.md`,
   `NATURE_DRAFT_v4.md`, `NATURE_DRAFT_v5.md`; check `scripts/pedigree/*.py` and
   report §14 for the same assumption.
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
