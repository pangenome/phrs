# Stocktake: repo evolution, the real review, and what to splice

Date: 2026-10-23. Prepared after discovering that the response material had been drafted against
`paper_prep/synthesis/SIMULATED_PEER_REVIEW_v1.md` (internal, agent-147, 2026-05-17) instead of the
real review, `PEER_REVIEW.md` (Nature editor + **five referees**).

---

## 1. What the real review actually asks

Referee #2 and #4 are co-reviewers submitting jointly with other reports (no independent text).
The substance is in **#1** (evolutionary genomics/recombination, 10 major + 13 minor),
**#3** (genome organization, 8 major + minors) and **#5** (human genomics, 6 specific).

### Referee #1

| # | ask | our status |
|---|---|---|
| M1 | **Concerted evolution is never tested.** Duplicon phylogenies (D4Z4/DUX4L, TUBB8/B, OR4F, SEPTIN14): do copies cluster by end or by chromosome, more than divergence time predicts? | **not addressed** (our module architecture is suggestive but not this test) |
| M2 | **Pedigree cannot distinguish a new exchange from an inherited haplotype.** Needs the tract absent from *both* paternal copies; three sub-points: (a) chr5q/chr1p called in Fig 5C *and* 5E, (b) all Fig 5C calls on the same paternal haplotype, (c) margin ≈1% with comparable dips outside the tract. Asks for both-parental-haplotype comparison, genome-wide margin distribution, permuted parent-child false-positive rate, unbroken-contig check. | **SUBSTANTIALLY ANSWERED — this session's core work.** Confirms their diagnosis: child-side claims were on the *paternal* copy (−73% of calls), chr5q/chr1p is ancestral, one event survives with heterozygous-donor evidence. **Still missing:** genome-wide margin distribution, permuted false-positive estimate, contig-breakpoint check |
| M3 | **Exchange rate conflicts with Lemmers 2010** (four 4q/10q transfers in all human evolution). | **partially, indirectly** — our conclusion (one ancestral arrangement, transmitted; not repeated de novo events) removes most of the conflict, but no explicit rate calculation |
| M4 | **Multi-mapping reads placed at random.** "The flank control is a good idea but its logic needs stating." | **answered in the response material, NOT in the manuscript** — dilution argument (uniform placement attenuates ~1/k per end, symmetric across classes → PHR enrichment is a lower bound) + the flank numbers. The flanking control itself has been deleted from `paper.tex` |
| M5 | **Generic segmental-duplication control**: non-subtelomeric pericentromeric/interstitial duplications matched on identity, length, copy number. | **not addressed** (roadmap; our Q2b) |
| M6 | **Mouse flanks are acro/telocentric and satellite-dominated**, so "p-arm" flanks sit at both the top of the Jaccard axis and the maximum of the multi-mapping problem; p-arms peak earlier, q-arms at diplotene. | **partially** — our arm-level Mantel and the p/q split reproduce their observation; the **satellite-composition confound is not addressed** |
| M7 | **A null preserving dependency**: permutation shuffling PHR labels holding distance-to-telomere and arm identity, bootstrap over chromosome ends; "Fig. 4C needs an explicit test, since 'strongest at the zygotene bouquet' appears in the abstract." | **ANSWERED** — arm-level Mantel with 10,000 row-and-column permutations; the zygotene claim was removed and replaced by a pachytene peak (Introduction, Results, Fig4c caption, Discussion all made consistent) |
| M8 | **500 kb flank right-censors the length distribution**; re-extract saturated ends at 2–5 Mb, state the censored fraction. | **not addressed** |
| M9 | **Gene enrichment is invisible-content**: 215 of 402 copies unannotated, K = 187, 25-fold spermatid enrichment rests on 19 copies of one SEPTIN14 duplicon. | **not addressed** |
| M10 | **Place the compartment in a comparative frame** (apes, yeast, Phaseolus, Plasmodium VAR). | **not addressed** |
| m2,4,5,6,7,8,11,12,13 | reference-free scope statement; Fig 1A vs 1B disagreement; acrocentric vs general PHRs in Fig 1A; "concentrated at ends" quantified vs a segdup-density null; Fig 4A vertical-axis points and per-assay B/W labelling; **Fig 5 number inconsistencies (PAR1 138 vs 92 kb; 28 kb vs a different length in 5E)**; **population-stratified descriptives**; D4Z4/D4Z4L recent literature; "unorthodox" in the title | mostly **not addressed**; the Fig 5 panels were *deleted* rather than reconciled |

### Referee #3

| # | ask | our status |
|---|---|---|
| M1 | Novelty beyond Linardopoulou 2005 / Guarracino 2023 | framing response needed |
| M2 | Fig 1A "same pattern recurs" is not quantitatively tested | not addressed |
| M3 | **Random assignment of near-identical PHR reads could create or redistribute contacts** — restrict pairs or show it does not | answered in response material (dilution), not in manuscript |
| M4 | **B/W 0.0013 (flanks) vs 0.028 (PHRs) described as strengthening** — the flank paradox | same as R1 M4; needs the manuscript fix |
| M5 | Somatic vs meiotic evidence; human meiotic claim | **partially** — mechanism paragraphs rewritten, but the human-meiosis gap stands |
| M6 | **41/48 arms depend on thresholds** (95% identity, ≥5 alignments, 5 kb, 500 kb flank, 4 failed windows) | not addressed (we did threshold sensitivity for the *pedigree* calls instead) |
| M7 | **Jaccard captures aggregate similarity, not order/position** — distinguish similarity from colinearity | **ANSWERED — and extended into a new result**: module order concordance (not conserved), telomere-relative orientation (98.6% same), module maps (ED3) |
| M8 | **Call pedigree findings candidate exchanges consistently**; not evidence for recurrence | **ANSWERED** — one structurally-supported event, recurrence withdrawn |

### Referee #5

| # | ask | our status |
|---|---|---|
| L86–90 | sensitivity of the 15,668 PHR count to window size | not addressed |
| L152–154 | stability of the copy-number-aware GO results across assemblies | not addressed |
| L189–191 | which Hi-C/Pore-C datasets, which cell types; mouse cell type | partially (Methods exist, clarity fix) |
| L193–195 | **why only 2,830 pairs from 15,668 PHRs** — what filtering? | not addressed |
| L229–230 | Fig 5 legend: what are the dark-grey inter-chromosomal exchanges? | **moot/needs rework** — Fig 5 panels were deleted |
| L348–351 | why sample 11.6% of pairs rather than 3.5% | not addressed |

**Tally.** Of ~25 substantive asks: **4 answered** (R1 M2 largely, R1 M7, R3 M7, R3 M8), **3 answered in the letter but not in the manuscript** (R1/R3 M4 multi-mapping; R1 M6 partially), **~18 not addressed**.

---

## 2. The mouse result: resolved, not wrecked

The referees raised the mouse analysis twice, and we have changed it in exactly the direction
they asked:

- **R1 M7** asked for "a null preserving that dependency: a permutation shuffling PHR labels while
  holding distance-to-telomere and arm identity, with bootstrap intervals from resampling
  chromosome ends", and for the abstract's "strongest at the zygotene bouquet" to be tested.
- **What we did:** arm-level Mantel test with 10,000 row-and-column permutations on arm-level
  matrices (the dependency-preserving test), giving ρ = 0.41 / 0.58 / 0.63 / 0.33 across leptotene
  → diplotene, all p < 1e-4, null centred at 0 ± 0.11. The peak moved from zygotene to **pachytene**,
  so the bouquet claim was removed rather than merely caveated. Consistency was propagated to the
  Introduction, Results, Fig4c caption and Discussion.
- **R1 M6**'s observation ("p-arms peak earlier, q-arms at diplotene") is reproduced in our text:
  p-arm pairs peak at zygotene (ρ 0.285→0.508), q-arm pairs at diplotene (0.352→0.420).

**Residual gap:** R1 M6's deeper warning — that mouse p-arm flanks are satellite-dominated, so the
Jaccard axis and the multi-mapping problem both maximise there — is **not** addressed. That needs a
satellite-masked re-analysis of the mouse flanks, not a wording fix. Verdict: the mouse analysis is
now defensible and answers M7; it does not yet answer M6.

---

## 3. Repo evolution and what a splice would involve

- **56 commits since 2026-10-15**, **14 of them touching `submission/paper.tex`**.
- `submission/paper.tex` went **1026 → 1166 lines**. Net growth hides real deletions: the flanking
  control, RPE-1, Dip-C/PBMC, the pedigree pattern categories (sandwich_same_hap/complex) and the
  numbered limitations list were all removed — **and referees quote several of those** (R3 M4 quotes
  B/W 0.0013/0.028; R5 L229 quotes the Fig 5 legend's grey exchanges; R1 m8 quotes Fig 5 tract
  lengths).
- The version the referees read is `1235d4b` ("as actually submitted. oops!"), with `c23c449`
  (continuous line numbering, which R1 asked for) and `cc7921d`/`ab32b41` (funding) after it.

**Recommended splice strategy — restore the submitted text, re-apply only our corrections:**

Base = `submission/paper.tex` at `cc7921d` (submitted + lineno + funding). Then re-apply as isolated
hunks, in this order:

1. Pedigree section correction — `8123abd`, then the copy-assignment/QC/support-text revisions
   (`298d7b6`, `a9e6edf`, `85e9a9b`, `e6d74fd`).
2. Mouse analysis replacement — `c27433d`, plus the consistency fixes `e5c8be9`, `9a71b24`.
3. Population differentiation subsection — `5349d3e`, `3c70f44`, plus the PCLAI paragraph `b362563`.
4. Statistical-controls Methods paragraph + character-bootstrap consequence — `c27433d`.
5. Fig 5 rebuild (A pedigree names, B zoom, C paternal ribbon, new D pairwise) — `159cfd9` onward.
6. **ED 3 module maps** — `159cfd9` (new, no conflict with the base).
7. **Reinstate, do not delete**: the flanking unique-sequence control with *both* B/W and Mantel
   numbers (R3 M4/R1 M4), the acrocentric-vs-general PHR distinction, the Fig 5 tract-length numbers
   reconciled rather than removed (R1 m8), and the pedigree pattern categories (R1 m7).

That keeps every piece of science we produced while restoring the material referees explicitly cite.

---

## 4. Genuinely new work still required (not fixable by splicing)

Ranked by how directly a referee asked for it:

1. **R1 M1** — duplicon phylogenies to test concerted evolution (D4Z4/DUX4L, TUBB8/B, OR4F, SEPTIN14).
2. **R1 M8** — 2–5 Mb re-extraction for saturated ends; report the censored fraction.
3. **R1 M5 / R3 M3 control** — non-subtelomeric pericentromeric + interstitial duplications, matched
   on identity/length/copy number, through the identical pipeline.
4. **R1 M9** — rework the gene enrichment around K = 187 and the SEPTIN14 dependence.
5. **R1 M2 completions** — genome-wide margin distribution, permuted parent-child false-positive rate,
   contig-breakpoint check.
6. **R1 M6** — satellite-masked mouse re-analysis.
7. **R3 M6** — threshold sensitivity for the 41/48 PHR call.
8. **R5 L193** — explain the 2,830-pair filter; **R5 L86** — window-size sensitivity.
9. **R1 M3** — explicit rate comparison against Lemmers 2010.
10. **R1 m11** — descriptive PHR length/community by superpopulation (cheap, and directly supported by
    the F_ST work we now have).
11. **Editor's compliance list** — Data availability, Code availability (DOI-minting repository),
    ethics statement, ED as `.jpg/.tif/.eps` (not PDF), Source Data per figure, word count in the
    cover letter, tracked changes, updated Reporting Summary and Code/software checklist.
