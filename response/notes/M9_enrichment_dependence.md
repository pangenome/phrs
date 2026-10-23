# M9 — Referee #1 major point 9: enrichment driven by invisible gene content and one amplified family

> Provenance for this note: see `response/PROVENANCE.md` (tool versions: `response/notes/tool_versions_and_pins.md`).

Status: COMPLETE
Repo: /home/erikg/phrs (branch: none — no git operations performed)
Started: (see journal below)

## Journal

- [step 1] Created notes skeleton.

## Task summary

(1) Verify the referee's arithmetic: 402 gene copies with midpoints in CHM13 PHRs; 215 unannotated
    and excluded; K = 187 across 31,235 term combinations.
(2) Family dependence: leave-largest-family-out recomputation per significantly enriched term;
    test the ~25-fold "spermatid development" claim resting on 19 SEPTIN14-associated duplicon copies.
(3) Report invisible (unannotated) PHR gene content fraction and give a corrected manuscript sentence.

## Inputs (to be verified present)

- AUTHORITY: PEER_REVIEW.md (read); response/POINT_BY_POINT.md (item I own)
- FORBIDDEN: paper_prep/synthesis/SIMULATED_PEER_REVIEW_v1.md (never read)
- submission/supporting_material/ED_Fig2_arm_functional_components.tsv
- gene_family_copy_counting_analysis.md
- data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz
- data/chm13.phrs.bed

## Findings

See "Verdict in one paragraph" and sections (1), (2), (3) below.

## Per-term table

See section 2.2 (top 25 + the named headline terms) and the complete 209-term table in
`scripts/cladistics/enrichment_family_dependence.tsv` (840 data rows = 209 midpoint source-family
+ 209 midpoint locus-family + 211 overlap source-family + 211 overlap locus-family).

## Outputs

- scripts/cladistics/enrichment_family_dependence.py
- scripts/cladistics/enrichment_family_dependence.tsv
- response/notes/M9_enrichment_dependence.md

---

## Verdict in one paragraph

Referee #1's arithmetic is **correct in every digit**: 402 PHR-midpoint gene copies, 215 of them
excluded, K = 187, 31,235 hypotheses, 209 supported terms, 19 spermatid copies at 25.4-fold
(the ledger's own 25.424). One phrase is wrong (`215 are unannotated` — all 402 carry their own
annotation identity; the 215 are copies with no *defensible ontology source*), and one phrase is
wrong in the other direction (`a single SEPTIN14-associated duplicon` — they are 19 copies at
**9 distinct loci on 10 different chromosomes**, i.e. one *family*, not one duplicon; the
submitted manuscript already says "19 SEPTIN14-associated PHR copies" and does not say
"duplicon"). The substantive allegation is **confirmed and stronger than stated**:
the 187-copy PHR annotation universe is drawn from only **39 functional-source families**, of
which **eight families supply 146/187 = 78.1%** of the copies (the eight that dominate the terms
supply 142/187 = 75.9%); **167/187 = 89.3%** of those copies
rest on *inherited* annotation (NCBI Gene "Related functional gene", tier B) rather than on their
own exact annotation; **124 of the 209 supported terms have 100% of their PHR copies from one
family**; and after removing each term's single largest contributing family, **0 of 209** survive
the frozen primary rule (BH <= 0.05 and global Holm <= 0.05), with only 30 even at nominal
p <= 0.05 (median p inflation **14.2 orders of magnitude**). `spermatid development` (GO:0007286)
goes from p = 2.49e-20 to **p = 1.0** because all 19 of its copies are one family. The interpretation
that survives is the one the Discussion already makes — *annotation burden across physical copies of
a few amplified families* — not functional enrichment of the PHR compartment.

---

## (1) The referee's arithmetic: verified / corrected

Everything below is recomputed by `scripts/cladistics/enrichment_family_dependence.py`, which
aborts on any mismatch. Verbatim gate lines from the final run (`/tmp/m9_final.log`):

```
  independent recount, gene midpoints inside PHR BED: 402 (referee: 402)
  independent recount, any overlap with PHR BED:     412
  VERIFIED: 402 raw PHR copies - 215 excluded = 187 retained (K); excluded = 215 (referee: 215)
  VERIFIED: 31,235 hypotheses, 209 supported terms (referee: 209)
  VERIFIED: all 209 supporting-table rows agree with the frozen midpoint rows (a,b,c,d,K,non-PHR,N,p)
```

| referee quantity | value | how verified | verdict |
|---|---|---|---|
| 402 gene copies with midpoints in CHM13 PHRs | 402 | **independent recount**: every `gene` record of `data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz` (61,312 records) tested against the 37 intervals of `data/chm13.phrs.bed`, midpoint = floor((start0+end0)/2), 0-based half-open containment | **VERIFIED** (and equals the pinned ledger `MAPPING_COVERAGE.tsv: phr_raw_physical_copy_cn`) |
| 215 unannotated / excluded | 215 | ledger `phr_ontology_ineligible_cn`; and independently 402-215 = 187 matches the 187 distinct PHR `copy_id`s in the contributor ledger | **VERIFIED as a count; CORRECTED as a description** (see below) |
| K = 187 | 187 | (i) ledger `phr_ontology_eligible_cn`; (ii) distinct PHR copies in `EXACT_TERM_CONTRIBUTORS.tsv.gz` midpoint cells = 187; (iii) constant `K_phr_total` over all 31,235 midpoint rows | **VERIFIED** |
| 31,235 term combinations | 31,235 | row count of `ontology_v7/TERM_RESULTS.tsv.gz` per assignment regime (midpoint and overlap each 31,235); matches `FINAL_V7_REPORT.md` | **VERIFIED** |
| 209 supported terms | 209 | recomputed BH within collection + global Holm **from the frozen p-values alone**, reproducing the frozen `bh_q_within_collection` and `holm_p_global` columns with max abs deviation **0.000e+00**, then counting `q<=0.05 & holm<=0.05` = 209 | **VERIFIED** |
| ~25-fold spermatid development, 19 copies | fold 25.424, a_T = 19 of K = 187 | frozen row GO:0007286/direct | **VERIFIED** |

Two corrections to the referee's wording, both material:

1. **"215 are unannotated" is not right.** `MAPPING_COVERAGE.tsv` reports
   `own_annotation_identity_cn = 402` for the PHR midpoint scope: all 402 copies have their own
   stable annotation identity. The 215 fail *source assignment* — genome-wide the ineligible
   29,346 copies break down (per `GENOMEWIDE_SOURCE_MAP_REPORT.md`) as
   `UNSUPPORTED_FAIL_CLOSED 28,407`, `TYPE_ONLY 890`, `AMBIGUOUS_FAIL_CLOSED 45`,
   `UNRESOLVED 4`. Per-copy PHR-side dispositions for the ineligible 215 are not retained in the
   repository, so the PHR split of those four categories cannot be quoted.
   What *is* exact: of the 215, **111 (51.6%) carry a `LOC*` placeholder symbol** versus 1 of 187
   (0.5%) among the eligible, and the 215 are biotypically lncRNA-heavy
   (lncRNA 103 / 47.9%, pseudogene 54 / 25.1%, miRNA 49 / 22.8%, protein_coding 7 / 3.3%).
   So "functionally unannotated" is defensible for about half of them; "unannotated" as a blanket
   statement is not.
2. **"a single SEPTIN14-associated duplicon" is not right.** All 19 copies map to the same
   functional source (`NCBIGene:346288`, symbol SEPTIN14, disposition
   `EXPLICIT_RELATED_FUNCTIONAL_GENE`, tier B), but they are 9 distinct `gene_name` loci
   (`SEPTIN14P22` x10, `SEPTIN14P13` x2, `SEPTIN14P2/P3/P5/P11/P16/P18/P19` x1 each) on
   10 chromosomes/seqids (chr1, chr2, chr3, chr5, chr7, chr9, chr11, chr16, chr19, chr20).
   The claim that stands is "one family, 19 copies"; the claim that fails is "one duplicon".

---

## (2) Family dependence of the enrichment

### 2.0 What "family" means here, and why it matters more than anything else in this note

Two definitions are carried through everything below, because the answer to the referee depends
entirely on which one is used:

* **source family** (primary) — `functional_source_id`, i.e. the family root through which the
  ontology annotation is actually inherited. This is the unit that the referee's phrase
  "a few amplified families" means, and the unit that makes copies redundant.
* **locus family** (sensitivity) — `gene_name`, the RefSeq Liftoff locus symbol
  (`DUX4L6`, `SEPTIN14P22`, ...). A locus is *not* a family: two loci of the same family are
  near-identical paralogues and are not independent.

**This instability is itself a result and it is large**: at source-family level **0/209** terms
survive leave-largest-family-out under the frozen primary rule; at locus-family level **152/209**
survive. The locus-level number looks reassuring and must not be used in the reply, because
removing `DUX4L6` (25 copies) leaves 40 other DUX4-family copies untouched. The instability is
reported, not resolved in our favour.

### 2.1 The composition of the PHR annotation universe (midpoint, K = 187)

Verbatim from the run:

```
  midpoint PHR universe by source family (39 distinct):
    NCBIGene:100288687    65  (34.8% of 187)     DUX4
    NCBIGene:346288       19  (10.2% of 187)     SEPTIN14
    NCBIGene:54838        14  ( 7.5% of 187)     WBP1L
    NCBIGene:6147         11  ( 5.9% of 187)     RPL23A
    NCBIGene:100287171    10  ( 5.3% of 187)     WASHC1
    NCBIGene:1663         10  ( 5.3% of 187)     DDX11
    NCBIGene:23152         9  ( 4.8% of 187)     CIC
    NCBIGene:2969          8  ( 4.3% of 187)     GTF2I
    NCBIGene:3581          4  ( 2.1% of 187)     IL9R      ... 39 families in total
    annotation confidence (disposition, tier letter) in the midpoint PHR universe:
        {('EXPLICIT_RELATED_FUNCTIONAL_GENE', 'B'): 167, ('EXACT_SELF', 'A'): 20}
```

* eight source families = **146/187 = 78.1%** of the entire PHR annotation universe (copy-count
  ranking, which includes GTF2I at 8 copies; GTF2I dominates no supported term -- see 2.2);
* DUX4 alone = **65/187 = 34.8%**;
* **167/187 = 89.3%** of the copies rest on *inherited* annotation, not their own.
* All 65 DUX4 copies are `pseudogene` (64) or `protein_coding` (1); all 19 SEPTIN14 copies and
  all 14 WBP1L copies are `pseudogene`. The enrichment universe is overwhelmingly
  pseudogene copies carrying inherited annotation.
* Referee's "near-identical members": **not tested here** (no sequence-identity measurement in
  this deliverable; that is the duplicated-segment/identity lane's job). The *amplified-family*
  half of his sentence is fully confirmed; the *near-identical* half is neither confirmed nor
  contradicted by this analysis.

### 2.2 Per-term table (source families, midpoint regime) — top 25 of 209 by frozen global Holm

Full 209-term table (and the overlap regime and the locus regime) in
`scripts/cladistics/enrichment_family_dependence.tsv`, 840 data rows, columns include
`copies_in_phr_a`, `largest_family`, `largest_family_copies_in_term`,
`largest_family_share_of_term_copies`, `p_before`, `p_after_leaveout`,
`bh_q_after_leaveout`, `holm_after_leaveout`, `survives_*`, `term_copies_exact_self_annotation`.

| term_id | term | a/K | largest family | copies from it | share of a | p before | p after leave-family-out | verdict |
|---|---|---|---|---|---|---|---|---|
| GO:0070317 | negative regulation of G0 to G1 transition | 65/187 | DUX4 | 65 | 1.00 | 2.87e-136 | 1.00e+00 | DIES |
| GO:0070316 | regulation of G0 to G1 transition | 65/187 | DUX4 | 65 | 1.00 | 1.70e-133 | 1.00e+00 | DIES |
| R-HSA-9819196 | Zygotic genome activation (ZGA) | 65/187 | DUX4 | 65 | 1.00 | 2.35e-118 | 1.00e+00 | DIES |
| R-HSA-9816359 | Maternal to zygotic transition (MZT) | 65/187 | DUX4 | 65 | 1.00 | 2.27e-81 | 1.00e+00 | DIES |
| GO:0000977 | RNA pol II transcription regulatory region binding | 74/187 | DUX4 | 65 | 0.88 | 8.45e-79 | 4.28e-04 | survives nominal only |
| GO:0001216 | DNA-binding transcription activator activity | 74/187 | DUX4 | 65 | 0.88 | 1.51e-76 | 6.94e-04 | survives nominal only |
| GO:0001228 | DNA-binding transcription activator act., RNA pol II | 74/187 | DUX4 | 65 | 0.88 | 1.16e-75 | 8.38e-04 | survives nominal only |
| GO:0031965 | nuclear membrane | 65/187 | DUX4 | 65 | 1.00 | 4.95e-75 | 1.00e+00 | DIES |
| GO:0005635 | nuclear envelope | 66/187 | DUX4 | 65 | 0.98 | 4.70e-74 | 8.25e-01 | DIES |
| GO:0010948 | negative regulation of cell cycle process | 65/187 | DUX4 | 65 | 1.00 | 3.38e-69 | 1.00e+00 | DIES |
| GO:0006915 | apoptotic process | 65/187 | DUX4 | 65 | 1.00 | 2.15e-68 | 1.00e+00 | DIES |
| GO:0045786 | negative regulation of cell cycle | 65/187 | DUX4 | 65 | 1.00 | 4.67e-68 | 1.00e+00 | DIES |
| GO:0000976 | transcription cis-regulatory region binding | 65/187 | DUX4 | 65 | 1.00 | 1.54e-66 | 1.00e+00 | DIES |
| GO:0045893 | positive regulation of DNA-templated transcription | 84/187 | DUX4 | 65 | 0.77 | 2.68e-62 | 2.15e-06 | **only flag that keeps BH** (q=0.041) |
| GO:0003700 | DNA-binding transcription factor activity | 75/187 | DUX4 | 65 | 0.87 | 3.49e-61 | 7.00e-03 | survives nominal only |
| GO:0008285 | negative regulation of cell population proliferation | 65/187 | DUX4 | 65 | 1.00 | 7.16e-61 | 1.00e+00 | DIES |
| GO:0010564 | regulation of cell cycle process | 75/187 | DUX4 | 65 | 0.87 | 1.39e-57 | 1.46e-02 | survives nominal only |
| GO:1990837 | sequence-specific double-stranded DNA binding | 66/187 | DUX4 | 65 | 0.98 | 2.28e-55 | 9.64e-01 | DIES |
| GO:0000981 | DNA-binding transcription factor act., RNA pol II | 83/187 | DUX4 | 65 | 0.78 | 6.09e-55 | 8.60e-05 | survives nominal only |
| GO:0051726 | regulation of cell cycle | 75/187 | DUX4 | 65 | 0.87 | 1.10e-52 | 3.72e-02 | survives nominal only |
| GO:0043565 | sequence-specific DNA binding | 84/187 | DUX4 | 65 | 0.77 | 3.07e-52 | 1.15e-04 | survives nominal only |
| GO:0012501 | programmed cell death | 65/187 | DUX4 | 65 | 1.00 | 2.45e-45 | 1.00e+00 | DIES |
| GO:0008219 | cell death | 65/187 | DUX4 | 65 | 1.00 | 9.06e-45 | 1.00e+00 | DIES |
| GO:0045944 | positive regulation of transcription by RNA pol II | 74/187 | DUX4 | 65 | 0.88 | 2.02e-44 | 2.15e-01 | DIES |
| GO:0006355 | regulation of DNA-templated transcription | 93/187 | DUX4 | 65 | 0.70 | 8.43e-43 | 2.88e-05 | survives nominal only |

Headline terms the manuscript actually names, with the same columns:

| term_id | term | a | dominant family | copies | share | fold before | fold after | p before | p after | tier-A copies in term |
|---|---|---|---|---|---|---|---|---|---|---|
| GO:0070317 | negative regulation of G0 to G1 transition | 65 | DUX4 | 65 | 1.00 | 789 | 0 | 2.87e-136 | 1.00e+00 | 1 |
| R-HSA-9819196 | Zygotic genome activation (ZGA) | 65 | DUX4 | 65 | 1.00 | 221 | 0 | 2.35e-118 | 1.00e+00 | 1 |
| GO:0031965 | nuclear membrane | 65 | DUX4 | 65 | 1.00 | 31.7 | 0 | 4.95e-75 | 1.00e+00 | 1 |
| GO:0045893 | positive regulation of DNA-templated transcription | 84 | DUX4 | 65 | 0.77 | 10.6 | 3.67 | 2.68e-62 | 2.15e-06 | 2 |
| GO:0038146 | CXCL12 signaling | 14 | WBP1L | 14 | 1.00 | 397 | 0 | 1.26e-27 | 1.00e+00 | 0 |
| GO:1990700 | nucleolar chromatin organization | 10 | DDX11 | 10 | 1.00 | 850 | 0 | 2.40e-21 | 1.00e+00 | 0 |
| **GO:0007286** | **spermatid development** | **19** | **SEPTIN14** | **19** | **1.00** | **25.4** | **0** | **2.49e-20** | **1.00e+00** | **0** |
| GO:0001764 | neuron migration | 19 | SEPTIN14 | 19 | 1.00 | 24.3 | 0 | 5.47e-20 | 1.00e+00 | 0 |
| GO:0005515 | protein binding | 169 | DUX4 | 65 | 0.38 | 1.28 | 1.2 | 8.73e-11 | 1.66e-04 | 9 |
| GO:0004919 | interleukin-9 receptor activity | 4 | IL9R | 4 | 1.00 | 227 | 0 | 3.92e-08 | 1.00e+00 | 1 |

Summary statistics (midpoint / source families, 209 terms):

```
  share==1.0 count: 124          (a_after==0 count: 124)
  p<=0.05 after: 30              bh after: 2              primary rule after: 0
  median p_before 4.495e-15      median p_after 1.000e+00
  median log10 drop: 14.24
```

Which family dominates how many of the 209 terms, and how many copy-edges it supplies
(copy-edges = sum of a_T over the terms in question; total over all 209 = 8,340):

| source family | PHR copies | terms it dominates | copy-edges | share of all copy-edges |
|---|---|---|---|---|
| DUX4 | 65 | 72 | 4,680 | 56.1% |
| DDX11 | 10 | 49 | 490 | 5.9% |
| SEPTIN14 | 19 | 37 | 703 | 8.4% |
| WASHC1 | 10 | 28 | 280 | 3.4% |
| WBP1L | 14 | 10 | 140 | 1.7% |
| IL9R | 4 | 6 | 24 | 0.3% |
| CIC | 9 | 5 | 45 | 0.5% |
| RPL23A | 11 | 2 | 22 | 0.3% |
| **these 8 total** | **142** | **209** | **6,384** | **76.5%** |

Note the two different "top eights", which are easy to conflate and must not be: the eight largest
families *by PHR copy count* are DUX4 65, SEPTIN14 19, WBP1L 14, RPL23A 11, WASHC1 10, DDX11 10,
CIC 9, GTF2I 8 = 146/187 (78.1%); the eight families that *dominate* the 209 terms are the ones in
the table above and hold 142/187 (75.9%). GTF2I, though the 8th largest by copy count,
**dominates no supported term**, and removing it costs nothing (209/209 survive).

Referee's sentence "The rest are near-identical members of a few amplified families, so one
family drives a term" is therefore correct term by term: **every one of the 209 supported terms
has its copies dominated by exactly one of eight families, 124 of them exclusively.**

### 2.3 Late-breaking: 89.3% of the universe is inherited annotation

```
  annotation confidence in the midpoint PHR universe:
    {('EXPLICIT_RELATED_FUNCTIONAL_GENE','B'): 167, ('EXACT_SELF','A'): 20}
```

84/209 terms contain **no** exact-self (tier A) copy at all; the remaining 125 contain between 1
and 9. The term `spermatid development` contains **0** tier-A copies: all 19 are tier-B inheritance
from SEPTIN14. This is a second, independent reason the enrichment cannot be read as evidence
about PHR gene function: the signal is inherited annotation propagated from 20 reference copies
onto 167 paralogues.

### 2.4 The referee's specific spermatid claim — tested exactly

```
  GO:0007286 spermatid development: a_T=19 / K=187 ; c_T=127 / non-PHR=31779 ; fold=25.424 ;
      p=2.493319e-20 ; holm=7.772e-16
  its 19 PHR copies: functional sources {'NCBIGene:346288': 19}
  symbols {'SEPTIN14': 19}
  loci (gene_name) {'SEPTIN14P22': 10, 'SEPTIN14P11': 1, 'SEPTIN14P16': 1, 'SEPTIN14P19': 1,
                    'SEPTIN14P13': 2, 'SEPTIN14P2': 1, 'SEPTIN14P3': 1, 'SEPTIN14P5': 1,
                    'SEPTIN14P18': 1}
  chromosomes/seqids {'chr11': 2, 'chr16': 2, 'chr19': 2, 'chr1': 2, 'chr20': 1, 'chr2': 2,
                      'chr3': 2, 'chr5': 2, 'chr7': 2, 'chr9': 2}
  distinct genomic loci: 9  -> the referee's 'single duplicon' is not literally correct
                              (9 loci, one family)
```

* 19 copies: **confirmed**. 25.4-fold: **confirmed** (ledger 25.424, referee "~25").
* one family (`SEPTIN14`/`NCBIGene:346288`), 100% of the copies: **confirmed**.
* removing that family removes **all 19** copies, so the test is degenerate: a' = 0, K' = 168 and
  the exact upper tail is **p = 1.0**. The term does not "weaken"; it ceases to exist.
* in the coherent re-analysis that removes SEPTIN14 from *every* term, 172/209 survive, i.e.
  SEPTIN14 removal costs exactly the 37 terms it dominates and nothing else (see 2.5);
  `spermatid development` is one of those 37.

### 2.5 Coherent re-analyses: remove one family, or the k largest, from *every* term

The per-term test above is a dependence diagnostic (each term loses its own largest family, so the
resulting hypothesis set is heterogeneous). The referee also asks "how many of the 209 terms
survive". That is answered by a coherent re-analysis: remove the same families from every term and
re-run BH + global Holm over all 31,235 tests each time. Verbatim (midpoint / source):

```
  single:NCBIGene:100288687 (DUX4)     : 137/209 supported terms survive
  single:NCBIGene:346288    (SEPTIN14) : 172/209
  single:NCBIGene:54838     (WBP1L)    : 199/209
  single:NCBIGene:6147      (RPL23A)   : 204/209
  single:NCBIGene:100287171 (WASHC1)   : 179/209
  single:NCBIGene:1663      (DDX11)    : 155/209
  single:NCBIGene:23152     (CIC)      : 204/209
  single:NCBIGene:2969      (GTF2I)    : 209/209
  cumulative_top1 (DUX4)             : 137/209
  cumulative_top2 (+SEPTIN14)        : 106/209
  cumulative_top3 (+WBP1L)           :  99/209
  cumulative_top4 (+RPL23A)          : 100/209   <-- NOT monotone, see 2.6
  cumulative_top5 (+WASHC1)          :  68/209
  cumulative_top6 (+DDX11)           :  13/209
  cumulative_top7 (+CIC)             :   8/209
  cumulative_top8 (+GTF2I)           :   6/209
```

Reading: DUX4 alone costs 72 terms (137 survive) and DUX4 + SEPTIN14 cost 103 (106 survive);
removing the eight largest families (146 of the 187 copies) leaves **6 of 209**. The single-family
lines also show the converse: GTF2I removal costs nothing (209/209), while CIC and RPL23A each cost
5 (204/209) although CIC dominates 5 terms and RPL23A only 2. The overlap regime gives
139 / 105 / 98 / 100 / 67 / 11 / 8 / 5 -- the same shape.

### 2.6 Instability of the conclusion under regime and definition (reported as a finding)

1. **Family definition dominates every other choice.** source-family removal: 0/209 survive the
   primary rule (30 nominal). locus-family removal: 152/209 survive (209 nominal, 198 BH).
   The referee's claim is about families and is confirmed; anybody quoting the locus number to
   rebut him would be answering a different question.
2. **Assignment regime barely matters.** midpoint (K=187, 209 supported) vs any-overlap
   (K=193, 211 supported): 0 vs 0 survive the primary rule; 30 vs 33 nominal; single-family
   removal 137 vs 139 (DUX4), 172 vs 172 (SEPTIN14); cumulative top-8 6 vs 5. The dependence is
   not a PHR-boundary artefact.
3. **Cumulative removal is not monotone** (99 -> 100 -> 68 as the ladder adds RPL23A, then
   WASHC1, in the midpoint regime's deterministic order). Removing a family changes both the numerator and the denominator K, and a
   small family can slightly *help* the remaining terms by shrinking the PHR arm. Any statement of
   the form "survivors decrease monotonically as families are removed" would be false; the honest
   statement is "6 of 209 survive removal of the eight largest families (5 in the overlap
   regime)".
4. **Cumulative order is fixed deterministically** by (-PHR copy count, family id), so the two
   runs with different hash seeds give byte-identical output.

### 2.7 What this analysis does NOT do (and why)

* The fully symmetric leave-family-out (family's copies removed from the **non-PHR background**
  as well) is **not computable from retained artefacts**. The non-PHR per-copy side is not in the
  repository: `EXACT_TERM_CONTRIBUTORS.tsv.gz` contains 16,763 midpoint PHR edges + 816
  overlap-only PHR edges and **no** non-PHR contributor rows; the genome-wide table
  (`PHYSICAL_COPY_TERM_EDGES.tsv.gz`, 2,929,709 closure edges per
  `ontology_v6/GENOMEWIDE_SOURCE_MAP_REPORT.md`) is not distributed, and the SHA-pinned source tables in
  `ontology_v6/INPUT_MANIFEST.tsv` are absent (`.../chm13_copy_enrichment/sources/` does not
  exist). So the background (c_T, M_T, N) is held at its observed, family-inclusive value.
  Direction of the resulting bias, stated as an argument and not as a computed fact: holding the
  removed family's non-PHR term burden in M keeps the expected PHR count at the higher value, so
  p_after is, if anything, *larger* than a symmetric removal would give; the counts in 2.2/2.5 are
  therefore strict rather than generous. (The companion N adjustment is negligible: N = 31,966
  and the largest family contributes 65 copies.)
* The referee's closing request — repeat against identity- and length-matched non-subtelomeric
  duplications, and with copies collapsed to families — is **not this deliverable** (it is the
  matched-control lane). What is supplied here that partially serves it: the coherent
  family-removal ladder in 2.5, which is the "collapsed to families" question asked from the other
  end (how much of the signal survives when families are removed rather than collapsed).
* No new annotation was needed: the pinned repository artefacts listed in the inputs section were
  sufficient; nothing outside the repository was read (no moosefs access at all).

---

## (3) The invisible fraction and what it does to the interpretation

Verbatim:

```
  PHR gene copies (midpoints): 402 ; with a defensible ontology source: 187 ; without: 215 (53.48%)
  genome-wide gene copies: 61312 ; eligible 31966 ; ineligible 29346 (47.86%)
  biotypes, all 402 PHR copies: pseudogene=204 (50.7%), lncRNA=104 (25.9%), miRNA=51 (12.7%),
      protein_coding=24 (6.0%), transcribed_pseudogene=18 (4.5%), misc_RNA=1 (0.2%)
  biotypes, 187 ontology-eligible: pseudogene=136 (72.7%), transcribed_pseudogene=17 (9.1%),
      protein_coding=17 (9.1%), miRNA=2 (1.1%), lncRNA=1 (0.5%)
  biotypes, 215 ontology-ineligible: lncRNA=103 (47.9%), pseudogene=54 (25.1%), miRNA=49 (22.8%),
      protein_coding=7 (3.3%), misc_RNA=1 (0.5%), transcribed_pseudogene=1 (0.5%)
  PHR bases: 6014981 ; gene-body bases covered by the 187 eligible copies: 727114 (12.1%) ;
      by the 215 ineligible copies: 1207003 (20.1%)
  of the 215 ineligible copies, 111 (51.6%) carry a LOC* placeholder symbol
  of the 187 eligible copies, 1 (0.5%) carry a LOC* placeholder symbol
```

* **53.5% of PHR gene copies are invisible** to the enrichment analysis (215/402). By gene-body
  bases the invisible share is **62.4%** (1,207,003 / 1,934,117 bp of gene body; the two sets are
  merged within themselves and are not mutually exclusive, so this is an upper-ish summary, not a
  partition of PHR bases).
* The genome-wide invisible fraction is 47.9% (29,346/61,312), so the PHR deficit is real but
  modest (**+5.6 points**, not a qualitative difference). This matters: the referee's stronger
  reading ("over half the PHR gene content is invisible") is true in absolute terms but is not
  evidence that PHRs are uniquely unscannable — 48% of the whole genome's gene copies are in the
  same position.
* What it does to the interpretation: the denominator of every fold is the annotated genome, and
  the numerator is the annotated 46.5% of PHR copies, which is *not a random 46.5%* — it is
  enriched 72.7% pseudogene (versus 25.1% among the invisible) and dominated by eight families.
  The enrichment therefore cannot be read as "PHRs are functionally specialised"; it can only be
  read as "the copies that happen to inherit annotation are disproportionately from a few
  amplified families". That is exactly the Discussion's existing framing
  (`paper.tex:405-410`: "describe annotation burden across physical gene copies rather than
  independent functional systems"), and it is the framing the Abstract and the enrichment
  subsection currently lose.

### Corrected manuscript sentences (drop-in replacement text)

For the enrichment subsection (`paper.tex`, after "rather than collapsed by gene name or family"):

> "Of the 402 PHR copies only 187 (46.5%) carry a defensible GO or Reactome source, and
> 167 of those 187 (89.3%) inherit that annotation from a related reference gene rather than
> carrying it themselves; the 187 copies come from 39 functional-source families, with eight
> families contributing 146 (78.1%) and DUX4 alone contributing 65 (34.8%). Each of the 209
> supported terms is dominated by a single one of these families (124 of 209 exclusively) and no
> term keeps its significance once its largest contributing family is removed (median p inflation
> 14 orders of magnitude); in a coherent re-analysis, removing the eight largest families leaves
> 6 of the 209 terms. The result is therefore a statement about the annotation burden of a few
> amplified families, not about the functional content of the PHR compartment."

For the Abstract (replaces "copy-number enrichment is concentrated in gene modules associated with
chromosome dynamics, intracellular organization and germ-cell development"):

> "copy-number enrichment of inherited annotation is concentrated in a small number of amplified
> families — chiefly DUX4, SEPTIN14, DDX11, WASHC1 and WBP1L — whose annotated copies fall into
> modules associated with chromosome dynamics, intracellular organization and germ-cell
> development."

---

## Headline claims: keep / reword / drop

| claim in the manuscript | verdict after this analysis | action |
|---|---|---|
| "209 significantly enriched terms" (`paper.tex:224`) | arithmetically exact (verified 209 at BH<=0.05 & Holm<=0.05) | **KEEP the number**, add the family-dependence caveat of (3) |
| "spermatid development ... enriched 25-fold ... carried by 19 SEPTIN14-associated PHR copies" (`paper.tex:237-239`) | **exactly right**, and the manuscript never said "duplicon" | **KEEP verbatim**; optionally add "all 19 are inherited annotations from SEPTIN14 and removal of that family removes the term entirely" |
| "the PHR compartment is functionally enriched for germ-cell development / chromosome dynamics" (abstract + subsection framing) | 0/209 survive leave-largest-family-out under the primary rule; 89.3% of copies are inherited annotation | **REWORD/DROP the functional reading** in Abstract and subsection; the Discussion already has the correct framing |
| "Compared with all annotated gene copies elsewhere in the genome" (`paper.tex:221-223`) | only 46.5% of PHR copies are in the comparison; 53.5% are invisible | **REWORD** to "compared with the 31,966 genome-wide copies carrying a defensible GO/Reactome source" |
| "distributed across multiple communities" for the WASH/DDX11 module (`paper.tex:232-236`) | DDX11 dominates 49 terms (removing it alone: 155/209 survive, i.e. 54 lost) and WASHC1 28 (179/209, i.e. 30 lost); they are ladder steps 5 and 6 | **KEEP the community distribution** (it is a separate, factual claim) but do not present it as independent systems |
| "A second distributed system is associated with male germ-cell development" | one family (SEPTIN14), 37 terms, all tier-B | **REWORD**: "a single amplified family (SEPTIN14) supplies the germ-cell, septin and cytokinesis terms" |
| Extended Data Fig. 2 "six functional classes" | classes are post-inference display summaries of terms dominated by ≤8 families | **KEEP as display**, add one sentence that the six classes are 8 family clusters re-described |

Claims that must **not** be made in the reply:
* "the enrichment is robust to removing the dominant family" — false at source-family level.
* "152/209 terms survive" — true only for locus-symbol removal, which is not family removal.
* "the PHR is uniquely annotation-poor" — 47.9% of genome-wide copies are equally invisible.

---

## Reproduction

```
python3 scripts/cladistics/enrichment_family_dependence.py --self-test      # exact-tail oracle
python3 scripts/cladistics/enrichment_family_dependence.py                  # writes the TSV
python3 scripts/cladistics/enrichment_family_dependence.py --output -       # TSV to stdout
```

Self-test output (tail): exhaustive sweep `N<=12, K<=4, M, a` all match a brute-force enumeration
oracle; upper tail monotone over a=0..K at N=31,966/K=187/M=146; all 12 named cases pass.
(Incidental cross-check: the monotonicity case uses M=146, which happens to be the spermatid row's
M = a_T + c_T = 19 + 127, and the sweep reproduces the frozen p = 2.493319e-20 at a=19.)
Full run: 840 data rows; walls recorded in the journal (98.2 s, 108.6 s in-script; the box ran at
load average 38-73 out of 256 cores due to other lanes, so the wall is an upper bound).

Design decision worth flagging to whoever integrates this: `ED_Fig2_arm_functional_components.tsv`
is a 6-class x 48-end copy matrix and contains no term-level or family-level information, so it
cannot support this analysis; the arithmetic lives in `ED_Fig2_supported_ontology_terms.tsv` plus
the pinned `ontology_v7` ledger. No new annotation was needed.

---

## Journal (appended as work proceeded)

- [step 2] Read PEER_REVIEW.md line 88 (verbatim point 9); response/POINT_BY_POINT.md:134-137
  (R1-M9 "rework around the 215 unannotated copies; quantify the SEPTIN14 dependence of the
  25-fold spermatid term; report enrichment conditional on the largest family"). NEVER opened
  SIMULATED_PEER_REVIEW_v1.md.
- [step 3] Located inputs. All present and readable (no moosefs reads needed):
  ED_Fig2_arm_functional_components.tsv (30,189 B), ED_Fig2_supported_ontology_terms.tsv
  (136,569 B, 209 rows), ontology_v7/{TERM_RESULTS.tsv.gz, EXACT_TERM_CONTRIBUTORS.tsv.gz,
  MAPPING_COVERAGE.tsv, NAMED_COHORT_AUDIT.tsv, FINAL_V7_REPORT.md},
  data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz, data/chm13.phrs.bed.
  MAIN FINDING about inputs: `ED_Fig2_arm_functional_components.tsv` (6 classes x 48 ends) carries
  NO term-level or family-level information; it cannot support this analysis. The term-level
  arithmetic lives in `ED_Fig2_supported_ontology_terms.tsv` + the pinned ontology_v7 ledger.
  The per-copy family attribution lives ONLY in `ontology_v7/EXACT_TERM_CONTRIBUTORS.tsv.gz`.
- [step 4] Pre-check with awk/python: independent GFF3+BED recount of gene midpoints in PHRs.
  First attempt returned 408 for any-overlap (wrong: early loop break keyed on midpoint order);
  fixed to a full interval sweep -> 412, matching MAPPING_COVERAGE.tsv exactly.
- [step 5] Wrote scripts/cladistics/enrichment_family_dependence.py. Self-test (`--self-test`)
  initially failed 90 exhaustive tiny-scale cases (K=0 with a>0). Root cause: guard ordering in
  the exact upper tail. Fixed; all cases now match the brute-force oracle.
- [step 6] First full run aborted loudly twice, both correct aborts:
  (a) "ledger overlap PHR 412 != my recount 408" (the loop-break bug above);
  (b) "overlap recomputed support count 211 != frozen 0" -- the overlap regime stores its verdict
      in `sensitivity_support`, not `primary_support`. Fixed per-regime support column.
- [step 7] Second look at the annotation-confidence columns: my first implementation compared
  `mapping_confidence` (values "A"/"B") to "EXACT_SELF" and therefore reported "209/209 terms are
  100% inherited". That was a BUG (silently all-False). Corrected to use
  `source_assignment_disposition` ("EXACT_SELF" / "EXPLICIT_RELATED_FUNCTIONAL_GENE"); the
  corrected answer is 125/209 terms contain at least one tier-A (exact-self) copy, and this was
  cross-checked by an independent standalone recount (125).
- [step 8] Determinism: PHR-universe family counting iterated a `set`, so Counter tie order was
  hash-seed dependent. Fixed by explicit `(-count, family)` sorting and a `tops()` helper for all
  "most common" printouts. Verified: two runs with PYTHONHASHSEED=1 and =7 produce byte-identical
  TSV and identical logs.
- [step 9] Walls (this box is shared; load averages recorded at run start):
  19:56:30 loadavg 38.86 57.83 47.02 -> 98.2 s total wall (full run).
  20:02:24 loadavg 46.59 52.45 47.95 -> 2 x full run (seeds 1, 7), deterministic.
  20:06:58 loadavg 73.11 65.48 54.55 -> 108.6 s in-script wall (final production run).
  No other job of mine on the box; the load is other lanes' work (256 cores, 1.0 TB RAM).
  Self-test (`--self-test`) is a few seconds and does no I/O beyond arithmetic.
- [step 10] Cross-checks that caught errors in my OWN intermediate output (kept here for audit):
  (a) an early version of the "which family dominates" summary listed GTF2I as dominating 8 terms.
      Re-deriving that table from the TSV showed GTF2I dominates **0** supported terms (it is the
      8th largest family by copy count, 8 copies, but its terms are dominated by other families).
      The table and every sentence built on it were corrected against the TSV, not against memory.
      Lesson recorded: the "8 largest families by copy count" (146/187 copies) and "the 8 families
      that dominate the 209 terms" (142/187 copies) are DIFFERENT sets; both are now stated.
  (b) a first pass compared `mapping_confidence` ("A"/"B") against "EXACT_SELF" and silently
      produced "209/209 terms 100% inherited". Corrected and then cross-checked against an
      independent standalone recount (125 terms contain >=1 tier-A copy, 84 contain none).
  (c) an early run piped through `head` died of SIGPIPE before writing the TSV, so a query read a
      stale file and reported 0 tier-A copies. The stale-read was caught by checking the TSV mtime
      against the run start time; all final numbers come from the completed run at 20:08:47.
- [step 11] Final gates, verbatim from the completed run (`/tmp/m9_final.log`):
  exhaustive exact-tail oracle PASS; independent PHR recount 402/412 == pinned ledger;
  31,235 x 2 regimes; 209 = 211 frozen support counts reproduced by my own BH/Holm;
  all 209 supporting-table rows agree with the frozen rows; contributor ledger reproduces a_T for
  every one of 31,235 hypotheses in both regimes (16,763 and 17,579 PHR term edges);
  my exact-tail implementation matches the frozen `p_exact_upper` to max relative deviation
  1.608e-10 (midpoint) / 1.902e-10 (overlap) over 1,456 / 1,664 rows with a_T > 0.
  Output: `scripts/cladistics/enrichment_family_dependence.tsv`, 33 columns x 840 rows,
  0 empty/nan cells, byte-identical across PYTHONHASHSEED=1 and =7.
- [step 12] Idempotence check before hand-off: a fresh full run (`--output /tmp/m9_verify.tsv`)
  at 20:17:21 (loadavg 121.89 115.06 85.92, in-script wall 109 s) reproduced the shipped
  `scripts/cladistics/enrichment_family_dependence.tsv` **byte for byte** (`diff -q` clean).
  Only these three files were created; no existing file was modified, and no git command was run.
