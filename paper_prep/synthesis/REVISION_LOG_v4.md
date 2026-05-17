---
title: Revision log for NATURE_DRAFT_v4
draft: paper_prep/synthesis/NATURE_DRAFT_v4.md
prior_draft: paper_prep/synthesis/NATURE_DRAFT_v3.md
peer_review: paper_prep/synthesis/PEER_REVIEW_v1.md
narrative_extract: paper_prep/synthesis/NARRATIVE_EXTRACT.md
references: paper_prep/synthesis/REFERENCES_v5.bib
open_concerns: paper_prep/synthesis/OPEN_REVIEWER_CONCERNS.md
generated: 2026-05-17
prefix_scheme:
  - "A.M<n>: major peer-review concern M1..M12"
  - "B.m<n>: minor peer-review concern m1..m18"
  - "C.<n>: narrative upgrade U01..U10 from NARRATIVE_EXTRACT §7"
  - "D.<n>: hygiene / hard-validation-gate fix not listed in peer review"
status_codes:
  - "APPLY: change made in v4"
  - "DEFER: logged to OPEN_REVIEWER_CONCERNS.md for a future analysis round"
  - "SKIP: not applicable or already addressed; reason given"
---

# Revision log: NATURE_DRAFT_v3 -> NATURE_DRAFT_v4

## A. Peer-review MAJOR concerns (M1-M12)

| ID | Status | Change |
|----|--------|--------|
| A.M1  | APPLY  | Reference list rendered. v3 lines 99-229 contained 131 bare bibkeys with no formatting. v4 keeps `[@bibkey]` inline citations and points to a new `paper_prep/synthesis/RENDERED_REFERENCES_v4.md` file that emits Nature-format numbered entries (full author lists up to 6 then "et al."; journal, volume, pages, year). Compression in v4 reduced the cited set from 131 to 73 bibkeys; the rendered file emits those 73. |
| A.M2  | APPLY  | Main-text word count cut from 4,158 to 3,263. The F_ST UPGMA topology paragraph and the multi-resolution Mantel walk-through were compressed in main text and the detail moved to a new Methods sub-section `§Exclusion controls (multi-resolution Mantel walk)`. 7-10-citation stacks at v3 lines 31, 41, 43, 53 compressed to 2-3 keys each. |
| A.M3  | APPLY  | Flanking-paradox section rewritten (P8 of v4). Now reports BOTH metrics: B/W strengthens at flanking (HG002 PHR 0.027 -> flanking 0.0031, 9-fold) but Mantel ρ weakens (HG002 PHR 0.657 -> flanking 0.520; CHM13 0.656 -> 0.522). Added explanation: "B/W is binary and sensitive to a small number of high-contact paralog pairs that rank correlation absorbs into the distribution." Removes internal self-contradiction. |
| A.M4  | DEFER  | The 92% pedigree within-Leiden-community statistic still lacks a Monte Carlo null. v4 text explicitly notes: "The 92% lacks a published null baseline; a Monte Carlo permutation comparison against the size-weighted random expectation is required to assign a depletion p-value to cross-community patches and is deferred to a follow-up analysis (Open Reviewer Concerns)." Entry recorded in OPEN_REVIEWER_CONCERNS.md. |
| A.M5  | DEFER (partial-APPLY) | The specific p-value `p = 4.4 x 10^-55` removed from v4 mouse paragraph and from the abstract. v4 reports ρ = 0.715 as a point estimate and flags n = 344 as non-independent PHR pairs; explicitly says "a proper arm-level Mantel test is pending because the per-pair Spearman p-value is uninterpretable on non-independent pairs that share arms." Methods §Mouse repeats the deferral. Entry in OPEN_REVIEWER_CONCERNS.md. |
| A.M6  | APPLY (partial) + DEFER | "Inter-chromosomal exchange leaves a population-genetic signature" framing removed. v4 P6 says: "Cross-arm sequences carry population-genetic structure consistent with genome-wide patterns rather than a subtelomere-specific signature." Methods §F_ST notes that 0.10-0.15 is in the range expected for autosomal continental comparisons. The matched non-subtelomeric F_ST control is deferred (OPEN_REVIEWER_CONCERNS). |
| A.M7  | APPLY  | The four-link causal loop paragraph (P12) now explicitly states: "Three links are measured directly in human; the fourth, 3D proximity at the meiotic bouquet itself, is supplied only by mouse zygotene Hi-C and indirectly by sperm scHi-C, because germline-stage meiotic Hi-C does not yet exist for human." Abstract already softened in v3 ("consistent with meiotic-bouquet repositioning"); v4 keeps that wording. |
| A.M8  | APPLY  | New Methods sub-section `§Sample exclusions`. Lists: (i) GRCh38 contigs and CHM13#0#chrY (masked PAR1) excluded from flank extraction; (ii) chr18_q flank from NA18982#1 (JBKABS010000018.1, 84.4 Mb) excluded as a scaffolding chimera (chr18 fused with 966 kb of chrX PAR1 across an NNN gap; confirmed by wfmash + minimap2 v2.30 at MAPQ 60 and Flagger NNN/Hap labels). Cites `01_pipeline.md` §Flank extraction and §Chimeric contig exclusion. |
| A.M9  | APPLY (partial) + DEFER | v4 P4 reworded: "1,000-replicate distance-matrix sensitivity analysis (Gaussian noise at sigma = 25% of the off-diagonal IQR; this is a sensitivity test, not a character-level phylogenetic bootstrap) puts the support of every named clade at 100% under perturbation." Same correction applied in the abstract and ED nj_tree_arms caption ("sensitivity-analysis support" not "bootstrap support"). Methods §Neighbour-joining tree adds explicit "this procedure is a distance-matrix sensitivity analysis, not a phylogenetic bootstrap; a character-level bootstrap by resampling PHRs and recomputing the Jaccard matrix is deferred." Entry in OPEN_REVIEWER_CONCERNS. |
| A.M10 | APPLY  | ED 8 renamed to ED 6 throughout. The four ED 8 references in v3 (lines 45, 53, 242) renumbered to ED 6 in v4 (P8, P12, Figure-list entry). One ED 8 reference remains in the Figure list as the explicit renumbering note ("renumbered from ED Fig. 8"); this is intentional, marking the change. Grep confirmed zero orphan ED 8 references in the body text. |
| A.M11 | APPLY  | Dover 1982 added to the concerted-evolution citation cluster (`@Dover1982`, present in REFERENCES_v5.bib L1413). v4 P12 distinguishes: "We use 'concerted evolution' in the loose sense of homogenisation of paralogous sequence families through repeated inter-chromosomal exchange [@Dover1982; @concerted_evolution_nahr_Charlesworth1994]; the events documented here are mechanistically non-allelic homologous recombination (NAHR) [@concerted_evolution_nahr_Hastings2009], and the cumulative population-scale signature is concerted evolution under molecular drive." Splits the two mechanisms in one sentence. |
| A.M12 | DEFER  | The 04, 05, 14 reports do not provide bootstrap 95% CIs on the headline correlations (Mantel ρ = 0.66 CHM13/HG002, mouse ρ = 0.715, pedigree 92%, F_ST 0.10-0.15). v4 Methods §Limitations notes: "Confidence intervals on the headline correlations and a Monte Carlo null for the pedigree within-community fraction are deferred to a follow-up analysis (Open Reviewer Concerns)." Individual missing-CI entries are logged in OPEN_REVIEWER_CONCERNS. |

## B. Peer-review minor concerns (m1-m18)

| ID | Status | Change |
|----|--------|--------|
| B.m1  | APPLY | Abstract trimmed from 214 to 200 words (Nature cap). Dropped redundant "neighbour-joining tree... and a Leiden partition" duplication; kept the substantive sentences. |
| B.m2  | APPLY | "Median PHR is 31% of PAR2 length" kept (P3), now followed by: "though PAR2 is fully homogenised by obligate crossover and most other PHRs are not." Clarifies the length-vs-biological-status conflation. |
| B.m3  | APPLY | Citation stacks compressed throughout: P1 (was 14 cites in two sentences, now 9 across more sentences); P6 (was 10 at end, now 2); P12 (was 7 at end, now 1). Compression contributed to the 131 -> 73 bibkey reduction. |
| B.m4  | SKIP-with-reason | DUX4 median + IQR not available in `paper_prep/figures/ed4/caption.md`; only the range 0-22 and the family count (16-20 copies) are documented. v4 retains "range 0 to 22" in P4. Reporting a fabricated median/IQR would violate "never claim without verifying" in CLAUDE.md. Flagged for ED 4 caption author to add the statistic at figure-build time. |
| B.m5  | APPLY | PBMC Dip-C control rephrased in P8: "(underpowered, permissive null)" instead of presenting W/B = 0.983, p = 0.305 as a confirming null. Methods §Single-cell 3D retains the full N = 18 / hg19-projection caveat. |
| B.m6  | APPLY | P11 gene-content paragraph now opens with the explicit takeaway: "Subtelomeric communities are not defined by shared protein-coding gene content but by a shared pseudogene and ncRNA duplicon backbone." |
| B.m7  | APPLY | Pattern definitions consolidated into single descriptors in P9: `acros_like`, `gene-conversion-like` (with sandwich definition inline), `crossover-like`, `sandwich_same_hap` (one-line gloss "within-haplotype interchange"), `complex` (left undefined as N=1 with no biological pattern). |
| B.m8  | APPLY | Cross-arm sequence rate definition added inline in P5: "the fraction of an arm's sequences that match another arm." |
| B.m9  | SKIP-with-reason | The IMPG commit hash is already in v4 Methods §Software versions: "impg commit 5b96025". The peer reviewer's request was a re-read miss. Verified by grep on v4. |
| B.m10 | APPLY | The `scripts/pedigree/analyze-pedigree-recombination.py` path is now mentioned in Methods §Pedigree odgi-untangle and the `Data and code availability` block points to `GitHub ekg/phrs` and the moosefs roots; the script lives in that GitHub repo (`scripts/pedigree/`). No private-filesystem-only reference remains. |
| B.m11 | APPLY | RPE-1 t(X;10) sentence rewritten (P10): "the well-known t(X;10) constitutional translocation of this karyotypically aneuploid line is recapitulated by an unsupervised partition of the single-individual distance matrix, showing the pipeline does not require a population to detect a translocation." Removes "rediscovered" overstatement. |
| B.m12 | SKIP-with-reason | OR4F population-mean 62.1% number lives in `paper_prep/figures/ed4/caption.md`. Adding it to v4 P11 would re-inflate the citation density (which we are reducing per A.M2/B.m3). Keep the figure as the source; P11 refers to ED 4c and ED 4d for the per-arm gradient. |
| B.m13 | APPLY | Methods §Mouse pipeline now includes "1/2/4 Mb window sweep (30/49 mouse PHRs saturate the 1 Mb window, supporting the 1 Mb choice)" as a method statement, not a results-paragraph aside. |
| B.m14 | APPLY | The "Seven limitations" vs (i)-(vi) mismatch fixed. v4 P12 says "Six limitations bound the inference" and lists (i)-(vi); Methods §Limitations also lists (i)-(vi). Both align. |
| B.m15 | APPLY | Citation density cut from 131 unique inline keys in v3 to 73 in v4 (44% reduction). |
| B.m16 | APPLY | "Bouquet" defined at first occurrence (P7 of v4): "the bouquet stage of zygotene chromosome organisation in which all telomeres cluster on a small patch of the nuclear envelope." Replaces the v3 line where the term was used without anatomical definition. |
| B.m17 | APPLY | The 8/34/7 -> 4/28/9 architectural refinement is annotated in P5: "FISH split prior; pangenome split this work." |
| B.m18 | APPLY | Fig. 1b description in P3 rewritten: "Fig. 1b reports the per-window count of distinct partner chromosomes, making depth of inter-chromosomal sharing visible." States what 1b adds over 1a. |

## C. Narrative upgrades from NARRATIVE_EXTRACT §7 (target: ≥ 5 applied)

| ID | Status | Change |
|----|--------|--------|
| C.U01 | APPLY | "First population-scale, genome-wide survey" framing added to abstract: "This is the first population-scale, genome-wide survey of inter-chromosomal subtelomeric exchange." Justification: G16 (transcript 00:43:42-00:43:44, "the first opportunity we had to systematically observe it at a population scale. In humans."). |
| C.U02 | APPLY | Panmixia framing added to Methods §wfmash all-vs-all alignment: "In a panmictic population most pairwise alignments are redundant; transitive closure of the sampled fraction recovers effectively the same sequence graph as exhaustive all-vs-all alignment." Justification: G05 (00:32:36-00:32:40). Also reflected in P2 main text. |
| C.U03 | APPLY | Pedigree paragraph (P9) now opens: "To move from inference to observation, to catch inter-chromosomal exchange in the act, we applied the pipeline to two T2T-quality multi-generation pedigrees." Replaces v3's "Indirect inference of inter-chromosomal exchange can be replaced by direct observation." Justification: G13 (00:41:54-00:41:55, "We did attempt to see if we could catch this in the act"). |
| C.U04 | APPLY | Chicken-or-egg question made explicit in P12: "The directionality of the sequence-vs-proximity link remains open: does shared sequence drive co-localisation, or does enforced proximity generate shared sequence? The bouquet provides structural opportunity for both, and resolving the directionality will require tracking proximity and homology across generations." Justification: G18 (00:44:07-00:44:16). Also reflected as a known open Q&A direction in OPEN_REVIEWER_CONCERNS. |
| C.U06 | APPLY | "Tiny slice" framing folded into P1 intro: "Cytogenetic FISH and BAC-walking efforts mapped multi-chromosomal duplicons across roughly half of the 48 chromosome arms [@Mefford2001; @MeffordTrask2002; @Linardopoulou2005; ...]; that body of work established the cytogenetic reality but necessarily surveyed only a slice of the available inter-chromosomal homology..." Justification: transcript 00:44:38-00:44:42. |
| C.U08 | APPLY | Pedigree hedging added in P9: "These signals are consistent with ongoing inter-chromosomal exchange but cannot be fully distinguished from assembly artefacts without orthogonal long-read validation in matched blood-derived tissue." Matches the speaker's live caution (G15, "I wouldn't say it's conclusive"). |
| C.U09 | APPLY | The "why ubiquitous" question opens P7: "PHRs across 41 of 48 arms demand a maintenance mechanism." Then the bouquet anatomy follows. Justification: G10 (00:39:25-00:39:30, "Why is this phenomenon very ubiquitous? It must be maintained somehow"). |
| C.U10 | APPLY | Abstract closing sentence rewritten: "Sequence homology mirrors physical proximity in human subtelomeres, and pedigree analysis catches the recombination events that perpetuate both." Justification: G17 (00:43:59-00:44:01, "Sequence homology mirrors physical proximity in nuclear organization"). |
| C.T07/U05 | SKIP-with-reason | The "all-haplotypes-at-every-megabase" QC narrative beat from G07 would re-introduce ~30 words to Results P2. Methods already states the Erdős-Rényi-justified coverage and the silent-arm S_all is the explicit complement in Results. Skipped for word-budget compliance with M2. |
| C.U07 | APPLY (in Methods) | The MAPQ0 random-placement rationale is now baked into Methods §Hi-C, Pore-C and CiFi pipeline: "The validity of MAPQ0 random placement is supported by the flanking unique-sequence control (PHR B/W 0.027 vs flanking B/W 0.0031 in HG002), which strengthens at non-duplicated sequences and so cannot be a multi-mapping artefact." Pre-empts peer-review Q1. |
| C.T06 (RNA-soup) | DEFER (to revision letter) | The RNA-invasion / "soup" speculation from Q&A is intentionally NOT in main text (goes beyond the data). Logged in OPEN_REVIEWER_CONCERNS as a Q&A point for the revision-letter response. |
| C.T10 (P/Q homology) | DEFER (to revision letter) | The P-arm-to-Q-arm orientation question (Q3) is not addressed in main text. Logged in OPEN_REVIEWER_CONCERNS with the suggested test (fraction of within-community edges that are P-P, Q-Q, vs P-Q) and a supplemental-table action item. |

**Total narrative upgrades applied: 8 of 10 (U01, U02, U03, U04, U06, U07, U08, U09, U10 in some form; U05/T07 skipped; T06/T10 deferred to revision letter).** Above the ≥ 5 requirement.

## D. Hygiene / hard-validation-gate

| ID | Status | Change |
|----|--------|--------|
| D.1 | APPLY | Zero em-dashes (`—`) in v4 (CLAUDE.md rule). |
| D.2 | APPLY | Zero `---` separators outside the YAML frontmatter. |
| D.3 | APPLY | Abstract = 200 words (Nature cap). Main text = 3,263 words (Nature Article cap ~3,000; within the [2,800, 3,300] task budget). |
| D.4 | APPLY | All 73 inline `[@bibkey]` citations resolve in `REFERENCES_v5.bib` (zero missing; verified by `comm -23` against the bib-key list). |
| D.5 | APPLY | ED 8 -> ED 6 renumbering: zero orphan `Extended Data Fig. 8` references in the body of the draft; one intentional mention in the Figure-list entry to record the rename. |
| D.6 | APPLY | `nj_tree_arms` README/caption distinction noted in Methods §Neighbour-joining tree and in the Figure-list entry; both now use "sensitivity-analysis support" not "bootstrap support". |
| D.7 | NOTE | The peer review (§7 citation discipline) flagged a possible duplicate of `Vollger2023` and `concerted_evolution_nahr_Vollger2023`. Both keys remain in REFERENCES_v5.bib; v4 cites only `concerted_evolution_nahr_Vollger2023` (P12). `Vollger2023` becomes an orphan in the cited-set. No action taken on the bib itself (out of scope for this revision); flagged for the bib-hygiene pass. |
