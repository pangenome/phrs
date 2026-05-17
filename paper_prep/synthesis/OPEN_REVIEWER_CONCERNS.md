---
title: Open reviewer concerns deferred from NATURE_DRAFT_v4
parent_draft: paper_prep/synthesis/NATURE_DRAFT_v4.md
parent_review: paper_prep/synthesis/PEER_REVIEW_v1.md
revision_log: paper_prep/synthesis/REVISION_LOG_v4.md
date: 2026-05-17
purpose: |
  Concerns from PEER_REVIEW_v1.md and NARRATIVE_EXTRACT.md §6 (Q&A) that
  NATURE_DRAFT_v4 cannot close in this revision pass. Each entry: concern
  verbatim or paraphrased, reason deferred, what would close it, rough effort.
---

# Open reviewer concerns

The following major concerns from PEER_REVIEW_v1.md are deferred to a future analysis round. NATURE_DRAFT_v4 acknowledges each in-line (Main text or Methods §Limitations).

## D-M4. Pedigree 92% within-Leiden lacks a null baseline

**Concern (peer review M4 verbatim).** "The '92% within Leiden community' pedigree statistic has no null baseline. §pedigree reports 494/538 = 92% of WashU HQ inter-chromosomal patches sit within a Leiden community (Fig. 4a), and frames this as evidence that the population partition predicts where new recombination is found. But the null is never given. With 15 arm-level communities of unequal size, the expected within-community fraction under random pairing depends on the arm-size distribution: if 70% of flank coverage sits in three large communities (C3, C6, C11, C7), the random null could easily be 50-70%. 92% over 50% is not nothing, but it is not what 92% over uniform-prior reads like."

**Why deferred.** Requires a new Monte Carlo permutation pipeline. Not runnable in a draft-revision pass.

**What would close it.** Monte Carlo permutation of the 538 WashU patch source/target arms to random arm pairs (10,000 reps), preserving the per-arm patch-count marginal distribution. Report:
- observed 92% vs the permutation null mean and 95% CI;
- depletion p-value for cross-community patches;
- per-community within-community enrichment with BH-corrected q;
- repeat for the 11 CEPH1463 cross-assembler features (candidate denominator was the number of parent x chromosome-pair feature calls before within-Leiden filtering; record that number for the same null).

**Effort.** ~1 day for the Monte Carlo. Code path: `scripts/pedigree/analyze-pedigree-recombination.py` + a new permutation wrapper. Output: an updated Fig. 4a panel with the null distribution as a dashed line, plus a Methods §Pedigree null block.

## D-M5. Mouse ρ = 0.715 is a non-independent-pair Spearman, not a Mantel

**Concern (peer review M5 verbatim).** "The mouse ρ = 0.715 is computed on 344 non-independent PHR pairs without distance permutation. ... Per-PHR-pair Jaccard vs zygotene Hi-C contact across 344 pairs is structurally non-independent (pairs share arms, share PHRs, and spatial autocorrelation along chromosomes is severe). Reporting an ordinary Spearman p-value on this is statistically inappropriate."

**Why deferred.** Replacing the per-pair Spearman with a Mantel on the arm-level mouse distance and Hi-C matrices requires regenerating the two matrices in the same row/column ordering and running 10,000 row-and-column permutations. The matrices exist; the wiring does not.

**What would close it.** Arm-level Mantel test on the mouse distance matrix vs the mouse Hi-C matrix (per stage: leptotene, zygotene, pachytene, diplotene), with 10,000 row-and-column permutations. Report ρ and permutation p per stage, plus a bootstrap 95% CI on ρ.

**Effort.** ~2-3 days (mouse arm-level matrices need to be rebuilt; per-stage Hi-C matrices need exporting; Mantel framework already exists for the human side and can be reused).

**v4 mitigations.** Abstract no longer reports p = 4.4 x 10^-55. Main-text P10 keeps ρ = 0.715 as a point estimate, flags n = 344 as non-independent PHR pairs, and explicitly says "a proper arm-level Mantel test is pending." Methods §Mouse pipeline reiterates.

## D-M6. F_ST matched-region control

**Concern (peer review M6 verbatim).** "F_ST 0.10-0.15 between AFR and non-AFR is indistinguishable from background and does not support 'Inter-chromosomal exchange leaves a population-genetic signature'. ... Compute the matched genome-wide F_ST on a control set of non-subtelomeric autosomal regions of equivalent length per superpopulation pair, and report the subtelomeric F_ST as a difference or ratio against that baseline."

**Why deferred.** Requires a new F_ST pipeline run on non-subtelomeric autosomal regions of equivalent length per superpopulation pair. Not in scope for a revision pass.

**What would close it.** For each superpopulation pair (AFR vs each of AMR, EAS, EUR, SAS), draw 18,827 random non-subtelomeric autosomal windows of equivalent length (matched to the per-arm flank-length distribution), compute Hudson F_ST per window, summarise per superpopulation pair, and report subtelomeric F_ST as a difference or ratio against this baseline. Add a panel to Fig. 2 or to ED 2.

**Effort.** ~2 days (window sampling, F_ST recomputation, bootstrap CI per pair).

**v4 mitigations.** "Population-genetic signature" claim removed from main text. P6 says "Cross-arm sequences carry population-genetic structure consistent with genome-wide patterns rather than a subtelomere-specific signature." Methods §F_ST notes that 0.10-0.15 is in the autosomal continental range and explicitly defers the matched control.

## D-M9. Character-level NJ bootstrap

**Concern (peer review M9 verbatim).** "The 100% bootstrap support claim on the NJ tree is not a real bootstrap. ... Either (a) run a true character-level bootstrap by resampling PHRs with replacement and recomputing the Jaccard matrix per replicate, or (b) call the procedure what it is (a distance-matrix sensitivity analysis) and report support values labelled 'sensitivity support' not 'bootstrap support'."

**Why deferred.** v4 applies (b) immediately. (a) is the right phylogenetic answer and is deferred.

**What would close it.** Character-level bootstrap: resample 15,668 PHRs with replacement to create B = 1,000 replicate matrices, recompute the per-replicate 41 x 41 arm-level Jaccard, run `ape::nj()` per replicate, summarise per-clade support over B replicates. Repeat for UPGMA.

**Effort.** ~3 days (pgg-graph re-run per replicate is the bottleneck; can be sped up with cached IMPG closure).

**v4 mitigations.** Abstract, main text P4, and Methods §Neighbour-joining tree now consistently say "sensitivity-analysis support under distance-matrix perturbation"; the term "bootstrap" is removed except where the distinction is being explained.

## D-M12. Missing confidence intervals on headline correlations

**Concern (peer review M12 verbatim).** "Mantel ρ = 0.66 is reported with permutation p-values but no 95% CI. Spearman ρ = 0.715 in mouse is reported with a degenerate p-value and no CI. The pedigree 92% has no CI. The F_ST 0.10-0.15 has no CI. For an N = 6 Hi-C study with strong claims, point estimates without intervals are inadequate."

**Why deferred.** The 04/05/14 source-of-truth reports (`end-to-end-report/report/04_heterogeneity.md`, `05_hic_validation.md`, `14_pedigree_recombination.md`) do not currently emit bootstrap CIs for any of the four target statistics. v4 cannot fabricate them.

**Missing CIs to add (priority order):**

1. **Mantel ρ = 0.66, CHM13 Hi-C and HG002 Hi-C.** Bootstrap CI via 10,000 row-and-column permutations of the underlying 41 x 41 matrices (the permutation distribution provides the bootstrap reference). ~0.5 day.
2. **Mantel ρ = 0.66 -> 0.85 trajectory under exclusion controls.** Same method per exclusion set, all five resolutions. ~1 day total.
3. **Pedigree within-Leiden fraction = 92%.** Wilson-score 95% CI on the 494/538 binomial proportion, plus a Monte Carlo null (D-M4 above). ~0.5 day for the binomial CI; D-M4 effort for the null.
4. **F_ST 0.10-0.15 per superpopulation pair.** Block-jackknife 95% CI over arms. ~1 day.
5. **HG002 mouse zygotene ρ = 0.715.** Bootstrap CI on the per-pair Spearman, plus the arm-level Mantel from D-M5. ~0.5 day.

**Effort total.** ~3 days of analysis time. Adds one Methods sub-section and updates each in-text point estimate.

**v4 mitigations.** Methods §Limitations notes the missing CIs explicitly. Each headline statistic is reported as a point estimate without invented CIs.

## D-PeerQ1. Hi-C MAPQ0 multi-mapping artefact (peer reviewer + audience member, Q1)

**Concern.** The Hi-C pipeline retains MAPQ0 reads and randomly assigns them to one of the equally scoring loci. Reviewer Q1 (NARRATIVE_EXTRACT §6) asks whether the within-community contact signal is inflated by this random assignment in PHR regions, which by definition harbour identical paralogous sequence.

**Why partially addressed in v4.** The flanking unique-sequence control (PHR B/W 0.027 -> flanking 0.0031 in HG002) is the existing falsification of the multi-mapping-artefact hypothesis (P8 and Methods §Hi-C, Pore-C and CiFi pipeline). v4 makes this argument explicit.

**What would also close it.** (i) Repeat B/W computation under a MAPQ-strict filter (drop MAPQ0) and report the ratio of strict vs random-placement B/W; (ii) ChIA-PET or Micro-C cross-validation in the same individual if available.

**Effort.** ~1-2 days for the MAPQ-strict re-run on the existing mcool files.

## D-PeerQ2. RNA-soup speculation (peer reviewer + audience member, Q2)

**Concern.** Q2 in NARRATIVE_EXTRACT §6: the speaker speculates that transcription-generated RNA invading DNA could create a "soup" encouraging physical proximity of similar-sequence regions, potentially explaining proximity outside meiosis. This is a Q&A speculation, not in the data.

**Why not in v4 main text.** Beyond the experimental scope; speculative. Including it in main text would weaken the manuscript's evidence-vs-speculation discipline.

**What to do.** Address in the revision letter (cover letter response). Cite Aguilera-style R-loop biology or Engreitz/Chen-style RNA-tethering literature if the editor pushes for a mechanism paragraph.

**Effort.** ~0.5 day for a paragraph in the revision letter.

## D-PeerQ2b. Segmental-duplication generality (audience-member follow-up to Q2)

**Concern.** "If you look at this across all segmental duplications, not just the sub-telomeres ... be a little more careful about the control of the Hi-C data" (NARRATIVE_EXTRACT 00:47:29-00:47:31). The question is whether the sequence-to-3D correlation is unique to subtelomeric PHRs or extends to all segmental duplications genome-wide.

**Why deferred.** A pangenome-scale SD analysis is a separate paper.

**What would close it.** Compute the equivalent Mantel correlation on all HPRC v2 segmental duplications (~50-100 Mb of non-subtelomeric SDs), per the same Hi-C pipeline. Report whether the subtelomeric ρ = 0.66 is comparable to or stronger than the all-SD ρ.

**Effort.** ~1-2 weeks.

**v4 action.** Address in the revision letter; flag that the comparison is on the roadmap.

## D-PeerQ3. P-arm-to-Q-arm homology asymmetry (audience-member Q3)

**Concern.** Q3 in NARRATIVE_EXTRACT §6: P-arm-to-P-arm and Q-arm-to-Q-arm relationships dominate the communities, but P-to-Q homology also exists at opposite orientations. The asymmetry is not quantified in v4.

**What would close it.** For the 15 arm-level communities, count the fraction of within-community edges that are P-P, Q-Q, and P-Q. Add a supplemental table; report the P-vs-Q-orientation asymmetry per community. If P-Q edges cluster in C1 (4q/10q) or C2 (10p/18p), call that out.

**Effort.** ~0.5 day (the 41 x 41 distance matrix already encodes arm orientation; just regroup the edges).

## D-Bib. REFERENCES_v5.bib hygiene: possible Vollger2023 duplication

**Concern (peer review §7 citation discipline).** Both `Vollger2023` and `concerted_evolution_nahr_Vollger2023` exist in REFERENCES_v5.bib. The peer reviewer asks: are these duplicate entries for the same paper?

**Why not addressed in v4.** Out of scope for the draft-revision pass. The bib is the citation universe; v4 cites only `concerted_evolution_nahr_Vollger2023` (P12). `Vollger2023` is orphan in the cited-set for v4.

**What would close it.** A bib-hygiene pass that DOI-dedupes REFERENCES_v5.bib. If `Vollger2023` and `concerted_evolution_nahr_Vollger2023` share a DOI, collapse to a single key and update every reference in the synthesis directory.

**Effort.** ~0.5 day.

## D-Bib2. Smith 1976 ectopic-recombination foundation missing

**Concern (peer review §7 citation discipline).** Smith 1976 ("Evolution of repeated DNA sequences by unequal crossover") is the canonical mechanism reference for ectopic recombination between paralogous sequence and is not in REFERENCES_v5.bib.

**What would close it.** Add Smith, G. P. (1976). Evolution of repeated DNA sequences by unequal crossover. *Science* **191**, 528-535. (DOI 10.1126/science.1251186) to REFERENCES_v5.bib and cite in P12 next to Dover 1982.

**Effort.** ~10 minutes.

## D-Bib3. Cech 2004 telomere-clustering review missing

**Concern (peer review §7 citation discipline).** "Cech 2004 (telomere-clustering review) is missing despite the bouquet emphasis."

**What would close it.** Search PubMed for Cech-authored 2004 telomere-clustering reviews and add the matching entry to REFERENCES_v5.bib. Cite in P7 next to `bouquet_KotaSUN1MAJIN2020` / `ZicklerKleckner2015`.

**Effort.** ~30 minutes (lookup + bib insertion).

## Summary

| Concern | Class | Effort | Target round |
|---|---|---|---|
| D-M4  | new analysis (Monte Carlo null) | 1 day | next analysis pass |
| D-M5  | new analysis (Mantel) | 2-3 days | next analysis pass |
| D-M6  | new analysis (matched F_ST) | 2 days | next analysis pass |
| D-M9  | new analysis (character bootstrap) | 3 days | next analysis pass |
| D-M12 | bootstrap CIs (5 statistics) | ~3 days total | next analysis pass |
| D-PeerQ1 | additional control (MAPQ-strict) | 1-2 days | next analysis pass |
| D-PeerQ2 | revision-letter response | 0.5 day | cover letter |
| D-PeerQ2b | new analysis (all-SD Mantel) | 1-2 weeks | separate paper |
| D-PeerQ3 | descriptive table (P/Q asymmetry) | 0.5 day | next analysis pass |
| D-Bib (Vollger dedup) | bib hygiene | 0.5 day | bib pass |
| D-Bib2 (Smith 1976) | bib addition | 10 min | bib pass |
| D-Bib3 (Cech 2004) | bib addition | 30 min | bib pass |

Headline-numbers most in need of CIs and nulls before the next submission round: pedigree 92%, mouse ρ = 0.715, F_ST 0.10-0.15, Mantel ρ = 0.66.

---

## Narrative-match items deferred (added v5, 2026-05-17)

Items below were flagged in `paper_prep/synthesis/NARRATIVE_MATCH_PLAN.md` (40-finding aggregation of 5 brutal reviews) but not fully applied in NATURE_DRAFT_v5.md. Each is logged here with status and reason. Full per-finding before/after for the 31 APPLIED items is in `paper_prep/synthesis/REVISION_LOG_v5.md`.

### F20 — DUX4 cancer / oncofetal-programme angle in P1 — DEFERRED (no cite)

The narrative-match plan asked for a clause "DUX4 reactivation is also an oncofetal programme in multiple cancers" alongside the FSHD framing. REFERENCES_v5.bib does not contain a DUX4-cancer / DUX4-oncofetal citation; the closest match (`dux4_d4z4_fshd_geng2012`) is FSHD-pathophysiology with a side mention of germline-gene activation but not explicit oncofetal. The v5 task instruction forbids adding new bibkeys ("do not add new bibkeys, only reuse existing ones"). Action: add a DUX4-cancer review citation (candidate: Yao et al. 2014 PNAS "DUX4-induced gene expression is the major molecular signature in FSHD skeletal muscle" or a more direct cancer-context review) to REFERENCES_v6 and re-apply F20 in v6. Effort: 0.5 day bib hygiene + 1-line edit.

### F21 — Hi-C rare-contact-regime justification — DEFERRED (word-budget)

NARRATIVE_REVIEW_qa-coverage Finding 3 flagged Q1's "Hi-C is notoriously bad for measuring interchromosomal contacts" first-half question as unanswered: the paper addresses the MAPQ0/telomeric-contact half (F07 applied) but not the intrinsic rare-contact-regime issue (~2-5% of read pairs, enriched for random ligation artefact). The plan's recommended one-sentence fix would justify observed-over-expected normalisation and cite the 14-dataset multi-individual consistency as evidence against correlated artefact. With main text at 3295/3300 words there is no room for a 25-30 word addition. The point is partially implicit in the existing observed-over-expected language and in the F07 limitation acknowledgement. Action: in v6, when other should-fix items have been compressed further, add one sentence to P9 or Methods §Hi-C explicitly justifying the O/E normalisation for the rare-contact regime. Effort: 15 minutes once budget exists.

### F30 — End P14 on directionality, not on experiments list — DEFERRED (structural)

NARRATIVE_REVIEW_narrative-arc Finding 9 flagged that the talk ended on the chicken-or-egg as the intellectual puzzle, while v5 P14 still ends on the experiments list ("Long-read recombination maps in trios, matched germline LAD data, and full CEPH1463 cross-assembler analysis will close the remaining open links..."). The directionality language (now with F17's chicken-and-egg label) sits mid-paragraph in v5. Re-ordering would require duplicating the chicken-and-egg sentence at the paragraph end or repositioning the limitations list, which interferes with the Nature limitations-section convention. DEFERRED to v6 author judgement on whether to restructure P14 to end on the directionality question.

### F32 — "Simulate the full graph without building it" closing clause — DEFERRED (low value)

The Erdős-Rényi P2 paragraph already concludes "so transitive closure recovers virtually every subtelomere in the dataset (Methods)" which conveys the same substance as the talk's "simulate the full graph without having to build it." The verbatim talk phrasing is rhetorically tighter but adds no information. DEFERRED indefinitely; not a substantive narrative-match gap.

### F34 — Per-meiosis per-Mb crossover rate — DEFERRED (input data missing)

NARRATIVE_REVIEW_qa-coverage Finding 8 noted that the talk's "ongoing and frequent recombination exchange" is unquantified in the paper. Converting the 16 crossover-like patches to a per-meiosis per-Mb rate requires N transmissions per pedigree and total PHR length surveyed per parent-child pair, neither of which is reported in v4 or v5. Computing this rate is straightforward (~1 day) but requires running the pedigree pipeline rate-calibration analysis. Action: add a per-meiosis crossover rate to a follow-up pedigree-rates analysis; report as "X crossovers per meiosis per Mb of PHR sequence, 95% CI [a, b]" in v6 main text. Effort: 1 day analysis + 1-line edit.

---

| Concern (v5 narrative-match) | Class | Effort | Target round |
|---|---|---|---|
| F20 — DUX4 cancer angle | bib addition + 1 clause | 0.5 day | v6 bib + edit pass |
| F21 — Hi-C rare-contact justification | 1-sentence add (word-budget gated) | 15 min | v6 if budget |
| F30 — End P14 on directionality | structural re-order | 0.5 day | v6 author judgement |
| F32 — "Simulate full graph" verbatim | low-value polish | 10 min | indefinite defer |
| F34 — Per-meiosis crossover rate | new analysis | 1 day | next analysis pass |

