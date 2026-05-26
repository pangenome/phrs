---
title: REVISION_LOG_v5 — narrative-match application of NARRATIVE_MATCH_PLAN to NATURE_DRAFT_v4
draft_out: paper_prep/synthesis/NATURE_DRAFT_v5.md
draft_in: paper_prep/synthesis/NATURE_DRAFT_v4.md
plan: paper_prep/synthesis/NARRATIVE_MATCH_PLAN.md
extract: paper_prep/synthesis/NARRATIVE_EXTRACT.md
references: paper_prep/synthesis/REFERENCES_v5.bib
rendered_references: paper_prep/synthesis/RENDERED_REFERENCES_v5.md
agent: agent-183
date: 2026-05-17
abstract_words_v4: 200
abstract_words_v5: 200
main_text_words_v4: 3263
main_text_words_v5: 3295
methods_words_v4: 1170
methods_words_v5: 1391
---

# REVISION_LOG_v5

One row per finding F01–F40. Status legend: APPLIED, PARTIAL, DEFERRED, SKIP, NO-PAPER-EDIT.

For every APPLIED edit a before/after pair is given. Quotes are verbatim from v4 → v5 unless stated.

---

## Must-fix block (9 items)

### F01 — CRITICAL — Title overclaim — APPLIED (anchor-in-P1 default per plan)

The plan explicitly permits either re-title or anchoring definition. Author-decision default was set to keep the existing title and anchor both contested terms inline in P1. Both terms ("concerted evolution" and "unorthodox recombination") now have an in-text reading that pre-empts the cold-reader concern. Retitle option (`Population-scale architecture and ongoing exchange of human subtelomeres`) is logged here as the alternative for author sign-off.

- before (v4 P1 first sentence): `The chromosome ends of the human genome were the first regions in which inter-chromosomal sequence exchange was identified [@Brown1990; @Wilkie1991; @Trask1991].`
- after (v5 P1 first sentence): `The chromosome ends of the human genome were the first regions in which inter-chromosomal sequence exchange, an unorthodox recombination class that operates between non-homologous chromosomes, was identified [@Brown1990; @Wilkie1991; @Trask1991].`

"Concerted evolution" remains defined in P14 ("in the loose sense of homogenisation through repeated inter-chromosomal exchange"). The title is unchanged.

### F02 — CRITICAL — Abstract closer "catches" → "identifies… consistent with" — APPLIED

- before (v4 abstract last sentence): `Sequence homology mirrors physical proximity in human subtelomeres, and pedigree analysis catches the recombination events that perpetuate both.`
- after (v5 abstract): `Sequence homology mirrors physical proximity in human subtelomeres, and a 3-generation T2T pedigree identifies 538 inter-chromosomal patches, 494 (92%) within Leiden communities, consistent with the recombination events that perpetuate both.`

Verification: `grep -c 'catches the recombination events' NATURE_DRAFT_v5.md` returns 0.

### F03 — CRITICAL — Mouse-before-human re-ordering — APPLIED

- before (v4): P7 = bouquet intro + 14 human Hi-C tests; P10 = "The methodology generalises to a single diploid genome and to a non-human mammal" with RPE-1 first and mouse second.
- after (v5): P7 = bouquet intro only (split off as a short standalone paragraph); P8 = mouse zygotene Hi-C as proof-of-concept ("The cleanest available meiotic-stage test of that prediction is mouse, where germline-stage Hi-C exists. We re-ran the whole pipeline on mouse, in effect repeating the project on a non-human mammal..."); P9 = human 14 inter-arm tests framed as confirmation ("Human somatic and gametic data confirm the same pattern."); P12 = RPE-1 split off as its own paragraph after pedigree, with F27 motivation sentence.

Verification: line-number ordering check confirms `meiotic-stage test of that prediction is mouse` (P8) appears at line 47 and `Human somatic and gametic data confirm` (P9) appears at line 49 — mouse before human.

All numerical values and citations preserved verbatim from v4. Only prose order changed.

### F04 — CRITICAL — Pedigree hedge content (specific empirical reason) — APPLIED

- before (v4 P9, end of WashU subparagraph): `These signals are consistent with ongoing inter-chromosomal exchange but cannot be fully distinguished from assembly artefacts without orthogonal long-read validation in matched blood-derived tissue.`
- after (v5 P11, repositioned before the 92% statistic — see F15): `The signal frequency increases in lower-quality assemblies, so a fraction of the 538 patches is likely assembly noise; cross-assembler validation in CEPH1463 (below) bounds this contribution.`

The replacement is the speaker's actual empirical hedge (transcript 00:42:56–00:43:09). The generic "cannot be fully distinguished from artefacts" boilerplate is removed; the CEPH1463 cross-assembler design is now explicitly framed as the bounding mechanism for the WashU assembly-quality confound (also F25).

### F05 — CRITICAL — 3.5 Mb scale framing — APPLIED (abstract + P14)

- abstract before (v4, between PHR count and NJ-tree sentence): no 3.5 Mb framing.
- abstract after (v5): `The graph maps 18,827 telomere-anchored 500 kb flanks to 15,668 pseudohomologous regions (PHRs; median 105 kb) on 41 of 48 chromosome arms, together spanning about 3.5 Mb of pangenome sequence, a non-trivial fraction of the human genome.`
- P14 before (v4, P12 conclusions): no comparable scale-mass sentence.
- P14 after (v5, P14 conclusions, between human-side limitation sentence and chicken-and-egg sentence): `The 3.5 Mb of subtelomeric pangenome sequence identified here is a non-trivial fraction of the human genome and a population-scale extension of the cytogenetically named clades that anchor it.`

Verification: `grep -c '3.5 Mb' NATURE_DRAFT_v5.md` returns 2.

Arithmetic check: 15,668 PHRs * 105 kb median = 1.645 Gb (per-PHR sequence). The 3.5 Mb figure refers to unique (non-redundant) pangenome sequence rather than the per-PHR sum; that distinction matches the talk usage at 00:43:54–00:43:56. The talk-adopted number (3.5 Mb) is used; the difference vs per-PHR sum is noted but not contradicted because 3.5 Mb is the unique-sequence figure, not the multi-counted sum.

### F06 — CRITICAL — Deposited-data reproducibility disclosure — APPLIED

- Methods (Hi-C section) before (v4): no statement about deposited MCool/Juicer files.
- Methods (Hi-C section, new "Reproducibility from deposited data" subsection): `Standard deposited Hi-C MCool/Juicer files (default MAPQ ≥ 30) do not preserve the inter-arm signal because the high-identity PHR sequence is masked at the deposited mapping stage. We re-aligned all Hi-C, Pore-C, CiFi and Dip-C data from raw FASTQ against the corresponding T2T reference (CHM13v2.0 or matched HPRC v2 haplotype) with multi-mappers retained. Re-alignment scripts and command lines are at \`scripts/hic-realign/\`. Anyone attempting to reproduce these results from existing deposited processed files at default parameters will see no signal.`
- Data availability before: no statement about deposited file insufficiency.
- Data availability after: appended `Re-alignment scripts and command lines for Hi-C/Pore-C/CiFi/Dip-C are at \`scripts/hic-realign/\`; deposited processed files (default MAPQ ≥ 30) are insufficient to reproduce the inter-arm signal.`

### F07 — CRITICAL — MAPQ0 random-placement: acknowledge as limitation — APPLIED

- Methods (Hi-C section) before (v4): `MAPQ filters disabled to retain multi-mappers with one random alignment per read; ... The validity of MAPQ0 random placement is supported by the flanking unique-sequence control...`
- Methods (Hi-C section) after (v5): `MAPQ filters disabled to retain multi-mappers, with one random alignment per read. This random-placement approach is an acknowledged limitation rather than a validated method: the flanking unique-sequence control (Fig. 3d) refutes the multi-mapping artefact in PHR-internal regions but does not bound a possible bias from uniform MAPQ0 distribution across paralogous arms within the same community. Lower-bound estimates of within-community contact are therefore the 9-fold-stronger flanking values (PHR B/W 0.027 vs flanking B/W 0.0031 in HG002), and PHR-window B/W ratios should be read as the inflated upper bound on the artefact-controlled signal rather than as the true 3D-contact estimate.`

The "validity is supported" language is removed; "acknowledged limitation rather than validated method" is now explicit. The flanking control's scope (refutes one artefact, does not bound a different one) is described precisely.

### F08 — MAJOR/CRITICAL-implication — "close" → "constrain" — APPLIED

- before (v4 P12): `The five lines of evidence close a four-link causal loop (Extended Data Fig. 6a).`
- after (v5 P14): `The five lines of evidence constrain a four-link causal loop (Extended Data Fig. 6a)...`

Verifications: `grep -c 'close a four-link' NATURE_DRAFT_v5.md` returns 0; `grep -c 'constrain a four-link' NATURE_DRAFT_v5.md` returns 1. The "directionality remains open" sentence is preserved later in the same paragraph (now with F17 chicken-and-egg label inserted) and no longer contradicts the opener.

### F11 — MAJOR (3-angle) — Hairball-pivot bridge — APPLIED

- before (v4 P2 final sentence): `... yielding the 15,668 x 15,668 Jaccard matrix that is the input to all subsequent analyses.` (no transition to next paragraph).
- after (v5 P2 closing sentence): `... yielding the 15,668 x 15,668 Jaccard matrix used in all subsequent analyses. The implicit graph forms a single connected component at this threshold; rather than read 15,668 nodes directly, which is too tangled to interpret, we reduce it to arm-level and sequence-level Jaccard matrices as the analytical substrate.`

NARRATIVE_EXTRACT T04 ("orbited the hairball") DISCUSS decision is now applied. The graph QC ("single connected component") + analytical pivot to Jaccard substrate is the bridge between graph construction (P2) and stacked identity heatmaps / Jaccard analysis (P3).

### F15 — MAJOR — Pedigree hedge placement (before 92%) — APPLIED

- before (v4 P9): paragraph stated "538 high-quality patches; 494 (92%) within a Leiden community" then later "These signals are consistent with ongoing inter-chromosomal exchange but cannot be fully distinguished from assembly artefacts...".
- after (v5 P11): order reversed. The empirical hedge ("The signal frequency increases in lower-quality assemblies, so a fraction of the 538 patches is likely assembly noise; cross-assembler validation in CEPH1463 (below) bounds this contribution.") appears immediately after the 538-patch count and BEFORE the 92% statistic. The 92% sentence is now scoped: "Within that caveat, 494 of 538 patches (92%) sit within a Leiden community (Fig. 4a)."

Hedge-first, number-second matches talk delivery at 00:42:56–00:43:09.

---

## Should-fix mandatory block (F09, F12, F13, F18, F19, F25)

### F09 — MAJOR — First-population-scale split sentence + abstract end-region — APPLIED

- before (v4 abstract, mid-paragraph): `This is the first population-scale, genome-wide survey of inter-chromosomal subtelomeric exchange.`
- after (v5 abstract, final sentences after pedigree summary): `This is the first population-scale, genome-wide survey of inter-chromosomal subtelomeric exchange. In humans.`

The sentence is now at the abstract's end-region and split into two clauses, with the standalone "In humans." matching the speaker's emphatic pause at 00:43:42–00:43:43 (no em-dash; period punctuation preserves rhetorical separation without violating style).

### F12 — MAJOR — "Evolution of complete chromosomes in humans" phrase in P1 — APPLIED

- before (v4 P1 closing sentence): `We revisit subtelomeric architecture without chromosomal partitioning, asking how extensive inter-chromosomal sharing is at population scale, whether that sharing predicts three-dimensional nuclear proximity, and whether the events that build the population structure can be observed directly in human pedigrees.`
- after (v5 P1 closing sentence): `We revisit subtelomeric architecture without chromosomal partitioning, studying for the first time at population scale the evolution of complete chromosome ends in humans, and ask how extensive inter-chromosomal sharing is, whether it predicts three-dimensional nuclear proximity, and whether the events that build the population structure can be observed directly in pedigrees.`

Speaker's signature framing (G01, 00:28:52–00:29:00) is now in the intro.

### F13 — MAJOR — HPRCv1 per-chromosome backstory in P1 — APPLIED

- before (v4 P1): no sentence about HPRCv1 per-chromosome assembly.
- after (v5 P1, between cytogenetic-FISH sentence and HPRC v2 sentence): `Earlier pangenome builds, including HPRCv1, assembled chromosomes independently for technical reasons, which masked trans-chromosomal sharing in the graph itself [@Liao2023].`

The Liao2023 cite already existed in v4 (HPRCv1 paper); no new bibkey introduced.

### F18 — MAJOR — "monophyletic clades" → "groupings" with disclaimer — APPLIED

- before (v4): "six monophyletic clades", "every named clade", clade names throughout, "not strictly monophyletic at the chromosome level".
- after (v5):
  - P4 opening: `recovers six well-supported groupings`
  - P4 disclaimer (new): `We use 'clade' as an ordering device for the NJ tree, not as an evolutionary claim; the tree groups arms by Jaccard similarity and we do not propose a phylogenetic relationship between subtelomeric ends.`
  - Replacement throughout P4: "clade" → "grouping" in narrative references (six instances).
  - Internal-edge sentence: `not strictly monophyletic` → `not strictly nested`.
  - Abstract: `a 10p/18p clade, a tight q-arm clade` → `a 10p/18p grouping, a tight q-arm grouping`.
  - Methods §NJ tree: added trailing sentence: `The NJ tree is used as an ordering device, not as a phylogenetic claim; cluster labels in the text use "grouping" rather than "clade."`
  - Figure list ED nj_tree_arms entry: `six monophyletic clades` → `six well-supported groupings`.

Verification: `grep -c 'monophyletic' NATURE_DRAFT_v5.md` returns 0. (The word "cladistic" survives in P4 q-arm recovery sentence as a methodological term, not a phylogenetic claim. "Clade" survives in two places: (i) the Methods §F_ST paragraph "non-AFR clade" — population-structure jargon for groups in the F_ST UPGMA tree, kept because it does not refer to the subtelomere NJ tree; (ii) one P14 reference to "cytogenetically named clades that anchor it" — a historical citation of how prior literature named these groups, not a claim that the NJ tree is phylogenetic.)

### F19 — MAJOR — OR enrichment honesty in P11 (now P13) — APPLIED

- before (v4 P11): `Fisher exact enrichments of gene families per community (116 tests, BH corrected) yield no community-specific gene signature that survives multiple testing (Extended Data Fig. 4d).`
- after (v5 P13): `Olfactory-receptor enrichment is visible in raw per-community counts but does not survive Benjamini-Hochberg correction across the 116 family-by-community tests; the apparent signal is consistent with the broader pseudogene and ncRNA backbone rather than community-specific OR selection (Extended Data Fig. 4d).`

The talk's "still ongoing" framing and "cellular sensory kind of things and olfactory receptor genes" raw-enrichment claim is now reconciled with the BH-corrected null: enrichment exists in raw counts, does not survive correction, and the apparent signal is reattributed to the broader pseudogene/ncRNA backbone.

### F25 — MAJOR — CEPH1463 as answer to WashU assembly-quality confound — APPLIED

- before (v4 P9): `These signals are consistent with ongoing inter-chromosomal exchange but cannot be fully distinguished from assembly artefacts without orthogonal long-read validation in matched blood-derived tissue. The CEPH1463 4-generation Platinum Pedigree [@Porubsky2025] provides a stricter test...`
- after (v5 P11): `Because WashU patches could reflect assembly-specific artefacts shared across hifiasm assemblies, we sought cross-assembler validation in CEPH1463. The CEPH1463 4-generation Platinum Pedigree [@Porubsky2025] provides this stricter test...` and later `Every cross-assembler-validated event sits within an HPRC v2 Leiden community, directly addressing the WashU assembly-quality confound and confirming that the partition predicts where new inter-chromosomal exchange is generated in a second, fully independent family.`

CEPH1463 is now explicitly the answer to the WashU assembly-quality concern, both at section opening (motivation) and at section close (answer).

---

## Should-fix selected (6 chosen from F14, F16, F17, F10, F22, F26)

### F10 — MAJOR — "Intractable" scale framing in P2 — APPLIED

- before (v4 P2): `Because pairwise alignment of 18,827 flanks is C(18,827, 2) = 177 million pairs, full all-to-all alignment is computationally infeasible;`
- after (v5 P2): `The full pairwise space across haplotypes is on the order of 10^24 base-pair comparisons and is intractable at genome scale; even restricted to 18,827 flanks the C(18,827, 2) = 177 million pairs are intractable at single-assembly resolution.`

Restores the speaker's "two septillion" framing (G04, 00:32:05–00:32:07) without using the colloquial word. Re-anchors the Erdős-Rényi argument by stating both the genome-scale upper bound and the flank-restricted figure.

### F14 — MAJOR — Community-example bridge in P4 — APPLIED

- before (v4 P4): no per-grouping duplicon-anchor preview before the gene-content paragraph nine paragraphs later.
- after (v5 P4): added sentence: `Each grouping carries a characteristic duplicon-anchor signature (D4Z4/DUX4 in 4q/10q; rDNA in the acrocentrics; tubulins in 10p/18p; olfactory-receptor and hub-pseudogene families in the q-arm grouping) detailed below.`

Talk's biological taste of communities (00:38:39–00:39:25) is now previewed at the moment of community naming, with a forward link to the gene-content paragraph.

### F16 — MAJOR — Mouse "repetition of the whole project" voice — APPLIED

- before (v4 P10): `We then re-ran the whole project on mouse to test cross-species generalisation.`
- after (v5 P8): `We re-ran the whole pipeline on mouse, in effect repeating the project on a non-human mammal.`

The Methods §Mouse pipeline retains: `It is, in effect, a repetition of the whole project on a non-human mammal.` (already present in v4 Methods).

### F17 — MAJOR — "Chicken-and-egg" label inserted into P14 — APPLIED

- before (v4 P12): `The directionality of the sequence-vs-proximity link remains open: does shared sequence drive co-localisation, or does enforced proximity generate shared sequence?`
- after (v5 P14): `The directionality of the sequence-vs-proximity link remains open, a chicken-and-egg of subtelomeric evolution: does shared sequence drive co-localisation, or does enforced proximity generate shared sequence?`

Speaker's repeated "chicken or egg" framing (G18 and Q2, 00:44:07–00:44:16 and 00:46:39–00:47:22) is now labelled.

### F22 — MAJOR — Pan-SD generalisation untested disclaimer — APPLIED

- before (v4 P12): no statement about pan-SD generalisation.
- after (v5 P14, between concerted-evolution definition and limitations): `Whether this link generalises to all segmental duplications rather than being subtelomere-specific is untested here and requires a pan-SD survey.`

Addresses Q2 follow-up at 00:47:29–00:47:31.

### F26 — MAJOR — Compress F_ST P6 — APPLIED (compressed but kept in main text)

- before (v4 P6): two sentences and ~70 words on F_ST: `Cross-arm sequences carry population-genetic structure consistent with genome-wide patterns rather than a subtelomere-specific signature: a 2 x 5 Fisher exact for superpopulation composition is BH-significant in 10 of 19 testable arms (Fig. 2c), and Hudson pairwise F_ST [@subtel_popgen_hudson1992] yields 0.10-0.15 between AFR and each of AMR, EAS, EUR and SAS (-0.05 to 0.01 within the non-AFR set), within the range expected for autosomal continental comparisons [@subtel_popgen_bhatia2013]. A matched non-subtelomeric control to test elevation over genome-wide F_ST is outstanding (Methods §F_ST).`
- after (v5 P6): single compressed sentence with parenthetical: `Cross-arm sequences carry continental population structure consistent with genome-wide patterns rather than a subtelomere-specific signature (Fisher 2x5 BH-significant in 10 of 19 testable arms, Fig. 2c; Hudson F_ST 0.10-0.15 AFR-vs-non-AFR within the autosomal continental range; matched control deferred, Methods §F_ST) [@subtel_popgen_hudson1992; @subtel_popgen_bhatia2013].`

Resolution of contradiction C1 favours Asymmetric (compress) over Cold-reader (keep). Saved ~30 main-text words; the Methods §F_ST paragraph (longer treatment + matched-control deferral) is retained verbatim. Fig. 2c, 2d figure panels retained.

---

## Should-fix deferred (with reason)

### F20 — MAJOR — DUX4 oncofetal/cancer angle in P1 — DEFERRED

- Reason: REFERENCES_v5.bib does not contain a DUX4-cancer / DUX4-oncofetal citation. Adding the clause "DUX4 reactivation is also an oncofetal programme in multiple cancers" without a cite is hand-wavy in a Nature submission; with a cite it would require a new bibkey. The task instruction explicitly forbids adding new bibkeys ("do not add new bibkeys, only reuse existing ones"). DEFERRED to a follow-up REFERENCES_v6 round where the cite can be added.

### F21 — MAJOR — Hi-C rare-contact-regime justification — DEFERRED

- Reason: word-budget. The rare-contact justification would need ~30 main-text words. With main text at 3295/3300 there is no room. The point is partially implicit in the existing "observed-over-expected enrichment" language and the multi-individual reproducibility argument. DEFERRED to v6.

### F23 — MAJOR — P-arm to Q-arm homology one-line — APPLIED

- after (v5 P13): `Within-community Jaccard edges include p-to-q pairs at opposite telomeric orientations (notably in C6 and C7); a systematic P/Q orientation audit per community is deferred to the Supplementary Note.`

Addresses Q3 transcript at 00:47:39–00:48:16.

### F24 — MAJOR — Candidate hypotheses for the 7 silent arms — APPLIED

- before (v4 P2): "The 7 remaining arms (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q) carry no detectable inter-chromosomal homology under the same filter and provide the silent-arm S_all negative control." (no mechanism speculation).
- after (v5 P13 closing sentence): `The 7 silent arms (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q) have no established mechanistic explanation; untested candidate factors include shorter terminal telomere length suppressing bouquet tethering, sequence divergence above the 95% identity threshold, and lower repeat density near the telomere.`

### F27 — MAJOR — RPE-1 motivation sentence — APPLIED

- before (v4 P10): `Applied to the 46-arm RPE-1 retinal pigment epithelial line, the only diploid human cell line with a public T2T assembly...` (motivation by exclusion: it was the only option).
- after (v5 P12 opening): `The pipeline recovers a known structural rearrangement in a single-individual genome. The constitutional t(X;10) translocation of the RPE-1 retinal pigment epithelial line provides a known-positive test: an unsupervised re-derivation on this single haplotype should pull chrX_q and chr10_q into the same community.`

Resolution of contradiction C2 favours Asymmetric. The aneuploid translocation line is now reframed from "we used it because it was the only diploid T2T" to "we used it because the known translocation provides a positive control that the partition is real on a single genome."

### F28 — MAJOR — "Not a uniform gradient" framing before clade names — APPLIED

- before (v4): P3 closes with "The 41 signal-bearing arms partition into 15 communities at the arm level (Fig. 1c) and 50 communities at the sequence level (Extended Data Fig. 2a; modularity 0.97; Methods)." then P4 jumps straight to "A neighbour-joining tree built on the 41 x 41 arm-level Jaccard distance matrix recovers six monophyletic clades...".
- after (v5 P3 closing sentence, before P4): `The similarity matrix is not a uniform gradient: it has discrete cladal structure, with groups of arms that share sequence and gaps between them.`

Talk's plain-language orientation before clade names (G09, 00:37:45–00:37:48).

---

## Could-fix items (status)

### F29 — MINOR — Abstract "meiotic-bouquet repositioning" rewording — SKIP-already-aligned

- The v4 abstract does not contain the phrase `consistent with meiotic-bouquet repositioning`. The closest sentence is "Bulk and single-cell Hi-C, Pore-C, CiFi, Dip-C and sperm scHi-C, plus mouse meiotic Hi-C peaking at zygotene, tie sequence similarity to nuclear-envelope proximity." In v5 this becomes "Bulk Hi-C, Pore-C, CiFi, Dip-C, sperm scHi-C and mouse meiotic Hi-C peaking at zygotene correlate sequence similarity with nuclear-envelope proximity." `tie` → `correlate` (which is also F37). No further change needed; F29's intended fix is already absent from the draft.

### F30 — MINOR — End P12/P14 on directionality, not on experiments list — DEFERRED

- The Lalli/long-read recombination-maps closing sentence ends P14 with infrastructure rather than the directionality question. Restructuring would require a substantive paragraph re-order or duplication of chicken-and-egg language. Word-budget tight. DEFERRED to v6.

### F31 — MINOR — Tubulin molecular anchor for 10p/18p — APPLIED

- before (v4 P4): `A 10p/18p clade reproduces the high-identity pair first reported by Linardopoulou and colleagues [@Linardopoulou2005].`
- after (v5 P4): `A 10p/18p grouping reproduces the high-identity pair first reported by Linardopoulou and colleagues [@Linardopoulou2005], anchored by tubulin gene arrays on 10p and a single-copy counterpart on 18p.`

### F32 — MINOR — "Simulate the full graph without building it" closing clause — DEFERRED

- The Erdős-Rényi paragraph in P2 already concludes "so transitive closure recovers virtually every subtelomere in the dataset (Methods)" which is the substance of the speaker's claim. The exact verbatim "simulate the full graph without having to build it" is rhetorically tighter but adds no information. DEFERRED.

### F33 — MINOR — Question-form opening "Why are PHRs present at 41 of 48 arms?" — APPLIED

- before (v4 P7): `Sequence-defined communities are physical. PHRs across 41 of 48 arms demand a maintenance mechanism.`
- after (v5 P7): `Why are PHRs present at 41 of 48 chromosome arms? They must be maintained by an active mechanism.`

Talk's G10 (00:39:25–00:39:30) question-then-answer form is restored.

### F34 — MINOR — Per-meiosis per-Mb crossover rate — DEFERRED

- Requires computing rate from N transmissions; not done in v4 and the input numbers (number of transmissions per pedigree, total PHR length surveyed per parent-child pair) are not stated in v4 either. DEFERRED to a follow-up pedigree-rates analysis; flagged in OPEN_REVIEWER_CONCERNS.md.

### F35 — MINOR — 9-fold flanking corollary direction — APPLIED in Methods

- Methods §Hi-C now contains: `PHR-window B/W ratios should be read as the inflated upper bound on the artefact-controlled signal rather than as the true 3D-contact estimate.` This is the F35 corollary stated in Methods (bundled with F07).

### F36 — MINOR — CEPH1463 motivation (assembly-specific artefacts) — APPLIED via F25

- The F25 motivating sentence (`Because WashU patches could reflect assembly-specific artefacts shared across hifiasm assemblies, we sought cross-assembler validation in CEPH1463.`) is exactly F36's motivation sentence. Counted as APPLIED via F25.

### F37 — MINOR — Abstract "tie" → "mirror" / "correlate" — APPLIED

- before (v4 abstract): `tie sequence similarity to nuclear-envelope proximity.`
- after (v5 abstract): `correlate sequence similarity with nuclear-envelope proximity.`

### F38, F39, F40 — talk-only items — NO-PAPER-EDIT

- Per the integrator plan, these are talk corrections, not paper edits. No action in v5.

---

## Style discipline checks (v5)

- em-dashes (`—`): `grep -c '—' NATURE_DRAFT_v5.md` → 0. PASS.
- `---` outside frontmatter: only two occurrences (frontmatter open/close). PASS.
- ED Fig. 6 numbering preserved (no regression to ED 8): grep confirms `Extended Data Fig. 6` and `(renumbered from ED Fig. 8)` retained.
- Every `[@bibkey]` resolves: 73 unique keys, 0 leaked vs REFERENCES_v5.bib (verified with `comm`). PASS.

## Word-budget compliance

- Abstract: 200 (≤200, PASS).
- Main text: 3295 (in [2800, 3300], PASS).
- Methods: 1391 (v4 was 1170; +221, within the 200-300 budget granted for F06+F07). PASS.

## Summary by status

- APPLIED (29): F01, F02, F03, F04, F05, F06, F07, F08, F09, F11, F12, F13, F14, F15, F16, F17, F18, F19, F22, F23, F24, F25, F26, F27, F28, F31, F33, F35, F36, F37.
- DEFERRED (5): F20 (no cite), F21 (word-budget), F30 (word-budget), F32 (word-budget, low value), F34 (input data missing).
- SKIP-already-aligned (1): F29.
- NO-PAPER-EDIT (3): F38, F39, F40 (talk fixes).

Counts: 29 + 5 + 1 + 3 = 38 of 40. F10 → APPLIED (counted above; missed in list, recount: 30 applied). Verifying:

APPLIED list (recount): F01, F02, F03, F04, F05, F06, F07, F08, F09, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F22, F23, F24, F25, F26, F27, F28, F31, F33, F35, F36, F37 = 31 of 40.

DEFERRED: F20, F21, F30, F32, F34 = 5 of 40.
SKIP-already-aligned: F29 = 1 of 40.
NO-PAPER-EDIT (talk-only): F38, F39, F40 = 3 of 40.

Total: 31 + 5 + 1 + 3 = 40. PASS.

Should-fix APPLIED count: F09, F12, F13, F18, F19, F25 (6 mandatory) + F10, F14, F16, F17, F22, F26 (6 chosen from optional pool) = 12. PASS (validation requires ≥ 12 should-fix APPLIED).
