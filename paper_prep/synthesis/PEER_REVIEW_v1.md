---
title: Peer review of NATURE_DRAFT_v3 (Concerted evolution and unorthodox recombination of human subtelomeres)
reviewer: Simulated Nature reviewer (peer-review-v1)
draft_under_review: paper_prep/synthesis/NATURE_DRAFT_v3.md
references: paper_prep/synthesis/REFERENCES_v5.bib
date: 2026-05-17
style: brutally direct, first-person, no hedge
---

# Peer review of "Concerted evolution and unorthodox recombination of human subtelomeres"

I read this manuscript cold against the v3 draft, the figure captions for Fig. 1-4, ED 1-5, ED 8 and nj_tree_arms, and the underlying methods notes in `end-to-end-report/report/04_heterogeneity.md`, `05_hic_validation.md`, `06_dipc_validation.md`, `10_limitations.md` and `14_pedigree_recombination.md`. The submission has real biological content but is currently not in shape to be sent out to the wider review pool as an Article. Detail follows.

## 1. Significance assessment

Is this a Nature Article? Maybe, conditional on a major revision. The substantive claim, that 15,668 PHRs across 41 chromosome arms partition into 15 arm-level (50 sequence-level) communities and that those communities co-organise in 3D nuclear space, is the first population-scale extension of the Linardopoulou 2005 / Mefford-Trask 2002 subtelomeric duplicon framework into the HPRC v2 era. The methodological move (no chromosomal partitioning, IMPG transitive closure over all-vs-all wfmash, justified by the Erdős-Rényi connectivity bound) is real and is the strongest novelty. Tying the communities to a four-link causal loop sequence -> 3D -> recombination -> sequence is ambitious and would, if cleanly demonstrated, justify an Article. As submitted, the loop is not cleanly demonstrated: the human side of the loop is somatic Hi-C in lymphoblastoid lines, the meiotic side is a Spearman correlation in mouse zygotene Hi-C, and the pedigree side rests on a 92% within-community statistic that is never compared against any null. The work is more significant than a Nature Communications or Genome Research paper but stops short of the cleanly closed causal loop a Nature Article requires. The authors must close that gap or downgrade the loop claim before this is competitive at this venue.

## 2. Headline summary (editor's bullet)

Authors use 465 near-complete assemblies (464 HPRC v2 haplotypes plus CHM13) to build an implicit pangenome graph of 18,827 telomere-anchored 500 kb flanks via wfmash all-vs-all alignment and IMPG transitive closure, calling 15,668 PHRs that partition into 15 arm-level and 50 sequence-level Leiden communities recovering PAR1, PAR2, the acrocentric short arms, the Linardopoulou 10p/18p pair, a tight q-arm clade and 4q/10q DUX4. They report Mantel ρ ≈ 0.66 between arm-level sequence similarity and inter-arm Hi-C contact in CHM13 and HG002, and a Spearman ρ = 0.715 between mouse subtelomeric Jaccard and zygotene Hi-C contact. A T2T 3-generation pedigree contributes 538 inter-chromosomal patches, 92% of which sit within Leiden communities, plus 16 crossover-like and 133 gene-conversion-like events. The submission has real biology but is over word limit, has an unformatted reference list, missing ED figures 6-7 without comment, and a Mantel-vs-B/W internal contradiction in the flanking-paradox claim. I would not send this to wider review without major revision.

## 3. Major concerns

**M1. The Reference list is unformatted.** Lines 99-229 of the draft are bare bibkeys (e.g. `acrocentric_Altemose2022`, `Linardopoulou2005`). A Nature submission cannot have an unrendered bibtex citation list as its Reference section. This is not "fix in production"; this is a basic submission-quality failure that signals the submission is unfinished. **What would fix it.** Render the bibkeys against `REFERENCES_v5.bib` into Nature numerical format, with full author lists up to six, journal, volume, page, year. There are 131 inline cites in the draft but 372 keys in v5; the rendered list must contain exactly the 131 cited.

**M2. Word count is over the Nature Article cap.** Front matter declares `main_text_words: 4158`. Nature Articles are capped at ~3000 words (Methods are extra). The draft is 39% over budget. The most expensive paragraphs are the population-genetics paragraph (line 41, ~1200 words) and the 3D-validation paragraph (line 43, ~1100 words). **What would fix it.** Cut to 3000. Move the F_ST UPGMA tree and the multi-resolution Mantel walk-through to Methods or Extended Data. Compress the citation stacks; many sentences carry 6-9 references (line 41 ends with 10 citations, line 53 has another 7-citation stack). One or two carefully chosen citations per claim.

**M3. The flanking-paradox claim contradicts itself depending on which statistic you read.** In §3D the authors state that PHR B/W = 0.027 in HG002 Hi-C falls to flanking B/W = 0.0031 (a 9× strengthening), arguing that the 3D signal is not multi-mapping artefact (Fig. 3d). But the underlying `05_hic_validation.md` shows the Mantel ρ moves in the opposite direction: HG002 PHR Mantel ρ = 0.657 vs flanking Mantel ρ = 0.520, CHM13 PHR ρ = 0.656 vs flanking ρ = 0.522. B/W strengthens at flanking; Mantel ρ weakens at flanking. The paper reports only the favourable statistic. A reviewer who reads the methods notes will catch this. **What would fix it.** Report both statistics for flanking. Explain why B/W and Mantel ρ disagree (almost certainly the binary thresholding in B/W is more sensitive to a few high-contact pairs than rank correlation). Either rewrite the multi-mapping argument or acknowledge that the two metrics tell different stories.

**M4. The "92% within Leiden community" pedigree statistic has no null baseline.** §pedigree reports 494/538 = 92% of WashU HQ inter-chromosomal patches sit within a Leiden community (Fig. 4a), and frames this as evidence that the population partition predicts where new recombination is found. But the null is never given. With 15 arm-level communities of unequal size, the expected within-community fraction under random pairing depends on the arm-size distribution: if 70% of flank coverage sits in three large communities (C3, C6, C11, C7), the random null could easily be 50-70%. 92% over 50% is not nothing, but it is not what 92% over uniform-prior reads like. **What would fix it.** Compute the expected within-community fraction by Monte Carlo permutation of patch sources to random arm pairs (10,000 reps), report the observed 92% vs the bootstrap null distribution mean and 95% CI, and give the actual p-value for the depletion of cross-community patches. The same applies to the 11 CEPH1463 cross-assembler features.

**M5. The mouse ρ = 0.715 is computed on 344 non-independent PHR pairs without distance permutation.** Methods state that Mantel was used for human (10,000 row-and-column permutations) but the mouse section (line 49) reports "Spearman ρ = 0.715, p = 4.4 x 10⁻⁵⁵, n = 344 inter-chromosomal pairs" as a flat correlation. Per-PHR-pair Jaccard vs zygotene Hi-C contact across 344 pairs is structurally non-independent (pairs share arms, share PHRs, and spatial autocorrelation along chromosomes is severe). Reporting an ordinary Spearman p-value on this is statistically inappropriate. **What would fix it.** Run a Mantel test on the arm-level mouse distance matrix vs the arm-level mouse Hi-C matrix; report ρ with 10,000 row-and-column permutations. The ρ may stay near 0.7 but the p-value will be larger by orders of magnitude. The current p = 4.4 x 10⁻⁵⁵ is uninterpretable as written.

**M6. F_ST 0.10-0.15 between AFR and non-AFR is indistinguishable from background and does not support "Inter-chromosomal exchange leaves a population-genetic signature".** Human autosomal F_ST between AFR and non-AFR continental groups sits at ~0.10-0.15 [@subtel_popgen_bhatia2013], which is exactly the reported range. The UPGMA tree of cross-arm F_ST that "recovers an out-of-Africa topology" (Fig. 2d) is therefore not a subtelomere-specific finding; it is what any cross-population genomic comparison gives. **What would fix it.** Compute the matched genome-wide F_ST on a control set of non-subtelomeric autosomal regions of equivalent length per superpopulation pair, and report the subtelomeric F_ST as a difference or ratio against that baseline. If subtelomeric F_ST is not elevated above genome-wide, drop the population-signature claim. If it is, that is the actually interesting effect size.

**M7. The "four-link causal loop" is not closed in human.** §integrated (line 53, ED Fig. 8a) advertises a closed loop: sequence -> 3D proximity at the bouquet -> ectopic recombination -> new shared sequence -> back to sequence. Of the four links, only three are measured in human: sequence sharing (PHRs), 3D proximity (bulk Hi-C in LCLs, Pore-C, CiFi, Dip-C, and 20-cell sperm scHi-C), and the recombination output (pedigree patches). The crucial fourth link, 3D proximity at the meiotic bouquet, is supplied only by mouse zygotene Hi-C. The 20-cell sperm scHi-C is post-meiotic and the bulk LCL Hi-C is somatic interphase. The authors acknowledge this in limitation (i) but the abstract and §integrated write as if the loop is closed. **What would fix it.** Soften the abstract claim from "consistent with meiotic-bouquet repositioning" to "with the meiotic-bouquet model supported in mouse and indirectly inferred in human via sperm scHi-C." Restate the loop in §integrated as four links, three measured in human and one in mouse, and explicitly call out the missing human meiotic Hi-C measurement.

**M8. Sample selection is opaque.** Line 33: "one HPRC haplotype assembly was excluded for quality." Which one, on what criterion, by what threshold? Methods do not specify. Line 95 (limitations) lists a chr18q chimera control (NA18982#1) but that is a different exclusion. A Nature submission must state every exclusion with sample ID and quantitative threshold. **What would fix it.** Add a Methods sub-section "Sample exclusions" listing each excluded haplotype, the QC metric, the threshold, and a one-line justification.

**M9. The 100% bootstrap support claim on the NJ tree is not a real bootstrap.** The nj_tree_arms README is explicit: "True character-level bootstrap was not possible because the input is a derived distance summary rather than an alignment." Instead the authors add Gaussian noise at sigma = 25% of the off-diagonal IQR to the distance matrix and call it a 1000-replicate "perturbation bootstrap". This is a sensitivity analysis to perturbations of the distance summary, not a phylogenetic bootstrap. The 100% MRCA support number sits in the main text (line 37) and the abstract (line 27) without that caveat. **What would fix it.** Either (a) run a true character-level bootstrap by resampling PHRs with replacement and recomputing the Jaccard matrix per replicate (this is what a phylogenetics reviewer expects), or (b) call the procedure what it is (a distance-matrix sensitivity analysis) and report support values labelled "sensitivity support" not "bootstrap support". Option (a) is the right one; the data exists to do it.

**M10. ED Figures 6 and 7 are missing without comment.** The Figure list (line 237-243) jumps from ED 5 to ED 8 with no acknowledgement. A reviewer will read this as accidentally dropped figures. CLAUDE.md notes the slots are intentionally empty during development, but the submission must not present a gap. **What would fix it.** Either fill ED 6 and ED 7 with the figures the storyline needs (the obvious candidates are a per-arm gene-content overlay and a per-community DUX4/D4Z4 copy-number panel that currently sits in ED 4c), or renumber so the ED set runs ED 1-6 with no gap. The current gap will be flagged by every reviewer.

**M11. "Concerted evolution" is loosely defined and the citation pool is one-sided.** The authors cite Arnheim 1980, Ohta 1984 and Charlesworth 1994 for the term (line 53), but classical concerted evolution under molecular drive (Dover 1982, Trends Genet) is the canonical reference and is missing. The text uses the term to mean "ongoing inter-chromosomal recombination that homogenises sequence between non-homologous chromosomes" which is closer to non-allelic homologous recombination (NAHR) than to concerted evolution sensu stricto. **What would fix it.** Either (a) add Dover 1982 and one or two NAHR-vs-concerted-evolution mechanism reviews and explicitly distinguish the two, or (b) drop "concerted evolution" from the title and abstract and use "non-allelic homologous recombination" throughout. Option (b) is honest; option (a) is defensible if the authors really want to invoke the term.

**M12. Reporting standards: no confidence intervals on key correlations.** Mantel ρ = 0.66 is reported with permutation p-values but no 95% CI. Spearman ρ = 0.715 in mouse is reported with a degenerate p-value and no CI. The pedigree 92% has no CI. The F_ST 0.10-0.15 has no CI. For an N = 6 Hi-C study with strong claims, point estimates without intervals are inadequate. **What would fix it.** Report bootstrap 95% CIs for every Mantel ρ, every Spearman ρ over PHR pairs, every B/W ratio, and the pedigree within-community fraction.

## 4. Minor concerns

**m1. Abstract is 214 words; Nature caps abstracts at 200.** Trim 14 words. The sentence "consistent with meiotic-bouquet repositioning" already needs softening per M7; do both at once.

**m2. The sentence "the median PHR (105 kb) is 31% of the length of PAR2 (334 kb)" (line 35) is a memorable framing.** Keep it but add the obvious caveat that PAR2 is fully homogenised by obligate crossover while most PHRs are not. As written it conflates length with biological status.

**m3. Cite-stack overflow at line 41 (10 citations end of paragraph) and line 53 (7 citations end of paragraph).** Cut to two or three each. The current density looks like padding and will be read as such.

**m4. "Wide copy-number diversity" for DUX4 (line 27) gives the range 0-22 in ED 4c but not the median.** Report the median and IQR in the abstract or first DUX4 sentence so the reader has a quantitative anchor.

**m5. The PBMC Dip-C negative control (line 45) is presented as supporting absence of artefact (W/B = 0.983, p = 0.305) but is N = 18 with hg19 coordinate projection noise.** This is statistically underpowered. Either present as a true negative or remove. The current presentation oversells a null result.

**m6. Line 51 "no community-specific gene signature that survives multiple testing" sits in the gene-content paragraph without a clear takeaway.** State the conclusion explicitly: "subtelomeric communities are not defined by shared protein-coding gene content but by shared pseudogene and ncRNA duplicon backbone". As written the negative finding floats.

**m7. The sandwich_same_hap (115 patches) and complex (1 patch) categories in §pedigree are not explained in the main text.** §14_pedigree_recombination.md defines them but the main draft only names them in passing. Either define in the main text (one short sentence each) or move the count to Methods.

**m8. Per-arm "cross-arm sequence rate" of chrX_q 99.7%, chr21_p 94.0%, chr11_p 74.1% (line 39) is reported without explaining how the "rate" was computed.** §04_heterogeneity.md explains it as (cross-arm sequences / total sequences for that arm), but the main text reader has no way to interpret 99.7% without that denominator. Add one parenthetical sentence.

**m9. The Methods software-versions list (line 91) omits the IMPG version.** Line 63 cites IMPG by author but the Methods version block lists only wfmash, impg commit hash, pggb, odgi, samtools, bedtools, bgzip, hicexplorer, R/ape/vegan/cluster. Make the IMPG commit hash visible alongside wfmash and pggb so the row reads consistently.

**m10. Methods "Pedigree odgi-untangle" subsection (line 87) mentions `scripts/pedigree/analyze-pedigree-recombination.py` by path.** Either include the script in a code-availability section or remove the path; a private filesystem path is not a reviewer-actionable reference.

**m11. The RPE-1 t(X;10) "rediscovered from sequence alone" framing (line 49) overstates novelty.** RPE-1 carries the t(X;10) by construction (this is a well-known karyotypic feature of the line). Reframe as: "the known t(X;10) is recapitulated by an unsupervised Leiden partition of the single-individual distance matrix, demonstrating that the pipeline does not require a population to detect a translocation."

**m12. ED Fig. 4d footnote on OR4F pseudogenisation (line 99 of `paper_prep/figures/ed4/caption.md`) cites a population mean of 62.1% but the main text (line 51) does not quote this number.** Either include the 62.1% in the main text or drop the OR4F per-arm range; otherwise the reader cannot connect ED 4d to the main text claim.

**m13. The "33% of mouse PHRs saturate the 1 Mb extraction window" sentence (line 49) is buried in the cross-species generalisation paragraph but is methodologically important.** Promote one short sentence to the Methods (mouse pipeline section) and add a window-sensitivity panel (currently described in §10 limitation 18) to ED 5 or a new mouse-specific ED.

**m14. Limitations list at line 95 has only 6 of the 7 promised in line 53 ("Seven limitations bound the inference").** Count the limitations in §10_limitations.md (18 numbered) vs the Methods list (i)-(vi) vs the main-text discussion (seven). Reconcile.

**m15. The internal-citation density (131 unique inline keys for a 4158-word main text) is high even by Nature standards.** Combined with the cite-stacks (M3 minor m3), expect the production editor to ask for compression. Aim for ~60-80 unique references in the final.

**m16. "Bouquet" is mentioned 5 times in the main text but the actual bouquet anatomy (TERB1-TERB2-MAJIN trimer, SUN1-KASH5 LINC) is described in one citation cluster on line 45.** A reader unfamiliar with meiotic chromosome biology will not parse what "bouquet" means at first mention (line 41 abstract context). Add one short clause defining bouquet at first occurrence in the main text, around line 43.

**m17. The 8/34/7 vs 4/28/9 architectural-split refinement (line 39) reframes a result without sufficient framing.** Make explicit which numbers come from prior FISH literature vs from this work; right now both pairs are presented as Mefford/Linardopoulou-derived.

**m18. "Stacking by chromosome contributor (Fig. 1b) makes the same landscape quantitative" (line 35) is unhelpful prose.** State what Fig. 1b actually shows that Fig. 1a does not: the per-position count of partner chromosomes. One short clarifying clause.

## 5. Statistical rigor check

For every quantitative claim, I asked: is the test appropriate, is the sample size adequate, is multiple-testing handled, is a confidence interval reported?

| Claim (line) | Test | N | Multiple testing | 95% CI | Issue |
|---|---|---|---|---|---|
| Allele vs paralog Wilcoxon p < 10⁻³⁰⁰ (line 41, Fig. 2a) | Wilcoxon paired signed-rank | 5,946 pairs | Per-community + overall; not stated whether BH applied | None | OK if N is right; need to confirm pairs are truly independent (one per individual-community) |
| C7 paralog-closer p = 2.0 x 10⁻⁷ (line 41) | Wilcoxon paired | 156 pairs | None | None | Underpowered relative to other communities (156 vs 5,946); still significant but report CI |
| Spearman ρ < 0 in 39/48 arms (line 41, Fig. 2b top) | Per-arm Spearman | varies | No correction across 48 arms | None | Multiple testing across 48 arms not reported; Bonferroni-adjusted ρ < 0 count would be lower than 39 |
| Piecewise vs linear F-test 39/41 arms (line 41) | F-test on piecewise regression | varies | None | None | 39/41 looks robust; report the median F and a Bonferroni count |
| ITS within 50 kb of breakpoint 16/19 arms (line 41) | Co-localisation count | 19 | None | None | What null? Need a permutation null (random breakpoint placement) |
| Fisher exact superpop 10/19 arms (line 41, Fig. 2c left) | 2x5 Fisher exact | 232 individuals per arm | BH applied | None | OK; report the q-value range |
| F_ST AFR vs non-AFR 0.10-0.15 (line 41, Fig. 2c right, 2d) | Hudson F_ST | 232 individuals | None | None | M6: indistinguishable from genome-wide; need a matched-region control |
| HG002 Hi-C B/W = 0.027, p = 4.0 x 10⁻⁶⁶ (line 43, Fig. 3a) | Mann-Whitney on bootstrap resamples | ~75 arm-haplotypes | None stated | None | Bootstrap inflates effective N; the p-value is wildly significant but the CI on B/W is what matters; report it |
| CHM13 Hi-C B/W = 0.071, p = 6.0 x 10⁻¹⁸ (line 43) | Mann-Whitney | 38 arms (haploid) | None | None | Same |
| HG002 Pore-C B/W = 0.056, p = 3.9 x 10⁻⁸⁵ (line 43) | Mann-Whitney | ~75 arm-haplotypes | None | None | Same |
| Mantel ρ = 0.66 CHM13/HG002 (line 43) | Mantel, 10,000 row-and-col permutations | 41 arms | None | None | OK procedure; need bootstrap CI on ρ |
| Per-individual sequence-pair ρ = 0.83 (line 43) | Spearman | unspecified | None | None | "Lowest-coverage samples where the long-range signal-to-noise is best" sounds cherry-picked; report all samples in a table, not just the favourable one |
| Mantel ρ 0.66 -> 0.80 after exclusion (line 43, ED 5a/b) | Mantel | varies | None | None | A move from 0.66 to 0.80 is real; report the 95% CI and explain why exclusion strengthens rather than weakens |
| O/E enrichment 5.9x - 18.4x (line 43, ED 5c) | Ratio | varies | None | None | Report CI per dataset |
| 15x11 reproducibility heatmap, BH q < 0.05 majority (line 43, ED 5d) | Permutation null with BH | varies | BH applied | None | OK; state the exact number of cells passing q < 0.05 and q < 0.001 |
| C4 (chr7_q-chr12_q) significant in 4/5 diploid samples (line 43) | Per-sample test | 5 samples | None | None | 4/5 is not a Nature-grade result; the alternative "gene-content-driven" is barely falsified at N = 5 |
| PHR B/W 0.027 -> flanking B/W 0.0031 (line 45, Fig. 3d) | Paired comparison | per-sample | None | None | M3: contradicts flanking Mantel ρ; report both metrics |
| Dip-C radial 0.504 vs 0.556, p = 1.6 x 10⁻³⁵ (line 45) | t-test or Mann-Whitney | 2,927 vs 7,267 particles | None | None | Test type unspecified; particles within a cell are not independent (intra-cell correlation must be modelled) |
| GM12878 16/16 C-cells W/B < 1, 0/16 S_all (line 45) | Sign test | 16 cells | None | None | Sign test on 16/16 is p = 1.5 x 10⁻⁵ exact; report the actual p |
| Sperm 20/20 C-cells W/B < 1, 1/20 S_all (line 45) | Sign test | 20 cells | None | None | Same |
| PBMC W/B = 0.983, p = 0.305 (line 45) | Unspecified | 18 cells | None | None | M5: underpowered null |
| Mouse zygotene Spearman ρ = 0.715, p = 4.4 x 10⁻⁵⁵, n = 344 pairs (line 49) | Ordinary Spearman on pairs | 344 PHR pairs | None | None | M5 (major): pairs non-independent, Mantel required |
| ρ 0.574 - 0.715 across four meiotic stages (line 49) | Spearman per stage | varies | None | None | Same |
| Mantel ρ = 0.66 vs 0.485 vs 0.715 (line 53) | Mantel | varies | None | None | OK procedure; report CI |
| Pedigree 494/538 = 92% within community (line 47) | Count | 538 patches | None | None | M4: no null baseline; needs Monte Carlo permutation |
| CEPH1463 11 features (line 47) | Cross-assembler intersection | 11 | None | None | What was the candidate denominator that yielded 11? Without it, 11 is uninterpretable |
| C1 D4Z4 Mann-Whitney p = 5.3 x 10⁻⁶ (line 49 of `paper_prep/figures/ed8/caption.md`) | Mann-Whitney | unspecified | None | None | Test groups undefined in main text |
| TAR1 PAR vs autosomal: chrXp 0.3%, chrYp 1.1%, autosomal >99% (line 41) | Descriptive | per-arm | None | None | No test for the difference; consider Fisher exact or chi-squared |
| Kruskal-Wallis H = 100.89, p = 3.2 x 10⁻¹⁵ telomere length per community (ED 3c) | KW | 13,668 | None | None | OK; ED only, not in main text |
| RPE-1 Mantel ρ = 0.548 self-discovered vs 0.457 transferred (line 49) | Mantel | 46 arms | None | None | Report CI on both |

**Summary statistical issues.**

- Five tests (M3, M4, M5, M6, m5) are statistically unsound as written.
- Confidence intervals are missing from every correlation in the paper.
- Multiple-testing correction is applied to some tests (Fisher superpop) but not to others (per-arm Spearman ρ, per-arm piecewise vs linear F-test, multi-sample Mantel).
- The pedigree 92%, the mouse ρ = 0.715, and the F_ST 0.10-0.15 are the three numbers in the headline that most need to be re-tested before any further round.

## 6. Figure and Extended Data adequacy

I read the captions for Fig. 1-4, ED 1-5, ED 8, and nj_tree_arms (the nj tree has no `caption.md`, only `README.md`).

**Fig. 1.** (a) Stacked identity heatmap, (b) chromosome-contributor heatmap, (c) 41x41 Jaccard heatmap with UPGMA dendrogram, (d) per-arm architecture bars. This is the right content for a landscape figure. **Issues.** Panel (c) bundles two analyses (heatmap + dendrogram) and panel letters in the main text (line 35-39) reference 1a, 1b, 1c, 1d consistently, which is good. Panel (d) "4/28/9" split is referenced (line 39) but the bar layout needs to make the three categories visually obvious; currently the caption says "red/blue/teal" which is colour-coded but a colourblind reviewer cannot distinguish red and teal without checking the PDF. Use a shape or pattern in addition to colour.

**Fig. 2.** (a) Wilcoxon paired allele-vs-paralog bars, (b) per-arm Spearman gradient + piecewise breakpoints, (c) cross-arm superpop Fisher + Hudson F_ST, (d) UPGMA F_ST out-of-Africa tree. **Issues.** Panel (c) packs two heterogeneous panels (significance heatmap on left, F_ST matrix on right) into one panel letter. Split into 2c and 2d, push the UPGMA tree to 2e or to ED 2. Currently the main text (line 41) references "2c, left" and "2c, right" which is bad form for a Nature panel layout.

**Fig. 3.** (a) HG002 Pore-C 50 kb contact matrix, (b) 14-test forest plot, (c) S_all negative control, (d) flanking paradox + Dip-C radial inset. The forest plot is the strongest single panel because it survives the data heterogeneity. **Issues.** Panel (d) carries both the PHR-vs-flanking comparison and the Dip-C radial inset, which are two distinct claims. Either move the radial inset to ED 5 or split into 3d and 3e. The flanking-paradox claim (M3) needs the Mantel ρ added; either as a second tile in 3d or in the caption.

**Fig. 4.** (a) WashU pedigree patches, (b) CEPH1463 11 features, (c) RPE-1 self + CiFi contact, (d) mouse zygotene Hi-C vs Jaccard. All four panel letters are used in the main text (line 47-49). **Issues.** Panel (a) needs the within-community null baseline (M4) shown as a dashed line or a side annotation. Panel (b) needs the candidate denominator for "11 features pass" (m11 above). Panel (d) needs Mantel rather than Spearman if M5 is taken seriously.

**ED 1.** Pipeline schematic, per-arm flank counts, PHR length distribution, chr18q chimera control. Standard, adequate.

**ED 2.** 50-community UMAP, within-community Jaccard bimodality, cross-arm chord, arm-vs-seq confusion. Adequate; the chord plot is referenced in caption but not in the main text (line 39 cites ED 2c but the caption (c) is the chord plot, so this is consistent).

**ED 3.** TAR1, ITS length and motif, telomere length per community, per-arm TAR1 position. ED 3a is cited in main text (line 41). ED 3b, 3c, 3d are not cited in main text. **Fix.** Either cite the remaining panels somewhere in the main text or move them to a Supplementary section.

**ED 4.** GO PHR-only, copy-weighted GO, top-15 high-copy gene families, OR4F pseudogenisation gradient. Main text cites ED 4c and 4d (line 51). ED 4a and 4b are uncited. The GO PHR-only re-run (4a) is methodologically important (corrects a 10x window-size bias in the original GSEA) and deserves one sentence in the main text.

**ED 5.** Multi-resolution W/B, Mantel exclusion, O/E within vs between, 15x11 reproducibility. All four panels cited in main text (5a, 5c, 5d at line 43; 5a and 5b at the same line). Adequate.

**ED 6 and ED 7.** Missing without comment in the Figure list. M10 above.

**ED 8.** Causal loop, D4Z4-CTCF-lamin model, Lalli null, compartment diagnostic. All four panels (8a, 8b, 8c, 8d) referenced in main text (line 45, 53). Adequate.

**nj_tree_arms.** Has no `caption.md`. The README is well-written but not a publication caption. **Fix.** Write a 150-200 word caption for the NJ tree, including the method (ape::nj on the 41x41 Jaccard matrix), the rooting (MRCA of acrocentric short-arm clade), the support method (M9: rename "perturbation bootstrap" to "distance-matrix sensitivity analysis"), and the six monophyletic clades with their one-to-one Leiden mapping.

## 7. Citation discipline

I checked the seven foundational works the editor flagged:

- **Linardopoulou 2005** (`Linardopoulou2005` in bib): cited at line 31 (cytogenetic FISH context), line 37 (10p/18p clade), line 39 (architectural split). Correct placements.
- **Mefford-Trask 2002** (`MeffordTrask2002`): cited at line 31 (BAC walking), line 41 (Flint-Mefford two-domain model, OR4F prediction), line 51 (OR4F prediction). Correct placements.
- **Altemose 2022** (`acrocentric_Altemose2022`): cited at line 31 (acrocentric rDNA recombination), line 37 (acrocentric clade), line 41 (C7 reversal). Correct placements.
- **Bellott 2024** (`sexchrompars_bellott2024`): cited at line 31 (PAR1/PAR2 obligate crossover), line 35 ("largest known fully homogenised inter-chromosomal region in the human genome"). Correct placements.
- **Garrison 2024 PGGB** (`Garrison2024pggb`): cited at line 33 (PGGB graph), line 33 (PGGB methods); cited in Methods line 67. Correct placements.
- **HPRC v2** (`hprc_hprcv2_2025`): cited at line 31 (HPRC v2 release), line 31 again ("subtelomere companion to the HPRC v2 reference assembly publication"), Methods line 57. Correct placements.
- **Vollger 2023** (`Vollger2023` and `concerted_evolution_nahr_Vollger2023`): cited at line 35 (segmental duplication context), line 53 (pedigree concerted evolution). Both keys exist in v5; the simultaneous citation of both (line 53) likely refers to the same paper indexed twice in v5. **Fix.** Verify whether `Vollger2023` and `concerted_evolution_nahr_Vollger2023` are duplicate entries; if so, collapse to one key.

**Missing or under-cited:**

- **Dover 1982 ("Molecular drive", Trends Genet)** is the canonical concerted-evolution citation and is missing. M11 above.
- **Smith 1976 (the "OR / unequal crossing-over" foundational paper)** is the canonical mechanism reference for ectopic recombination between paralogous sequence and is missing.
- **Mefford-Trask 2002** is cited correctly for the two-domain model, but the original Flint et al. 1997 paper is cited (`Flint1997`); good. However, the "Stong 2014" paper (`Stong2014`) on 22q11 subtelomeric organisation is cited at line 35 and line 39 but the subtelomere-specific Stong/Eichler 2014 work is the primary reference for inter-chromosomal subtelomeric organisation in the post-FISH era and is under-emphasised relative to its weight.
- **Cech 2004 (telomere-clustering review)** is missing despite the bouquet emphasis.
- **Garrison 2018** (`Garrison2018`) for variation graphs is cited at line 33, correct.
- **Nurk 2022 / Logsdon 2021 / logsdon2025hgsvc** are cited at line 31 and line 57 for T2T and HPRC context, correct.

**Citation density issues** (already in M3 and m3): line 41 paragraph ends with 10 stacked citations for population genetics; line 53 ends with 7; line 31 carries 14 citations across two sentences. Compression required.

**Orphan keys.** 372 keys in REFERENCES_v5.bib vs 131 inline cites. 241 keys are uncited orphans. This is not a peer-review concern but a hygiene concern: prune REFERENCES_v5.bib to the cited 131 (or, if some are intentionally retained for revision, document the carry-over).

## 8. Recommendation

**MAJOR-REV** (reject in current form, encourage resubmission of a major revision with new analyses).

**Justification.** The biology is novel and the methodological contribution (no-chromosomal-partitioning IMPG transitive closure over wfmash all-vs-all, justified by Erdős-Rényi connectivity, yielding 15,668 PHRs in a single integrated analysis across 465 near-complete assemblies) is real and is the kind of contribution Nature is interested in. The arm-level Leiden partition and its one-to-one mapping onto cytogenetically named clades (PAR1, PAR2, acrocentrics, 10p/18p Linardopoulou pair, q-arm clade, D4Z4) is genuinely striking and is the strongest single result. However, the submission is currently not at the level required by Nature: (i) the Reference list is unformatted (M1); (ii) word count is 39% over Article cap (M2); (iii) the flanking-paradox claim contradicts itself between B/W and Mantel statistics (M3); (iv) the pedigree 92% has no null baseline (M4); (v) the mouse ρ = 0.715 is statistically inappropriate as an ordinary Spearman on 344 non-independent pairs (M5); (vi) F_ST 0.10-0.15 is indistinguishable from genome-wide and does not support the population-signature claim (M6); (vii) the four-link causal loop is not closed in human (M7); and (viii) the 100% NJ bootstrap is not a real bootstrap (M9). Major concerns M1, M2, and M3 are submission-quality failures that I would not personally re-review; the others (M4-M12) are analysis-quality concerns that require new computation. A revised submission that closes M1-M12 and reduces the 18 minor concerns (m1-m18) is, in my judgement, competitive at this venue. The current submission is not.
