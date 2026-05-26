---
title: Narrative Review -- Q&A Coverage Audit
angle: Q&A coverage (talk vs draft)
draft_reviewed: paper_prep/synthesis/NATURE_DRAFT_v4.md
transcript: /home/guarracino/Dropbox/grants/r21-cancer-pangenomics/notes/Session7-PopulationGenomics.en.srt
narrative_extract: paper_prep/synthesis/NARRATIVE_EXTRACT.md
peer_review_ref: paper_prep/synthesis/PEER_REVIEW_v1.md (do not duplicate)
revision_log: paper_prep/synthesis/REVISION_LOG_v4.md
reviewer: narrative-match-qa (agent-167)
date: 2026-05-17
---

# Narrative Review: Q&A Coverage Audit

## 1. My angle

I am auditing whether the paper answers the three Q&A items that an audience member forced into the open at the BoG conference talk: Q1 (Hi-C signal validity at 00:45:12), Q2 (chicken-or-egg causality at 00:46:16), and Q3 (P-arm-to-Q-arm homology at 00:47:39). I am also looking for hard questions the paper invites that the audience did not get to ask. I am NOT re-doing PEER_REVIEW_v1.md's statistical or structural concerns; where I overlap with M3 (flanking paradox), M4 (pedigree null), or M7 (causal loop not closed in human), I note the overlap in one sentence and add only what the peer reviewer missed from the Q&A angle.

---

## 2. Findings

---

### Finding 1 -- CRITICAL

**Q1 answer gap: "We couldn't see any of the signal when we took the deposited data sets" is absent from the paper.**

- **Draft cite:** Methods §Hi-C, Pore-C and CiFi pipeline (draft line ~84, single paragraph). The relevant sentence: "MAPQ filters disabled to retain multi-mappers with one random alignment per read."
- **Transcript cite:** 00:46:06-00:46:11 (Q1 response, subtitles 348-349).
- **Verbatim transcript:** "We couldn't see any of the signal when we took the deposited data sets. We had to realign them in order to do this."
- **Verbatim draft:** "MAPQ filters disabled to retain multi-mappers with one random alignment per read; within-vs-between (W/B) ratio computed by bootstrap on 10,000 permutations."
- **Divergence:** Eric's spoken admission is that the signal is completely invisible in standard deposited Hi-C files. Anyone who downloads existing MCool or Juicer-format files from GEO and runs standard analysis will see B/W ratios near 1 across all PHR communities. The paper omits this entirely. "MAPQ filters disabled" reads as a methodological choice, not as a mandatory re-alignment that invalidates all prior processed-data releases. This is a reproducibility crisis embedded in a conference Q&A that the paper does not disclose.
- **Concrete edit:** Add to Methods §Hi-C: "All source BAM/CRAM files were re-aligned from scratch; the deposited MCool and Juicer files from these experiments systematically exclude MAPQ0 reads and recover no within-community contact enrichment (B/W approaches 1 when MAPQ0 reads are removed). Re-aligned BAM files are deposited at [repository]." Add the same point to the Data and Code availability section. Without this, the paper is irreproducible from deposited data.

---

### Finding 2 -- CRITICAL

**Q1 answer gap: Paper promotes random MAPQ0 assignment as a validated method; Eric said "we probably should try to correct for" it.**

- **Draft cite:** Methods §Hi-C (line ~84): "The validity of MAPQ0 random placement is supported by the flanking unique-sequence control (PHR B/W 0.027 vs flanking B/W 0.0031 in HG002), which strengthens at non-duplicated sequences and so cannot be a multi-mapping artefact."
- **Transcript cite:** 00:45:36-00:45:44 (Q1 response, subtitles 337-339).
- **Verbatim transcript:** "So we've done something that we probably should try to correct for, but we randomly assign the reads to a location. If they map exactly the same between different locations, we randomly assign them to one."
- **Verbatim draft:** "The validity of MAPQ0 random placement is supported by the flanking unique-sequence control."
- **Divergence:** Eric framed random assignment as an acknowledged limitation requiring future correction. The paper promotes it as a validated approach supported by a control. The flip is not defended anywhere in the draft. The flanking control demonstrates that unique-sequence proximity is real; it does not demonstrate that random MAPQ0 assignment produces unbiased within-community contact estimates at the PHR loci themselves. When a read maps equally to chr4q and chr10q (both in community C1), random assignment distributes 50% of that read's contact to each arm -- which inflates within-community B/W for C1 regardless of true proximity. The flanking control refutes the "purely self-mapping within one chromosome" artefact, but not the "uniformly distribute MAPQ0 reads across community members" artefact. Eric knew this and said so; the paper does not.
- **Concrete edit:** Change "The validity of MAPQ0 random placement is supported by..." to: "MAPQ0 reads were assigned to one randomly chosen equally-scoring locus, an approach that distributes multi-mapper contacts uniformly across paralogous loci rather than discarding them. This approach has a known limitation: reads that map equally to two arms within the same community will inflate within-community contact counts regardless of true proximity. The flanking unique-sequence control (B/W = 0.0031 at 100 kb centromere-ward of PHR boundaries, where MAPQ0 cannot apply) provides evidence that the co-localisation signal is not purely an artefact of this distribution, but does not bound the inflation factor at PHR loci themselves."

---

### Finding 3 -- MAJOR

**Q1 answer gap: "Notoriously bad for interchromosomal contacts" gets a one-sided answer.**

- **Draft cite:** P7 (line ~44): "We assembled 14 inter-arm tests of three-dimensional proximity across six HPRC v2 individuals and CHM13: bulk Hi-C in CHM13 and five HPRC samples, Pore-C and CiFi in HG002..." and P8: "The 3D signal is not a multi-mapping artefact."
- **Transcript cite:** 00:45:12-00:45:32 (subtitles 329-335, the questioner's framing).
- **Verbatim transcript (questioner):** "Hi-C is notoriously bad for measuring interchromosomal contacts. And also really bad at finding telomeric contacts. So I'm wondering how you are able to get that signal."
- **Verbatim draft:** "The 3D signal is not a multi-mapping artefact. B/W ratios computed on the unique-sequence 100 kb regions centromere-ward of the PHR boundaries are stronger than those on the PHRs themselves."
- **Divergence:** The questioner raised two distinct concerns: (a) interchromosomal contacts in Hi-C are rare (~1-5% of all read pairs) and enriched for artifactual random ligations; (b) telomeric contacts are specifically bad because repetitive sequence causes MAPQ0 mapping. The paper answers (b) with the MAPQ0 re-alignment. It does not answer (a) at all. Interchromosomal contacts are discarded as noise in most Hi-C pipelines precisely because their frequency is so low relative to background ligation noise. The paper reports B/W ratios that are statistically significant (HG002 B/W = 0.027, p = 4.0e-66) but never explains why its observed-over-expected normalization is adequate for the specific noise regime of rare inter-chromosomal contacts. A reviewer familiar with Hi-C noise modeling will ask this question; it is not answered.
- **Concrete edit:** Add one sentence to P7 or to Methods §Hi-C: "Interchromosomal contacts are intrinsically rare in bulk Hi-C (~2-5% of read pairs) and are enriched for artifactual random ligation; we address this by (i) retaining MAPQ0 reads via random placement (overcoming the telomeric-repeat MAPQ0 filtering problem), (ii) applying observed-over-expected inter-chromosomal normalization, and (iii) demonstrating within-community enrichment across 14 independent datasets and six individuals, making correlated noise across all datasets unlikely."

---

### Finding 4 -- MAJOR

**Q2 follow-up: SD generalizability completely unaddressed in paper.**

- **Draft cite:** Nowhere in NATURE_DRAFT_v4.md is there a sentence testing or even noting whether the sequence-similarity-to-3D-proximity correlation holds for all segmental duplications rather than just subtelomeric PHRs.
- **Transcript cite:** 00:47:29-00:47:31 (subtitle 369, second audience member comment after Q2).
- **Verbatim transcript:** "If you look at this across all segmental duplications, not just the sub-telomeres and in response to Flores' point, be a little more careful about the control of the Hi-C data."
- **Verbatim draft:** (no corresponding passage; the paper discusses only PHRs and subtelomeric communities throughout)
- **Divergence:** The audience member made two points in one sentence: (a) this result may not be subtelomere-specific -- it might hold for all SDs; (b) another audience member ("Flores") raised a Hi-C control concern during the session. The paper does not address (a) anywhere. If the Mantel correlation between sequence similarity and 3D proximity holds for all segmental duplications genome-wide (a plausible null, since SDs by definition share sequence and many known SDs are co-localised), then the subtelomeric PHR result is not novel -- it is one instantiation of a general SD principle. NARRATIVE_EXTRACT §6 Q2 noted this: "The result may not be subtelomere-specific if it holds for all segmental duplications." PEER_REVIEW_v1 did not flag this specific gap. The paper never tests it.
- **Concrete edit:** Add one sentence to Discussion or Limitations: "Whether the sequence-similarity-to-nuclear-proximity correlation is specific to subtelomeric PHRs or is a general property of segmental duplications genome-wide has not been tested; a pangenome-scale SD analysis is required to determine whether PHRs are a subset of a broader co-localisation principle or are structurally distinct."

---

### Finding 5 -- MAJOR

**Q2 causality: chicken-or-egg posed in paper but three mechanistic scenarios not distinguished.**

- **Draft cite:** P12 (final paragraph): "The directionality of the sequence-vs-proximity link remains open: does shared sequence drive co-localisation, or does enforced proximity generate shared sequence? The bouquet provides structural opportunity for both, and resolving the directionality will require tracking proximity and homology across generations."
- **Transcript cite:** 00:44:07-00:44:16 (G18, subtitles 310-313) and Q2 response 00:46:39-00:46:44 (subtitle 357).
- **Verbatim transcript (conclusion, 00:44:07):** "It's a chicken or egg question. Is it that the sequence homology is actually driving the physical proximity, or is it the homology just a result of the fact that the proximity is there?"
- **Verbatim draft:** "The directionality of the sequence-vs-proximity link remains open: does shared sequence drive co-localisation, or does enforced proximity generate shared sequence?"
- **Divergence:** The paper's phrasing is better than the talk's -- it is more precise and symmetrically stated (REVISION_LOG_v4 C.U04 correctly applied this). However, the paper then frames the four-link causal loop in ED Fig. 6a as though it "closes" the argument, while the directionality disclaimer sits in a different paragraph. These two elements contradict each other. A four-link loop that cycles as sequence -> 3D -> recombination -> sequence is consistent with all three scenarios: (a) homology drove proximity initially and now reinforces it; (b) proximity created homology initially and now reinforces it; (c) both operated simultaneously from the start. The loop cannot distinguish between these. The paper should explicitly state that closing a feedback loop does not resolve the founding directionality, and that the loop is equally consistent with all three causal origins. It currently implies closing the loop answers the question when it does not.
- **Concrete edit:** In P12, after the directionality sentence, add: "The four-link loop in Extended Data Fig. 6a is consistent with all three mechanistic origins -- homology drove proximity, proximity drove homology, or both operated simultaneously -- and resolves to the same stable state regardless of which link was the founding event. The directionality question requires either an ancestral reconstruction from species pairs that diverged before the bouquet stage evolved, or a perturbation experiment that disrupts proximity without altering sequence."

---

### Finding 6 -- MAJOR

**Q3 answer: "Yeah, I do actually" (P-to-Q arm homology confirmed in Q&A) vs. complete silence in paper.**

- **Draft cite:** No passage in NATURE_DRAFT_v4.md addresses P-arm-to-Q-arm homology, the arm-orientation asymmetry, or the fraction of within-community edges that are P-P, Q-Q, or P-Q.
- **Transcript cite:** 00:47:39-00:48:16 (subtitles 375-392, Q3 and Eric's response).
- **Verbatim transcript (Q):** "You showed some homology between the P and the P and the Q and the Q. Do you have any examples of homology between a P arm and a Q arm?"
- **Verbatim transcript (Eric):** "Yeah, I do actually... It is true, and we probably should make a test to see if the Qs tend to be closer to Qs and P's tend to be closer to P's. I think it's bound to something actually because you tend to see the stronger relationships between the same arm, the shorter or longer one. To me, it's a bit of a mystery how that would be organized."
- **Divergence:** Eric publicly confirmed P-to-Q homology exists and acknowledged the need to quantify it. The paper contains no quantification. REVISION_LOG_v4 C.T10 explicitly deferred this: "DEFER (to revision letter) -- not addressed in main text. Logged in OPEN_REVIEWER_CONCERNS with the suggested test (fraction of within-community edges that are P-P, Q-Q, vs P-Q) and a supplemental-table action item." A Nature reviewer who watched this talk or reads NARRATIVE_EXTRACT will ask why the publicly acknowledged analysis was not done. Worse, community C7 (acrocentric p-arms) is described as containing five acrocentric p-arms -- but Eric in Q&A pointed to "a cluster you have P and Q inside of it" during Q3, which means even C7 may contain P-Q edges. This is not addressed. Deferring to a revision letter is not a substitute for an answer in the paper.
- **Concrete edit:** Add one sentence to Results (P4 or P5) or to Limitations: "Within-community Jaccard edges include arm-orientation-mismatched P-to-Q pairs; the fraction of P-P, Q-Q, and P-Q edges within each of the 15 communities has not been computed. Preliminary inspection during conference Q&A identified P-to-Q examples in at least one community; a systematic arm-orientation audit is deferred to a supplemental analysis."

---

### Finding 7 -- MAJOR

**Paper invites question: abstract "catches the recombination events" contradicts body's hedge and talk's spoken caution.**

- **Draft cite (abstract):** Abstract final sentence: "Sequence homology mirrors physical proximity in human subtelomeres, and pedigree analysis catches the recombination events that perpetuate both."
- **Draft cite (body):** P9: "These signals are consistent with ongoing inter-chromosomal exchange but cannot be fully distinguished from assembly artefacts without orthogonal long-read validation in matched blood-derived tissue."
- **Transcript cite:** 00:42:56-00:43:08 (subtitles 288-290, G15).
- **Verbatim transcript:** "But I wouldn't say it's conclusive, because unfortunately, when you look at much worse quality assemblies, you see a lot more of this. So it is compatible, but maybe not definitive proof of an actual recombination event."
- **Verbatim abstract:** "pedigree analysis catches the recombination events that perpetuate both."
- **Divergence:** "Catches" asserts detection of real events. "Cannot be fully distinguished from assembly artefacts" and "I wouldn't say it's conclusive" both deny that same detection. The abstract is inconsistent with the body and with the speaker's own live caution. REVISION_LOG_v4 C.U08 correctly added the hedge to P9 but left "catches" in the abstract. PEER_REVIEW_v1 did not flag this internal contradiction. A journal editor will read the abstract first and will flag the discrepancy when they reach P9.
- **Concrete edit:** Change abstract final sentence to: "Sequence homology mirrors physical proximity in human subtelomeres, and pedigree analysis identifies inter-chromosomal exchange signals consistent with the recombination events that perpetuate both." Replace "catches" with "identifies...consistent with" to match the body's hedged language.

---

### Finding 8 -- MAJOR

**Paper invites question: why exactly 7 arms lack PHRs? Zero mechanistic account provided.**

- **Draft cite:** P2 (line ~34): "The 7 remaining arms (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q) carry no detectable inter-chromosomal homology under the same filter and provide the silent-arm S_all negative control."
- **Transcript cite:** 00:36:19-00:36:24 (subtitles 161-162).
- **Verbatim transcript:** "One, two, three, four, five, yeah, seven exceptions. We note they don't seem to be mixing like the other ones."
- **Verbatim draft:** "The 7 remaining arms... carry no detectable inter-chromosomal homology under the same filter and provide the silent-arm S_all negative control."
- **Divergence:** Eric counted the seven exceptions live with zero mechanistic explanation. The paper does not provide one either. The seven silent arms are used only as a negative control; the paper never asks why these specific arms lack PHRs. Are they shorter? Do they have lower telomere-repeat content? Are they telomere-deficient in a subset of individuals? Do they lack the bouquet-anchoring proteins at their tips? Is their subtelomeric sequence diverged above the 95% identity threshold used here? Any of these would be a falsifiable hypothesis. The paper's silence on mechanism is a glaring omission that a reviewer working on meiotic chromosome biology will not miss. This is also a question that was not asked in the Q&A but is immediately invited by the data.
- **Concrete edit:** Add one sentence to Discussion: "The 7 arms that carry no detectable inter-chromosomal homology (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q) lack a mechanistic explanation; candidate factors include shorter average telomere length suppressing bouquet tethering, sequence divergence above the 95% identity threshold, or lower repeat density that reduces opportunity for ectopic exchange, none of which have been tested in this dataset."

---

### Finding 9 -- MINOR

**Paper invites question: de novo inter-chromosomal exchange rate per meiosis is not computed despite pedigree data existing.**

- **Draft cite:** P9 (line ~48): "The filter yields 538 high-quality inter-chromosomal patches; 494 of 538 (92%) sit within a Leiden community (Fig. 4a)... 16 crossover-like single reciprocal exchanges... 133 gene-conversion-like sandwiches."
- **Transcript cite:** 00:41:54-00:41:55 (G13, subtitle 267): "We did attempt to see if we could catch this in the act."
- **Verbatim transcript (conclusion segment, 00:44:22-00:44:25, subtitle 315):** "a completely assembled human pedigree suggests that there may be ongoing and frequent recombination exchange in these regions."
- **Verbatim draft:** "16 crossover-like single reciprocal exchanges (the largest spans 27.97 kb on PAN028 maternal chr3q); 115 sandwich_same_hap and 1 complex. Thirteen of sixteen crossover-like events are in PAN028."
- **Divergence:** The talk says "ongoing and frequent." The paper reports raw counts (16 crossover-like, 133 gene-conversion-like) but never converts these to a per-meiosis per-arm rate. The WashU pedigree is 4 individuals across ~3 parent-child transmission events. A reviewer who works on meiotic recombination will immediately ask: what is the rate of inter-chromosomal exchange per meiosis per 100 kb of PHR? That number is computable from the data (16 crossover-like events / [number of parent-child transmission events] / [PHR length surveyed]) and is the quantity that tells you whether "frequent" is accurate or oversold. Without this rate estimate, "frequent" is an unmeasurable claim.
- **Concrete edit:** Add to P9 or to Methods §Pedigree: "The 16 crossover-like events across [N] parent-child haplotype pairs yield an estimated rate of [X] inter-chromosomal crossovers per meiosis per Mb of PHR (95% CI by Poisson: [L, U]). At the per-generation scale, this rate implies [Y] new inter-chromosomal crossovers per meiosis genome-wide across the 41 signal-bearing arms." Compute and insert. If rate computation is infeasible with current pedigree depth, state that explicitly rather than reporting counts without context.

---

### Finding 10 -- MINOR

**Paper invites question: the 9-fold B/W strengthening at flanking implies PHR-level B/W is inflated, not validated, by the control.**

- **Draft cite:** P8: "HG002 Hi-C PHR B/W = 0.027 falls to flanking B/W = 0.0031 (9-fold strengthening, Fig. 3d). Multi-mapping of reads to paralogous sequence cannot inflate a within-community signal in regions that contain no paralogous sequence."
- **Transcript cite:** NARRATIVE_EXTRACT §6 Q1: Eric acknowledges random MAPQ0 assignment as something "we probably should try to correct for." Transcript 00:45:36-00:45:44 (subtitle 339).
- **Verbatim transcript:** "So we've done something that we probably should try to correct for."
- **Verbatim draft:** "Multi-mapping of reads to paralogous sequence cannot inflate a within-community signal in regions that contain no paralogous sequence."
- **Divergence:** The paper uses the 9-fold strengthening at flanking as evidence against multi-mapping artefact. But the logical corollary -- never stated in the paper -- is that the 9-fold difference means PHR-level B/W (0.027) is 9x weaker than the true proximity signal at flanking (0.0031). Put differently: if the flanking unique-sequence B/W = 0.0031 represents "true" 3D proximity (no multi-mapping) and PHR B/W = 0.027 represents "multi-mapping-inflated" proximity, the multi-mapping is making the PHR signal look WEAKER, not stronger. Note (partial overlap with PEER_REVIEW_v1 M3): M3 flagged that Mantel ρ weakens at flanking while B/W strengthens, and the paper now reports both and explains why. But neither M3 nor v4 states the corollary that PHR B/W = 0.027 is the inflated (weaker-than-true) estimate, and flanking B/W = 0.0031 is the artefact-controlled (stronger-than-PHR) estimate. The paper's framing inverts the logical direction of the control.
- **Concrete edit:** Add one clause to P8: "The 9-fold strengthening at unique-sequence flanking (B/W = 0.0031 vs PHR B/W = 0.027) is consistent with multi-mapping distributing PHR inter-chromosomal reads across paralogous loci on both within- and between-community arms, reducing the apparent within-community enrichment at PHR loci relative to what is measured at unique flanking sequence. The flanking B/W therefore provides the lower-bound artefact-controlled estimate of true 3D co-localisation, not a ceiling."

---

## 3. Recommendation summary

The paper's narrative match to the Q&A is uneven: Q2 (chicken-or-egg) is addressed adequately in v4 (though the causal-scenarios gap in Finding 5 weakens it), Q1 is answered only on one of its two concerns and contains a critical reproducibility disclosure gap (Finding 1), and Q3 is entirely absent from the paper with no acknowledgement. The single most important fix is Finding 1: disclose explicitly in Methods and in Data Availability that the signal is invisible in standard deposited Hi-C files and requires complete BAM re-alignment with MAPQ0 reads retained. Every other gap is fixable with one sentence; Finding 1's omission makes the paper irreproducible from its cited data sources.

---

## 4. Audit trail

**Tools used:**
- Read (full transcript at absolute path `/home/guarracino/Dropbox/grants/r21-cancer-pangenomics/notes/Session7-PopulationGenomics.en.srt`, all 1829 lines)
- Read (NATURE_DRAFT_v4.md, all 119 lines)
- Read (NARRATIVE_EXTRACT.md, all 280 lines)
- Read (PEER_REVIEW_v1.md, all 188 lines)
- Read (REVISION_LOG_v4.md, all 93 lines)

**Lines of transcript read:** First timestamp 00:27:59,850 (subtitle 1, line 1) through last timestamp 00:48:36,640 (subtitle 395, line 1829). Q&A section: subtitles 326-395, approximately lines 1515-1829.

**Transcript Q&A timestamp ranges parsed:**
- Q1 (Hi-C validity): subtitles 327-349 (00:45:10-00:46:11)
- Q2 (chicken-or-egg): subtitles 350-368 (00:46:14-00:47:22); follow-up comment subtitles 369-371 (00:47:23-00:47:31)
- Q3 (P-to-Q arm): subtitles 372-392 (00:47:32-00:48:33)

**Lines of draft read:** All main text and Methods (119 lines of NATURE_DRAFT_v4.md). Specific sections consulted: Abstract, P2 (methods/PHR definition), P7 (3D validation opening), P8 (anti-artefact argument), P9 (pedigree), P12 (causal loop and limitations), Methods §Hi-C Pore-C and CiFi pipeline, Methods §Pedigree odgi-untangle.

**Sections of NARRATIVE_EXTRACT.md consulted:** §1 (metadata), §2 (narrative arc), §3 (verbatim gold G01-G18), §4 (in-talk-not-in-draft T01-T10), §6 (Q&A: Q1, Q2, Q3 in full), §7 (suggested upgrades U01-U10).

**Overlap with PEER_REVIEW_v1.md noted and avoided:** M3 (flanking paradox, partially overlaps Finding 10 -- I add only the logical corollary the peer reviewer missed), M4 (pedigree null -- not repeated), M7 (causal loop not closed in human -- not repeated, only the scenarios-not-distinguished gap in Finding 5 is new). Findings 1, 2, 3, 4, 6, 7, 8, 9 are all new gaps not in PEER_REVIEW_v1.md.
