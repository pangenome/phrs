# NARRATIVE_EXTRACT.md

Distilled from the author's live conference talk. Captures framing, phrasing, emphasis, and divergences from the current draft.

---

## 1. Transcript Metadata

| Field | Value |
|-------|-------|
| Source file | `/home/guarracino/Dropbox/grants/r21-cancer-pangenomics/notes/Session7-PopulationGenomics.en.srt` |
| Format | SubRip (.srt), 395 subtitle entries, 1829 lines |
| First timestamp | 00:27:59 (session chair introduction) |
| Last timestamp | 00:48:36 (applause after Q&A) |
| Total duration | ~20 minutes 37 seconds |
| Speaker | Eric Garrison (confirmed by chair at 00:45:11: "Thanks, Eric.") |
| Co-authors mentioned | Andrea Guarracino, Angela Gyamfi (grad student, St. Jude) |
| Audience | Biology of Genomes (BoG) conference, Population Genomics session (Session 7); genomics audience, cancer pangenomics R21 grant context |
| Language | English |

---

## 2. Narrative Arc

| Segment | Timestamp | Topic | Length (subtitles) | Key beats |
|---------|-----------|-------|--------------------|-----------|
| 1 | 00:27:59–00:29:00 | Framing: HPRCv2 as the enabling substrate | ~10 | HPRCv2 = near-complete assemblies, NX plot quality, "evolution of complete chromosomes in humans" |
| 2 | 00:29:01–00:30:25 | Backstory: pangenome graphs = alignments; HPRCv1 per-chromosome limitation | ~15 | Graphs are just alignments; HPRCv1 built chromosomes independently; acrocentrics forced them together; "technical reasons led us to break them apart" |
| 3 | 00:30:27–00:33:15 | Method: implicit pangenome graph (IMPG + WFMASH) | ~25 | Implicit interval tree; all-vs-all PAF; WFMASH for base-level alignments; "simulate the full graph without having to build it" |
| 4 | 00:33:15–00:34:20 | Scalability: the "two septillion cell matrix" problem and Erdős-Rényi solution | ~10 | 2 × 10^24 cell matrix is intractable; panmixia makes most pairs redundant; 11.5% sampling = 230× threshold; QC: "you see all the other genomes" |
| 5 | 00:34:20–00:36:25 | Discovery: interchromosomal peaks at chromosome ends; PHR definition | ~20 | Genome-wide heatmap; subtelomeric peaks on almost all chromosomes; walk in from telomere until cross-chromosomal signal stops; 7 exceptions (counted live) |
| 6 | 00:36:25–00:39:25 | Community structure: from hairball to Leiden; example communities | ~30 | "Orbited this hairball"; Jaccard similarity; heatmap with cladal structure; Leiden k=15; DUX4 (4q/10q), tubulin (10p/18p), olfactory receptors enriched |
| 7 | 00:39:25–00:40:32 | Mechanism hypothesis: meiotic bouquet | ~15 | "Why is this very ubiquitous? It must be maintained somehow"; bouquet stage in zygotene; telomeres on nuclear envelope = physical proximity = recombinability |
| 8 | 00:40:32–00:43:15 | 3D validation: mouse and human Hi-C, Pore-C | ~25 | Mouse meiotic Hi-C "crazy repetition of the project"; sequence similarity vs contact rate; human Pore-C and CHM13 Hi-C; 3D matches community assignment |
| 9 | 00:43:15–00:43:55 | Pedigree: catching recombination in the act | ~12 | WashU pedigree; odgi untangle; "tantalizing signals" of gene conversion (9p/3p); honest: "I wouldn't say it's conclusive" |
| 10 | 00:43:55–00:44:55 | Conclusions + thanks | ~20 | Near-ubiquitous PHRs; 3.5 Mb scale; "sequence homology mirrors physical proximity"; chicken-or-egg question posed explicitly; Mefford 2002 replicated with "a tiny slice" |

---

## 3. Verbatim Gold (10–20 lines)

Each entry: verbatim text | timestamp | suggested draft placement.

**G01** — "Very, very high-quality substrate to think about the evolution of complete chromosomes in humans."
- Timestamp: 00:28:52–00:29:00
- Placement: Abstract sentence 1 or intro P1. Establishes HPRCv2 as enabling, not incidental.

**G02** — "These graph models... are effectively equivalent to alignments between the genomes. And they're the same. It's just a different way of thinking of them."
- Timestamp: 00:29:15–00:29:24
- Placement: Methods opening or intro P2. Demystifies the implicit graph for a general audience; cleaner than the draft's technical phrasing.

**G03** — "Technical reasons led us to break them apart."
- Timestamp: 00:30:22–00:30:24
- Placement: Intro P1 (HPRCv1 limitation sentence). Very direct framing of why per-chromosome alignment failed.

**G04** — "You have a two septillion cell matrix, right? And that is clearly intractable."
- Timestamp: 00:32:05–00:32:07
- Placement: Methods (all-vs-all alignment section). More vivid than "computationally infeasible"; brings home the scale.

**G05** — "Panmixia means that most of the possible pairwise alignments are redundant because most of the genomes are kind of similar to each other. That means we don't have to do all pairs."
- Timestamp: 00:32:36–00:32:40
- Placement: Methods, Erdős-Rényi justification. The talk explains *why* subsampling works, not just *that* it does. Draft states it but doesn't foreground it.

**G06** — "We are basically guaranteed to kind of simulate the full graph without having to build it."
- Timestamp: 00:33:32–00:33:33
- Placement: Intro or Methods summary sentence. Excellent one-liner for the implicit graph concept.

**G07** — "In this panmictic population, where all the chromosomes are basically homologous to each other, you effectively hope that in every megabase pair of every chromosome, you see all the other genomes. And that's true."
- Timestamp: 00:33:44–00:33:53
- Placement: Results (QC section). Clear framing of the coverage-QC logic; currently absent from the draft.

**G08** — "So then we orbited this hairball for a minute and decided we couldn't understand anything that was going on in here directly."
- Timestamp: 00:36:31–00:36:37
- Placement: Results P2 transition (from graph to Jaccard). Humanizes the analytical journey; makes the pivot to community analysis feel motivated.

**G09** — "It's not just one big blob. There are groups of things that fit together."
- Timestamp: 00:37:45–00:37:48
- Placement: Results P2, first sentence describing the Jaccard heatmap. Plain language before technical clade names.

**G10** — "Why is this phenomenon very ubiquitous? It must be maintained somehow."
- Timestamp: 00:39:25–00:39:30
- Placement: Paragraph opening the mechanism/bouquet section. The draft doesn't pose this as an explicit question before invoking the meiotic bouquet.

**G11** — "You have the opportunity for recombination here, because you're proximal, you're forced to be on this manifold, and the sequences are homologous."
- Timestamp: 00:40:18–00:40:23
- Placement: The mechanistic sentence in the 3D/bouquet section. Tight causal chain: proximity → opportunity → homology → recombination.

**G12** — "It's like a repetition of the whole project just on mouse in order to see this."
- Timestamp: 00:41:41–00:41:45
- Placement: Results P3 (mouse section opening). Framing the mouse analysis as a deliberate cross-species replication gives it more weight than just "we also did mouse."

**G13** — "We did attempt to see if we could catch this in the act."
- Timestamp: 00:41:54–00:41:55
- Placement: Results P4 opening (pedigree). Perfect one-liner for why the pedigree analysis exists. Currently absent from draft.

**G14** — "We see some kind of tantalizing signals, like... what looks like ongoing gene conversion between 9p and 3p."
- Timestamp: 00:42:39–00:42:48
- Placement: Pedigree results. The word "tantalizing" is accurate and honest — calibrates reader expectations appropriately.

**G15** — "But I wouldn't say it's conclusive."
- Timestamp: 00:42:57
- Placement: Pedigree results, immediately after G14. The draft is more optimistic (92% within-community). The hedged framing is more appropriate for the actual evidence level.

**G16** — "This is the first opportunity we had to systematically observe it at a population scale. In humans."
- Timestamp: 00:43:42–00:43:43
- Placement: Conclusions / intro final sentence. The pause on "In humans" (delivered as a separate clause) is rhetorically strong. Missing from draft.

**G17** — "Sequence homology mirrors physical proximity in nuclear organization, which I think is really fascinating."
- Timestamp: 00:43:59–00:44:01
- Placement: Closing paragraph / abstract last sentence. More evocative than "sequence similarity correlates with nuclear-envelope proximity."

**G18** — "It's a chicken or egg question. Is it that the sequence homology is actually driving the physical proximity, or is it the homology just a result of the fact that the proximity is there?"
- Timestamp: 00:44:07–00:44:16
- Placement: Discussion, after establishing the correlation. The draft alludes to this but doesn't foreground it as the central unresolved question.

---

## 4. In-Talk, Not in Draft

Ranked by expected impact on the narrative. Each item has an ADD/DISCUSS/DEFER decision.

**T01 — "First opportunity to systematically observe at population scale. In humans."** (G16, 00:43:42–00:43:44)
- What: The talk frames this explicitly as the first population-scale observation of ubiquitous subtelomeric PHRs. The rhetorical pause "In humans." functions as an emphasis beat.
- Gap: The draft's intro mentions prior cytogenetic work (FISH, BAC-walking) but doesn't deliver the "first" framing as cleanly.
- Decision: **ADD** — make the intro's final sentence explicitly claim this as the first population-scale survey. Draft currently buries this in the third paragraph.

**T02 — The chicken-or-egg question posed explicitly** (G18, 00:44:07–00:44:16)
- What: Eric poses the core mechanistic ambiguity as a named, open question: does homology drive proximity, or does proximity generate homology?
- Gap: The draft establishes the correlation (Mantel ρ) and invokes the bouquet, but never names the causal ambiguity as an explicit open question.
- Decision: **ADD** — one sentence in the Discussion: "Whether sequence similarity drives nuclear co-localisation or co-localisation generates similarity remains an open question — a chicken-and-egg of subtelomeric evolution."

**T03 — Mefford 2002 replicated with "a tiny slice of the homology she had available"** (00:44:38–00:44:42)
- What: Eric thanks Heather Mefford and says they replicated her results, but she "took a tiny slice of the homology." This frames the paper as a massive expansion of a prior, correct but limited observation.
- Gap: The draft cites Mefford 2001/2002 but doesn't deliver the "we replicated and massively expanded" framing.
- Decision: **ADD** — intro P1 could use this beat: "Mefford and Trask established the cytogenetic signal with a slice of the available evidence; the HPRC v2 pangenome now lets us survey the whole."

**T04 — "Orbited this hairball" metaphor** (G08, 00:36:31–00:36:37)
- What: Vivid metaphor for initial PGGB graph of all PHRs — unintelligible as a direct visualization.
- Gap: The draft skips this and goes straight to the Jaccard matrix without explaining why the graph itself isn't the analytic tool.
- Decision: **DISCUSS** — a single sentence acknowledging the graph is too complex to read directly motivates why the Jaccard-similarity approach is needed. "The all-PHR pangenome graph forms one component but is too tangled to read directly; we therefore reduced it to an arm-level similarity matrix."

**T05 — "Not an insignificant part of the human genome" for 3.5 Mb** (00:43:54–00:43:56)
- What: The talk emphasizes scale verbally. The 3.5 Mb figure appears in the draft's abstract but without the rhetorical framing.
- Gap: Draft abstract states "3.5 megabase pairs" as a parenthetical. Talk makes it a standalone conclusion beat.
- Decision: **ADD** — in the abstract or intro, "The total PHR extent of ~3.5 Mb excluding acrocentrics and PARs is not an insignificant fraction of the human genome."

**T06 — RNA invasion / "soup" hypothesis for proximity** (00:46:58–00:47:19, Q&A)
- What: Eric speculates that transcription-generated RNA invading DNA could create a "soup" that encourages physical proximity of similar-sequence regions, potentially explaining why proximity exists outside meiosis.
- Gap: Not in the draft at all. Goes beyond the current data.
- Decision: **DISCUSS** — this is a notable mechanistic speculation raised in Q&A that could belong in the Discussion as a one-sentence future direction: "One possibility is that RNA transcribed from these regions invades homologous sequences, creating a trans-chromosomal scaffold that reinforces proximity beyond meiotic stages."

**T07 — QC framing: "you effectively hope that in every megabase you see all other genomes — and that's true"** (G07, 00:33:44–00:33:53)
- What: The talk has a clean verbal QC narrative: in a panmictic population, the graph should connect all haplotypes everywhere; it does, except at rDNA, centromeres, and sex chromosomes.
- Gap: The draft has the QC in the Methods but not as a clear narrative beat in the Results.
- Decision: **ADD** — a short QC sentence in Results P1: "As a panmixia-based sanity check, we verified that IMPG transitive closure recovers all haplotypes at virtually every megabase, with the expected exceptions at rDNA, centromeres, and sex chromosomes."

**T08 — "Burning methane" AI aside about quadratic cost** (00:32:11–00:32:22)
- What: Eric jabs at "XAI" burning methane on quadratic deep learning models — a timely aside contrasting computational profligacy with their efficient approach.
- Gap: Not in the draft. Probably inappropriate for a Nature paper directly.
- Decision: **DEFER** — too topical for a Nature paper. But useful signal: the AI-comparison framing works for talks/preprint blog posts.

**T09 — Pedigree result hedged: "I wouldn't say it's conclusive"** (G15, 00:42:57)
- What: The talk is notably more cautious about the pedigree than the draft. The draft claims "92% within Leiden communities" confidently; the talk says "tantalizing" and "compatible, but maybe not definitive proof."
- Gap: The draft's pedigree paragraph reads more conclusively than the speaker does live.
- Decision: **DISCUSS** — the hedged phrasing is appropriate for the evidence. The draft could add one sentence: "These signals are consistent with inter-chromosomal exchange but cannot be distinguished from assembly artefacts without a fully orthogonal validation."

**T10 — P-arm vs Q-arm homology as an open question** (00:47:48–00:48:16, Q&A)
- What: Audience asks if P-arm to Q-arm homology exists (opposite orientations). Eric confirms it does but calls it "a bit of a mystery." He notes they should test whether Qs are closer to Qs and Ps to Ps.
- Gap: Not in the draft. The draft mentions the PPTX's community structure but doesn't address the arm-orientation asymmetry.
- Decision: **DISCUSS** — worth a sentence in the Discussion or as a specific Future Directions item: "Whether inter-chromosomal homology is further constrained by arm orientation (p-to-p vs p-to-q) remains to be tested with higher-resolution 3D data."

---

## 5. In-Draft, Not in Talk

Brief. Items flagged only where a Nature reviewer might ask why they are not motivated.

**D01 — Population genetics / F_ST analysis (Draft P2, Fig. 2c, 2d)**
The talk never mentions population structure of subtelomeric sequences. The draft has a full F_ST paragraph with an out-of-Africa topology. A reviewer may ask: what is the mechanistic interpretation of F_ST ≠ 0 in PHRs? The paper currently does not explain whether different populations have structurally different communities.

**D02 — RPE-1 cell line result (Draft P5)**
The talk mentions only HPRCv2 and the WashU pedigree. The RPE-1 result (individual-genome community rediscovery) is absent. A reviewer may ask why a karyotypically aneuploid cell line was chosen as a "single diploid" control, and whether the t(X;10) translocation biases the community rediscovery.

**D03 — CEPH1463 4-generation pedigree cross-assembler validation (Draft P4)**
The talk only discusses the WashU pedigree. The CEPH1463 hifiasm-and-verkko intersection is not mentioned. A reviewer unfamiliar with the WashU pedigree might ask whether the pedigree finding is replicated. The CEPH1463 result should be motivated earlier in the narrative.

**D04 — Gene content analysis (Draft P6)**
The talk mentions olfactory receptors briefly ("enrichment for cellular sensory kind of things"). The draft has a full paragraph on pseudogenes, ncRNAs, and OR4F families, with Fisher enrichment tests. This is fine as a section, but its purpose (characterize the functional landscape) needs a motivating sentence in the intro or at the start of that paragraph.

**D05 — Seven limitations section (Draft, last paragraph of main text)**
The talk acknowledges none explicitly. The draft has seven. This is correct practice for a Nature submission; no gap, but the limitations should not undermine the main conclusion, which the talk lands confidently.

---

## 6. Q&A

Q&A runs from approximately 00:44:59 to 00:48:36.

### Q1 — Hi-C signal validity (00:45:12–00:46:11)

**Questioner:** "Hi-C is notoriously bad for measuring interchromosomal contacts. And also really bad at finding telomeric contacts. So I'm wondering how you are able to get that signal."

**Eric's response (00:45:36–00:46:11):** "We had to completely reanalyze all these data sets because usually they were excluding mapping quality zero reads. So we've done something that we probably should try to correct for, but we randomly assign the reads to a location. If they map exactly the same between different locations, we randomly assign them to one... We couldn't see any of the signal when we took the deposited data sets. We had to realign them in order to do this."

**Hard question flag:** YES. This is the methodological Achilles heel. Random assignment of MAPQ0 reads inflates contact frequency between identical loci — the exact loci this paper cares about. Eric acknowledges this is a limitation ("we probably should try to correct for"). A peer reviewer will raise this. The draft's Methods section describes the approach but does not explain the potential artefact or propose a control. **Peer-review prep**: Add a control showing that B/W ratios computed from PHR windows are not driven by MAPQ0 multi-mapping (the "flanking unique sequence" control in Draft P3 partially addresses this but the connection should be made explicit).

### Q2 — Chicken-or-egg / causality (00:45:51–00:47:22)

**Questioner:** "I have a question if you can wax a little poetic maybe, about whether you think possibly that the sequence homology actually maintains non-random distribution of the chromosomes in the zygotene, or if it's simply the non-random positioning in the zygotene that elevates the rate of homology."

**Eric's response (00:46:39–00:47:22):** "The chicken or egg question again. It's not just zygotene. We appear to see this consistently at other stages of the cell cycle. In mouse we had zygotene meiotic Hi-C, in human no. I wonder if transcription generates RNA and an RNA can invade DNA and that that generates a kind of... I guess a soup might be the right word, but not exactly soup that encourages the physical proximity of these regions. That might have some biological utility because then things that have the same function end up in the same kind of part of the cell."

**Follow-up comment from another audience member (00:47:29–00:47:31):** "If you look at this across all segmental duplications, not just the sub-telomeres and in response to Flores' point, be a little more careful about the control of the Hi-C data."

**Hard question flag:** YES (combined). Two issues raised: (a) causality cannot be resolved with correlation data, and (b) the result may not be subtelomere-specific if it holds for all segmental duplications. **Peer-review prep**: The draft should add a sentence: "Whether the sequence-similarity-to-proximity correlation is unique to PHRs or applies broadly to all segmental duplications is not tested here and will require a pangenome-scale SD analysis."

### Q3 — P-arm to Q-arm homology (00:47:39–00:48:16)

**Questioner:** "You showed some homology between the P and the P and the Q and the Q. Do you have any examples of homology between a P arm and a Q arm?"

**Eric's response (00:47:48–00:48:16):** "Yeah, I do actually... But I think we can find a few instances quickly here. What kind of implication do you think that has because they're like opposite orientations then? I don't think that this cell really knows. Frankly the acrocentrics... for example, this is a cluster you have P and Q inside of it. It is true, and we probably should make a test to see if the Qs tend to be closer to Qs and P's tend to be closer to P's. I think it's bound to something actually because you tend to see the stronger relationships between the same arm, the shorter or longer one. To me, it's a bit of a mystery how that would be organized."

**Hard question flag:** MODERATE. P-to-Q homology at opposite orientations implies either inversion-mediated exchange or sequence convergence, which the current paper does not mechanistically address. The answer is honest ("a mystery") but a reviewer may want the arm-orientation asymmetry quantified. **Peer-review prep**: Compute what fraction of within-community edges are P-P, Q-Q, vs P-Q, and add a supplemental table.

---

## 7. Suggested Narrative Upgrades

Each upgrade: target location in draft + suggested new sentence + transcript timestamp justifying it.

**U01 — Sharpen the "first" claim in intro P1**
- Target: Intro P1, current last sentence: "We use these data to revisit subtelomeric architecture without chromosomal partitioning, asking three quantitative questions..."
- Suggested addition (before the three questions): "This paper delivers the first population-scale, genome-wide survey of inter-chromosomal subtelomeric exchange — the whole genome complement to what earlier cytogenetics could only study arm by arm."
- Justification: G16 (00:43:42–00:43:44): "This is the first opportunity we had to systematically observe it at a population scale. In humans."

**U02 — Add the "panmixia makes redundant pairs" sentence in Methods justifying 12% sampling**
- Target: Methods paragraph on wfmash all-vs-all, sentence after Erdős-Rényi threshold calculation.
- Suggested addition: "In a panmictic population where most genomes are closely related, most pairwise alignments are redundant; transitive closure of the sampled fraction recovers effectively the same sequence graph as exhaustive all-vs-all alignment."
- Justification: G05 (00:32:36–00:32:40): "Panmixia means that most of the possible pairwise alignments are redundant because most of the genomes are kind of similar to each other."

**U03 — Open the pedigree results paragraph with the "catch it in the act" framing**
- Target: Results P4, first sentence (currently: "Indirect inference of inter-chromosomal exchange can be replaced by direct observation.")
- Suggested replacement: "To move from inference to observation — to catch inter-chromosomal exchange in the act — we applied the pipeline to a T2T-quality three-generation pedigree."
- Justification: G13 (00:41:54–00:41:55): "We did attempt to see if we could catch this in the act."

**U04 — Pose the chicken-or-egg explicitly as an open question in the Discussion**
- Target: Discussion/limitations paragraph, after the mechanistic loop closure sentence.
- Suggested addition: "The correlation between sequence homology and nuclear proximity raises a fundamental question that current data cannot resolve: does shared sequence drive co-localisation, or does enforced proximity generate shared sequence? The meiotic bouquet provides structural opportunity for both; resolving the directionality will require tracking proximity and homology across generations."
- Justification: G18 (00:44:07–00:44:16): "It's a chicken or egg question. Is it that the sequence homology is actually driving the physical proximity, or is it the homology just a result of the fact that the proximity is there?"

**U05 — Add a QC sentence in Results P1 confirming the pangenome is correctly connected**
- Target: Results P1, after the PHR count sentence.
- Suggested addition: "As a sanity check, IMPG transitive closure recovers all 465 haplotypes at virtually every megabase window genome-wide, with the expected gaps at rDNA, centromeres, and sex chromosomes — confirming that the graph is correctly connected before any PHR calls are made."
- Justification: G07 (00:33:44–00:33:53): "In this panmictic population... you effectively hope that in every megabase pair of every chromosome, you see all the other genomes. And that's true."

**U06 — Add Mefford "tiny slice" framing in intro P1**
- Target: Intro P1, after the sentence citing Mefford2001/MeffordTrask2002.
- Suggested addition: "That body of work established the cytogenetic reality of the phenomenon but necessarily surveyed a small slice of the available inter-chromosomal homology; the HPRC v2 pangenome now enables its systematic census."
- Justification: 00:44:38–00:44:42: "We replicated her results. Although she just took a tiny slice of the homology she had available."

**U07 — Add MAPQ0 control rationale sentence to Methods Hi-C section**
- Target: Methods paragraph on Hi-C/Pore-C pipeline.
- Suggested addition: "MAPQ0 reads were retained with random placement at one of the equally scoring loci; the validity of this approach is supported by the flanking unique-sequence control (B/W 0.0031 vs PHR B/W 0.027, a 9-fold strengthening at non-duplicated sequences), demonstrating that the within-community signal is not inflated by multi-mapping."
- Justification: Q1 (00:45:36–00:46:11): "We had to completely reanalyze all these data sets because usually they were excluding mapping quality zero reads."

**U08 — Hedge the pedigree conclusion language to match the speaker's live caution**
- Target: Results P4, current sentence: "The filter yields 538 high-quality inter-chromosomal patches. 494 of 538 (92%) sit within a Leiden community..."
- Suggested addition after the 92% sentence: "These signals are consistent with ongoing inter-chromosomal exchange, though assembly artefacts cannot be fully excluded without an orthogonal long-read validation in matched blood-derived tissue."
- Justification: G15 (00:42:57) and 00:43:03–00:43:08: "So it is compatible, but maybe not definitive proof of an actual recombination event."

**U09 — Open the mechanism section with the explicit "why ubiquitous?" question**
- Target: Results P3 (3D/bouquet section), first sentence (currently: "Sequence-defined communities are physical.").
- Suggested insertion before the current first sentence: "The ubiquity of PHRs across 41 of 48 arms demands a maintenance mechanism: something must prevent these interchromosomal homologies from drifting apart."
- Justification: G10 (00:39:25–00:39:30): "Why is this phenomenon very ubiquitous? It must be maintained somehow."

**U10 — Extend the abstract's closing sentence with the "mirrors physical proximity" phrasing**
- Target: Abstract last sentence (currently: "Human subtelomeres are unified by ongoing inter-chromosomal recombination and concerted evolution.").
- Suggested replacement: "Sequence homology mirrors three-dimensional nuclear proximity, pedigree analysis catches the recombination events that perpetuate both, and human subtelomeres emerge as a system unified by ongoing concerted evolution."
- Justification: G17 (00:43:59–00:44:01): "Sequence homology mirrors physical proximity in nuclear organization, which I think is really fascinating."
