# Voice review plan for `submission/paper.tex`

## Scope and source basis

This is a planning artifact for harmonizing the active manuscript with the
revised Nature-style abstract and the BoG-2026 talk voice. It does not prescribe
changes to `submission/paper.tex` itself.

Primary inputs:

- `submission/paper.tex`, active single-file manuscript.
- `paper_prep/synthesis/ABSTRACT_nature.md`, target compressed Nature abstract.
- `paper_prep/synthesis/ABSTRACT_BoG.md`, longer talk-aligned anchor.
- `submission/notes/Session7-PopulationGenomics.en.srt`, talk transcript and
  cadence source.

The manuscript already follows the BoG story order: method enabling
population-scale subtelomere comparison, genome-wide PHR extent, community
structure, representative community architecture, 3D proximity, and pedigree
exchange. Preserve that order unless a later review finds a hard internal
inconsistency.

## Reusable voice rubric

Apply this rubric to every section-review task.

1. Start from a biological problem, then reveal the technical move.
   The target abstract begins with the unresolved biology: subtelomeric
   relationships resisted systematic analysis. The talk similarly motivates
   the method by explaining that chromosome-by-chromosome graph building broke
   apart shared acrocentric sequence (`Session7-PopulationGenomics.en.srt:179`
   to `Session7-PopulationGenomics.en.srt:232`) before introducing the implicit
   data structure (`Session7-PopulationGenomics.en.srt:245` to
   `Session7-PopulationGenomics.en.srt:290`). Avoid method-first paragraphs
   unless the section is explicitly Methods.

2. Use direct, measured claims: "we construct", "we find", "we observe",
   "we ask", "we define", "we use".
   The abstract voice is active and compact (`ABSTRACT_nature.md:9`), not
   promotional. The talk cadence uses plain transitions such as "the next step
   is to try to understand first the scale" (`Session7-PopulationGenomics.en.srt:671`
   to `Session7-PopulationGenomics.en.srt:682`) and "to make this a little more
   quantitative" (`Session7-PopulationGenomics.en.srt:1126` to
   `Session7-PopulationGenomics.en.srt:1133`). Keep verbs concrete.

3. Treat figures as the argument spine.
   Each Results subsection should make one figure-scale claim, then give the
   evidence needed to read that figure. The talk moves slide by slide: length
   scale (`Session7-PopulationGenomics.en.srt:685` to
   `Session7-PopulationGenomics.en.srt:760`), graph/community structure
   (`Session7-PopulationGenomics.en.srt:767` to
   `Session7-PopulationGenomics.en.srt:938`), 3D evidence
   (`Session7-PopulationGenomics.en.srt:1031` to
   `Session7-PopulationGenomics.en.srt:1243`), and pedigree evidence
   (`Session7-PopulationGenomics.en.srt:1246` to
   `Session7-PopulationGenomics.en.srt:1357`). Do not create side quests for
   reviewer-era analyses that have no current figure.

4. Explain technical novelty only to the level needed for the claim.
   The target abstract defines the implicit pangenome graph in one appositive:
   "a reference-free, all-to-all alignment..." (`ABSTRACT_nature.md:9`). In the
   body, the method section can expand this, but Results should avoid overlong
   algorithmic digressions that delay the biological result.

5. Make known positives and new discoveries a paired pattern.
   The abstract pairs expected relationships with novel structure
   (`ABSTRACT_nature.md:9`). The Results should repeatedly use this move:
   known PAR/acrocentric/10p-18p/4q-10q systems validate the approach, while
   six-q-arm and broader community structure are the new finding.

6. Use community language as descriptive structure, not phylogeny.
   The talk explicitly says the tree-like ordering is "just for grouping" and
   does not claim a direct phylogenetic relationship
   (`Session7-PopulationGenomics.en.srt:849` to
   `Session7-PopulationGenomics.en.srt:870`). The manuscript should preserve
   that caution wherever it uses "clade", "tree", "phylogeny", "community", or
   "cladistic".

7. Keep causality framed as a model unless directly measured.
   The talk's 3D section frames meiosis as opportunity: telomeres are proximal
   on the nuclear envelope and homologous sequence has a chance to recombine
   (`Session7-PopulationGenomics.en.srt:1091` to
   `Session7-PopulationGenomics.en.srt:1123`). It closes with the chicken-or-egg
   uncertainty (`Session7-PopulationGenomics.en.srt:1439` to
   `Session7-PopulationGenomics.en.srt:1458`). Use "predicts", "is consistent
   with", "mirrors", "tracks", and "suggests" for causal direction unless the
   experiment directly observes the event.

8. Let caveats sharpen, not deflate, the claim.
   The talk says the pedigree signal is positive but not conclusive
   (`Session7-PopulationGenomics.en.srt:1335` to
   `Session7-PopulationGenomics.en.srt:1357`) and later says the pedigree
   "suggests" ongoing exchange (`Session7-PopulationGenomics.en.srt:1460` to
   `Session7-PopulationGenomics.en.srt:1468`). Manuscript caveats should be
   placed near the claim they qualify, then the section should restate the
   defensible take-home.

9. Prefer exact headline numbers where stable.
   Use 466 near-complete haplotype assemblies when describing the HPRC v2 plus
   CHM13-scale dataset in the abstract-level framing, and use 465 HPRC v2
   haplotypes plus CHM13 when describing the arm-flank census. Use 3.51 Mb in
   Results and figure legends; 3.5 Mb is acceptable only in broad abstract or
   concluding prose.

10. Keep the final synthesis compact and unresolved where the data are
    unresolved.
    The talk conclusion is a sequence of compact claims: near-ubiquitous PHRs,
    population-scale observation, about 3.5 Mb outside acrocentrics/PARs,
    sequence homology mirroring proximity, unresolved directionality, and
    pedigree evidence suggesting ongoing exchange
    (`Session7-PopulationGenomics.en.srt:1373` to
    `Session7-PopulationGenomics.en.srt:1468`). The manuscript close should
    preserve this order and caution level.

## Current voice match and divergence

### Strong matches

- The active abstract at `submission/paper.tex:52` to `submission/paper.tex:54`
  is already close to `ABSTRACT_nature.md:9`: compact problem statement,
  method, genome-wide result, communities, 3D evidence, pedigree evidence, and
  final synthesis in one paragraph.
- The introduction preserves the old-evidence-to-new-capability move. Lines
  `submission/paper.tex:60` to `submission/paper.tex:94` summarize
  cytogenetic/clone evidence and the unresolved genome-wide question, then
  `submission/paper.tex:96` to `submission/paper.tex:106` turns to CHM13 and
  HPRC v2.
- The Methods-enabling Results section has an effective talk-like first
  sentence: "We treat every haplotype as its own reference"
  (`submission/paper.tex:110`). This is concise and conceptually grounded.
- The PHR extent section follows the "known positives plus ubiquitous new
  signal" pattern: acrocentrics and PARs appear as positive controls, then the
  signal recurs at almost every chromosome end (`submission/paper.tex:138` to
  `submission/paper.tex:151`).
- The community section correctly warns that ordering is a grouping device, not
  a phylogenetic claim (`submission/paper.tex:171` to
  `submission/paper.tex:173`), matching the talk's caution.
- The 3D section starts from the biological model and then tests it
  (`submission/paper.tex:196` to `submission/paper.tex:204`), matching the talk
  transition into the meiotic bouquet.
- The pedigree section contains the right caution: "positive but not
  definitive" and assembly-artifact risk (`submission/paper.tex:240` to
  `submission/paper.tex:243`), consistent with the talk.
- The closing paragraph preserves the chicken-or-egg uncertainty
  (`submission/paper.tex:254` to `submission/paper.tex:257`) and therefore
  avoids overstating causality.

### Divergences and review risks

- Abstract-level directness may overstate pedigree evidence relative to the
  body. `submission/paper.tex:53` says "we directly observe inter-chromosomal
  exchange", while the body says "positive but not definitive"
  (`submission/paper.tex:240` to `submission/paper.tex:243`) and the talk says
  the pedigree "suggests" exchange (`Session7-PopulationGenomics.en.srt:1460`
  to `Session7-PopulationGenomics.en.srt:1468`). Review should decide whether
  "directly observe" is reserved for the untangle patches as observations, or
  softened to "find direct evidence consistent with".
- The conclusion may over-convert a model into mechanism. Lines
  `submission/paper.tex:245` to `submission/paper.tex:253` state a
  self-reinforcing loop and maintenance by ongoing recombination; lines
  `submission/paper.tex:254` to `submission/paper.tex:257` then reopen
  directionality. A rewrite should make the loop explicitly a model before the
  causal uncertainty sentence.
- The 3D section contains a possible mismatch between body-level confidence and
  Methods-level caveat. The body says the unique-sequence control "confirms" the
  signal is not a multi-mapping artefact (`submission/paper.tex:221` to
  `submission/paper.tex:225`), while Methods says random placement is an
  acknowledged limitation and PHR-window B/W ratios are an inflated upper bound
  (`submission/paper.tex:512` to `submission/paper.tex:526`). Harmonize to a
  precise claim: controls support the signal and argue against a simple
  multi-mapping artefact, but do not estimate the exact PHR-internal contact
  magnitude.
- Methods still include reviewer-era control material that the project
  instructions say was cut from body and Methods unless needed for consistency:
  exclusion controls (`submission/paper.tex:556` to `submission/paper.tex:565`),
  single-cell 3D controls (`submission/paper.tex:566` to
  `submission/paper.tex:576`), and limitations mentioning FST
  (`submission/paper.tex:636` to `submission/paper.tex:640`). These may be
  necessary if figures or body still reference them, but they are high-risk for
  narrative drift and should be audited before prose rewriting.
- Figure 5 legend references cross-assembler validation and RPE-1 positive
  control as "described in the main text" (`submission/paper.tex:403` to
  `submission/paper.tex:404`), but the body no longer describes CEPH1463 or
  RPE-1. This is an internal consistency problem; fix the legend or restore only
  the minimal needed cross-reference after deciding whether these analyses are
  part of the current figure story.
- The manuscript alternates "interchromosomal" and "inter-chromosomal":
  abstract `submission/paper.tex:53`, questions at `submission/paper.tex:103`,
  community text at `submission/paper.tex:164`, figure legends at
  `submission/paper.tex:296` and `submission/paper.tex:363`. Standardize before
  section tasks fan out.
- The abstract says 3.5 Mb (`submission/paper.tex:53`) while Results says
  3.51 Mb (`submission/paper.tex:149`) and the Fig. 1 legend says 3.510 Mb
  (`submission/paper.tex:306`). This is not wrong, but the style rule should be
  explicit so reviewers do not "fix" it inconsistently.
- `submission/paper.tex:415` says "233 HPRC v2 v1.1 individuals", which is
  awkward and likely should be checked against the intended HPRC v2 naming
  convention during the Methods review.

## Global narrative changes before section rewriting

Do these checks before assigning independent section rewrites.

1. Decide the paper-wide evidence-strength ladder.
   Recommended ladder:
   - "We find" for sequence homology, PHR extent, graph connectivity, and
     community structure.
   - "We observe" for measured contacts and observed untangle patch patterns.
   - "We find evidence for" or "suggests" for ongoing germline exchange in
     humans.
   - "The model predicts" or "is consistent with" for the meiotic-bouquet
     mechanism and directionality between proximity and homology.

2. Resolve the reviewer-era Methods drift.
   Audit Methods subsections `submission/paper.tex:556` to
   `submission/paper.tex:576` and the FST/CI sentence at
   `submission/paper.tex:636` to `submission/paper.tex:640`. Keep only material
   needed to support current Fig. 4, Extended Data Fig. 1, or explicit body
   claims. Do not reintroduce body discussion of within-community heterogeneity,
   popgen/FST, the 14-test 3D forest, CEPH1463, RPE-1, gene enrichment, or other
   reviewer-era analyses unless required to fix an internal inconsistency.

3. Align Fig. 5 legend with the current body.
   `submission/paper.tex:403` to `submission/paper.tex:404` currently promises
   main-text discussion of CEPH1463 and RPE-1. This should be resolved before a
   pedigree rewrite, because it affects whether the section remains a clean
   talk-aligned "catch this in the act" endpoint.

4. Standardize terminology and numerals in a pre-pass.
   Use the standardization table below. This avoids each fan-out reviewer making
   local choices about PHR, interchromosomal, 465/466, 3.51 Mb, communities, and
   evidence verbs.

5. Preserve the BoG order, but sharpen handoffs.
   The current order is good. Rewriting should mostly improve paragraph
   transitions: unresolved biology to method, method to genome-wide scale, scale
   to structure, structure to representative loci, structure to 3D proximity,
   and proximity to pedigree exchange.

## Proposed section-review fan-out

Each task should produce a marked-up review memo, not edited manuscript text,
unless a later implementation round explicitly authorizes edits. Every memo
should include concrete examples with `submission/paper.tex` line references,
proposed replacement direction for risky sentences, and any cross-section terms
that must be coordinated.

### Task 1: Abstract, title, keywords, and paper-level promise

- Scope: `submission/paper.tex:42` to `submission/paper.tex:56`.
- Inputs: `ABSTRACT_nature.md`, `ABSTRACT_BoG.md`, transcript conclusion
  `Session7-PopulationGenomics.en.srt:1373` to
  `Session7-PopulationGenomics.en.srt:1468`, standardization table below.
- Review questions: Does the active abstract match the Nature target voice?
  Does it overstate direct observation or causal mechanism? Are 466/465, 3.5/3.51
  Mb, PHR naming, and interchromosomal spelling handled consistently?
- Expected deliverable: A one-page abstract voice memo with sentence-by-sentence
  risk labels: keep, tighten, soften, or coordinate globally.

### Task 2: Introduction and framing questions

- Scope: `submission/paper.tex:60` to `submission/paper.tex:106`.
- Inputs: `ABSTRACT_BoG.md:5` to `ABSTRACT_BoG.md:7`, transcript background
  `Session7-PopulationGenomics.en.srt:37` to
  `Session7-PopulationGenomics.en.srt:232`, current Fig. 1 and Fig. 2 legends
  `submission/paper.tex:295` to `submission/paper.tex:331`.
- Review questions: Does the introduction move from cytogenetic history to the
  HPRC/CHM13 opportunity without becoming a review article? Are the three
  questions at `submission/paper.tex:101` to `submission/paper.tex:106` the
  right paper-wide contract?
- Expected deliverable: A framing memo identifying paragraphs that match the
  talk/Nature cadence and paragraphs that should be compressed, especially
  where historical detail delays the current study.

### Task 3: Implicit pangenome graph method as Results

- Scope: `submission/paper.tex:108` to `submission/paper.tex:134`.
- Inputs: transcript method setup `Session7-PopulationGenomics.en.srt:236` to
  `Session7-PopulationGenomics.en.srt:380`, Methods lines
  `submission/paper.tex:435` to `submission/paper.tex:465`,
  `ABSTRACT_nature.md:9`.
- Review questions: Does the Results section explain the implicit pangenome
  graph at the right level, or does it import too much Methods detail? Are the
  465 HPRC haplotypes plus CHM13 and 466-assembly framings consistent?
- Expected deliverable: A technical-voice memo specifying which sentences are
  essential for the Results claim and which should be left to Methods.

### Task 4: Genome-wide PHR extent and Fig. 1 argument

- Scope: `submission/paper.tex:136` to `submission/paper.tex:151` plus Fig. 1
  legend `submission/paper.tex:291` to `submission/paper.tex:309`.
- Inputs: transcript PHR-scale section `Session7-PopulationGenomics.en.srt:640`
  to `Session7-PopulationGenomics.en.srt:760`, `ABSTRACT_nature.md:9`.
- Review questions: Does the section clearly establish PHR ubiquity and scale
  without overstating biology beyond the 500 kb flank limit? Are 3.51 Mb,
  median/mean, PAR2 comparison, and seven no-signal arms used consistently?
- Expected deliverable: A Fig. 1 review memo with a preferred number style and
  any needed caveat language for flank-size truncation.

### Task 5: Communities and representative architectures

- Scope: communities `submission/paper.tex:153` to
  `submission/paper.tex:173`, representative loci `submission/paper.tex:175` to
  `submission/paper.tex:192`, Fig. 2 legend `submission/paper.tex:311` to
  `submission/paper.tex:332`, and Fig. 3 legend `submission/paper.tex:334` to
  `submission/paper.tex:354`.
- Inputs: transcript graph/community section
  `Session7-PopulationGenomics.en.srt:767` to
  `Session7-PopulationGenomics.en.srt:960`, `ABSTRACT_nature.md:9`.
- Review questions: Is "community" consistently defined as sequence-similarity
  grouping? Are "clade", "phylogeny", and "tree" appropriately qualified? Do
  the representative communities feel like examples of the global result rather
  than a disconnected locus catalog?
- Expected deliverable: A community-language memo with a list of terms to
  standardize and any sentences that imply unsupported phylogeny.

### Task 6: 3D proximity, meiotic bouquet, and Fig. 4/ED1 caution

- Scope: `submission/paper.tex:194` to `submission/paper.tex:225`, Fig. 4
  legend `submission/paper.tex:356` to `submission/paper.tex:388`, ED1 legend
  `submission/paper.tex:647` to `submission/paper.tex:665`, and Methods support
  `submission/paper.tex:510` to `submission/paper.tex:593`.
- Inputs: transcript bouquet/3D section
  `Session7-PopulationGenomics.en.srt:1031` to
  `Session7-PopulationGenomics.en.srt:1243`, Q&A caution
  `Session7-PopulationGenomics.en.srt:1527` to
  `Session7-PopulationGenomics.en.srt:1671`.
- Review questions: Does the section distinguish human somatic 3D evidence from
  mouse meiotic-stage evidence? Is multi-mapping caveat language aligned between
  Results and Methods? Are "predicts", "mirrors", "tracks", and "organized by"
  used at appropriate caution levels?
- Expected deliverable: A 3D evidence memo with a claim-strength ladder for
  every 3D sentence and a recommendation on whether reviewer-era controls remain
  Methods-only.

### Task 7: Pedigree exchange and final synthesis

- Scope: `submission/paper.tex:227` to `submission/paper.tex:261`, Fig. 5
  legend `submission/paper.tex:390` to `submission/paper.tex:406`, and Methods
  support `submission/paper.tex:595` to `submission/paper.tex:608`.
- Inputs: transcript pedigree and conclusion
  `Session7-PopulationGenomics.en.srt:1246` to
  `Session7-PopulationGenomics.en.srt:1468`, `ABSTRACT_nature.md:9`.
- Review questions: Does the section preserve the "catch this in the act"
  narrative while retaining the talk's non-conclusive caution? Does the final
  paragraph present the self-reinforcing loop as a model rather than proven
  mechanism? Does the Fig. 5 legend incorrectly refer to analyses no longer in
  the body?
- Expected deliverable: A pedigree/conclusion memo with specific softening or
  tightening recommendations for "directly observe", "ongoing exchange",
  "maintained by", and "confirmation".

### Task 8: Methods and back-matter consistency audit

- Scope: all Methods and back matter, `submission/paper.tex:408` to
  `submission/paper.tex:667`.
- Inputs: all task memos above, project instruction that reviewer-era analyses
  without figures were intentionally cut from the body and Methods, current
  figure legend references.
- Review questions: Which Methods subsections are required by current figures
  and body claims? Which Methods lines preserve cut reviewer-era analyses? Are
  figure legends self-contained and consistent with the body?
- Expected deliverable: A Methods consistency memo listing keep/remove/relocate
  candidates with line references. This task should run after Tasks 1-7, because
  it integrates their claim-support needs.

## Standard terms, claims, and caution levels

| Topic | Standard | Avoid or check | Rationale and examples |
| --- | --- | --- | --- |
| PHR term | Use "pseudo-homolog region (PHR)" on first definition, then "PHR". In broad prose, "subtelomeric PHRs" is preferred. | Do not alternate with "pseudohomologous region" unless quoting or harmonizing historical literature. | Current definition is `submission/paper.tex:129` to `submission/paper.tex:131`; keywords use "pseudo-homolog region" at `submission/paper.tex:56`. |
| Interchromosomal spelling | Prefer "interchromosomal" as the paper-wide adjective. Use "inter-chromosomal" only if the journal style or an existing compound construction requires it. | Mixed spelling within abstract/body/legends. | Current mixed forms appear at `submission/paper.tex:53`, `submission/paper.tex:103`, `submission/paper.tex:164`, `submission/paper.tex:296`, and `submission/paper.tex:363`. |
| 466 vs 465 | Use "466 near-complete haplotype assemblies" for abstract-level HPRC v2 plus CHM13 scale. Use "465 HPRC v2 haplotypes plus CHM13" for the arm-flank census and figure legends. | Do not imply 466 HPRC haplotypes if CHM13 is the extra anchor; check "233 HPRC v2 v1.1 individuals" wording. | Active abstract uses 466 at `submission/paper.tex:53`; Results and Methods use 465 plus CHM13 at `submission/paper.tex:111` and `submission/paper.tex:415` to `submission/paper.tex:418`. |
| 3.51 Mb | Use 3.51 Mb in Results and legends; 3.5 Mb is acceptable in abstract/conclusion as a rounded reader-facing value. | Do not mix 3.5, 3.51, and 3.510 within the same section unless the precision difference is intentional. | Results: `submission/paper.tex:149`; Fig. 1 legend: `submission/paper.tex:306`; abstract and close: `submission/paper.tex:53`, `submission/paper.tex:250`. |
| Community language | "Sequence-similarity communities", "arm-level communities", "Leiden communities", or "community block structure". | "Phylogeny" without scare quotes or qualification; "clade" where it implies species-like descent of chromosome arms. | Manuscript correctly qualifies ordering at `submission/paper.tex:171` to `submission/paper.tex:173`; Fig. 2 legend uses "phylogeny" in quotes at `submission/paper.tex:323` to `submission/paper.tex:324`. |
| Direct observe vs suggest | "We observe patch patterns"; "the pedigree provides direct evidence consistent with ongoing exchange"; "suggests ongoing exchange" for the biological conclusion. | "We directly observe interchromosomal exchange" if used without the assembly-artifact caveat. | Abstract `submission/paper.tex:53` is stronger than body caution at `submission/paper.tex:240` to `submission/paper.tex:243` and transcript caution at `Session7-PopulationGenomics.en.srt:1335` to `Session7-PopulationGenomics.en.srt:1357`. |
| 3D evidence | "Sequence similarity predicts/mirrors/tracks 3D proximity"; "human somatic contact increases with PHR similarity"; "mouse meiotic Hi-C peaks at zygotene." | "Organized by meiosis" or "meiotic architecture maintains PHRs" as a proven human mechanism. | Human data are somatic in the body at `submission/paper.tex:211` to `submission/paper.tex:220`; mouse supplies bouquet-stage data at `submission/paper.tex:205` to `submission/paper.tex:210`; limitations note this at `submission/paper.tex:254` to `submission/paper.tex:257`. |
| Pedigree evidence | "538 high-quality patches; 494/538 (92%) within sequence communities; positive but not definitive." | Presenting every patch as a validated germline recombination event. | Body lines `submission/paper.tex:233` to `submission/paper.tex:243`; Methods lines `submission/paper.tex:595` to `submission/paper.tex:608`. |
| Known vs novel relationships | Pair "recover known relationships" with "reveal novel groupings" in the same paragraph where possible. | Listing novel communities before validating positives. | Abstract target uses this pattern in `ABSTRACT_nature.md:9`; manuscript uses it at `submission/paper.tex:164` to `submission/paper.tex:170`. |
| Reviewer-era analyses | Keep out of body and Methods unless needed for current figures or consistency. | Reintroducing FST, CEPH1463, RPE-1, gene enrichment, 14-test forests, or within-community heterogeneity as narrative claims. | Potential drift appears in Methods/back matter at `submission/paper.tex:556` to `submission/paper.tex:576`, `submission/paper.tex:636` to `submission/paper.tex:640`, and Fig. 5 legend `submission/paper.tex:403` to `submission/paper.tex:404`. |

## High-risk language and claim mismatches to resolve

1. Pedigree directness:
   - Current strong version: `submission/paper.tex:53`, "we directly observe
     inter-chromosomal exchange".
   - Current cautious version: `submission/paper.tex:240` to
     `submission/paper.tex:243`, "positive but not definitive".
   - Talk target: `Session7-PopulationGenomics.en.srt:1335` to
     `Session7-PopulationGenomics.en.srt:1357`, compatible but not definitive;
     `Session7-PopulationGenomics.en.srt:1460` to
     `Session7-PopulationGenomics.en.srt:1468`, pedigree suggests ongoing
     exchange.
   - Review action: choose a single paper-wide phrasing ladder before abstract
     and conclusion edits.

2. Mechanism and maintenance:
   - Current claim: `submission/paper.tex:245` to `submission/paper.tex:253`
     describes a self-reinforcing loop and says PHRs are maintained by ongoing
     non-allelic homologous recombination.
   - Current caveat: `submission/paper.tex:254` to `submission/paper.tex:259`
     leaves causal direction unresolved and notes human 3D data are somatic.
   - Review action: frame the loop as a model supported by sequence, 3D, and
     pedigree evidence; reserve "maintained by" for a qualified final sentence.

3. 3D control strength:
   - Current body: `submission/paper.tex:221` to `submission/paper.tex:225`
     says the control confirms within-community contact is not a multi-mapping
     artefact.
   - Current Methods: `submission/paper.tex:512` to `submission/paper.tex:526`
     says random placement is a limitation and PHR-window B/W ratios are an
     inflated upper bound.
   - Review action: align wording to "supports", "argues against", or "bounds"
     rather than "confirms" if the Methods limitation remains.

4. Fig. 5 legend promises absent text:
   - Current legend: `submission/paper.tex:403` to `submission/paper.tex:404`
     says CEPH1463 and RPE-1 are described in the main text.
   - Current body: no matching CEPH1463/RPE-1 discussion in Results.
   - Review action: remove that promise or add only the minimal text necessary
     if the figure truly requires it; do not revive the old reviewer-era
     narrative.

5. Methods cut-list drift:
   - Current Methods include exclusion controls, single-cell 3D controls, and
     an FST CI mention (`submission/paper.tex:556` to
     `submission/paper.tex:576`; `submission/paper.tex:636` to
     `submission/paper.tex:640`).
   - Review action: audit these against current figures before any voice polish,
     because polishing text that should be removed wastes effort and risks
     pulling the paper away from the BoG story.

6. Precision and nomenclature:
   - Current examples: 3.5/3.51/3.510 Mb precision split
     (`submission/paper.tex:53`, `submission/paper.tex:149`,
     `submission/paper.tex:306`); 466 vs 465 plus CHM13
     (`submission/paper.tex:53`, `submission/paper.tex:111`,
     `submission/paper.tex:296`); interchromosomal vs inter-chromosomal
     (`submission/paper.tex:53`, `submission/paper.tex:103`,
     `submission/paper.tex:164`, `submission/paper.tex:296`,
     `submission/paper.tex:363`).
   - Review action: run one terminology pre-pass before fan-out, or require the
     integration reviewer to normalize these after fan-out.

## Recommended execution order

1. Run Task 1 and Task 8 as gating reviews: Task 1 fixes the paper-wide promise;
   Task 8 identifies any text that should not be polished because it no longer
   belongs in the manuscript.
2. Run Tasks 2-7 in parallel after the terminology/evidence ladder is agreed.
3. Join the memos in an integration pass that checks cross-section consistency:
   abstract promise, introduction questions, figure spine, evidence verbs,
   terminology, and Methods support.
4. Only after the integration pass, open an implementation round to edit
   `submission/paper.tex` as a single LaTeX file.
