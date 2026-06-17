# Voice review 06: 3D proximity and meiotic bouquet

## Scope and voice target

Reviewed `submission/paper.tex` lines 194-225, 356-388, 510-593, and
647-665 against the house abstract voice in
`paper_prep/synthesis/ABSTRACT_TEXTURE_SYNTHESIS.md`. The relevant house rule is
to keep the human and mouse evidence in separate lanes: human Pore-C/Hi-C shows
that sequence-similar subtelomeres contact one another more often in three
dimensions, while mouse meiotic Hi-C shows that this relationship peaks at
zygotene, when telomeres cluster in the bouquet. The section should support a
model, not promote nuclear architecture into a demonstrated maintenance
mechanism.

In current `paper.tex`, Fig. 4c is the mouse meiotic-bouquet panel and ED1 is
the CHM13 human Hi-C companion to Fig. 4a (`submission/paper.tex:647-665`).
That distinction is correct in the scoped manuscript text and should be
preserved.

## 1. 3D evidence claim-strength ladder

| Line(s) | Sentence or caption claim | Evidence lane | Current strength | Voice recommendation |
| --- | --- | --- | --- | --- |
| 196-201 | "Why should pseudo-homolog regions persist..." through "opportunity for ectopic recombination." | Mechanistic model premise from meiotic biology | Appropriate if read as motivation. "Would bring" and "provide the opportunity" are model verbs, not direct observations. | Keep the opportunity framing. Consider "could bring" only if the surrounding section is made more conservative; otherwise "would" is acceptable because it describes the bouquet geometry rather than a measured exchange event. |
| 202-204 | "The model predicts that arms sharing a sequence community should make more three-dimensional contact, most strongly at the bouquet stage." | Model-level prediction | Mostly right. "Predicts" is acceptable because the subject is "the model." "Should make" slightly implies contact causation rather than observed association. | Prefer "The model predicts that arms sharing a sequence community should show more three-dimensional contact, with the strongest relationship at the bouquet stage." |
| 205-210 | "The cleanest test is in mouse..." through "only direct bouquet-stage measurement..." | Mouse meiotic timing | Strong but defensible if "direct" modifies bouquet-stage 3D measurement, not recombination. It correctly does not call the mouse result human evidence. | Keep mouse as the germline-stage timing test. Avoid "cleanest test" if it reads as the cleanest test of the whole human model; use "The direct timing test is in mouse..." or "The stage-resolved test is in mouse..." |
| 211-217 | "In human, three-dimensional contact likewise increases..." through "recovers the community block structure directly." | Human somatic/contact-map association | The evidence lane is right: HG002 Pore-C and CHM13 Hi-C are per-PHR-pair human contact associations. "Likewise" usefully links direction of association, but can blur mouse meiotic timing into human somatic maps. "Recovers ... directly" is a little strong but acceptable for a matrix-ordering observation. | Use "In human somatic/contact-map data, three-dimensional contact also increases..." if space allows. Replace "recovers ... directly" with "shows" or "reveals" only if the caption keeps the exact B/W statistic; "shows" is plainer and safer. |
| 218-220 | "The same within-community enrichment holds..." | Human and human-derived single-cell contact maps | Appropriate as supporting replication if kept separate from the mouse timing claim. | Keep in Results only if these datasets remain represented in Methods and are not reviewer-era side quests. Otherwise move to Methods-only as supporting contact-map replication. |
| 221-225 | "Because subtelomeric PHRs are paralogous..." through "confirms ... not a multi-mapping artefact." | Mapping caveat/control | Too strong relative to Methods lines 514-523. The Methods say random placement is an acknowledged limitation, PHR-window B/W is an inflated upper bound, and the flanking control does not bound all within-community MAPQ0 bias. | Replace "confirms" with "argues against a simple multi-mapping artefact" or "supports the signal while bounding the mapping concern." Do not say the control fully refutes all multi-mapping bias in Results. |
| 361 | Fig. 4 title: "Sequence similarity predicts three-dimensional proximity." | Figure-level synthesis | "Predicts" is strong but acceptable as a statistical figure title if the text repeatedly says "association" and "supports a model." Risk: readers may infer directionality from sequence to architecture. | Safer title: "Sequence similarity tracks three-dimensional proximity" or "Sequence similarity is associated with three-dimensional proximity." If retaining "predicts," avoid other causal verbs nearby. |
| 362-372 | Fig. 4a caption: HG002 Pore-C scatter and CHM13 ED1 note. | Human Pore-C association; ED1 human Hi-C replicate | Good. "Make more contact" and "contact rises with similarity" are measured association verbs. The caption correctly says each dot is one inter-chromosomal PHR pair within one genome. | Keep. If revising globally, use "contact one another more often" to match the abstract exactly. |
| 373-376 | Fig. 4b caption: community-ordered Pore-C matrix. | Human Pore-C community structure | "Within-community blocks dominate" is a strong visual/statistical claim but grounded by B/W and p-value. | Keep if B/W direction is documented for readers. Consider "within-community contacts are enriched" if the B/W metric could be misread because lower values are stronger elsewhere in Methods. |
| 377-386 | Fig. 4c caption: mouse meiotic Hi-C and zygotene peak. | Mouse meiotic timing | Correctly distinguishes non-human mammal and stage-resolved meiotic Hi-C. "Same no-averaging test as panel (a)" is useful but risks making mouse and human lanes sound identical. | Keep the mouse lane explicit: "the analogous per-pair test." Retain "peaks at the zygotene bouquet" because this is the timing evidence the house abstract needs. |
| 655-663 | ED1 caption: CHM13 Hi-C companion to Fig. 4a. | Human Hi-C replicate | Appropriate. ED1 is human CHM13 Hi-C, not mouse. "Contact rises with similarity" is measured and restrained. | Keep. If manuscript-level ED numbering changes, protect this distinction: CHM13 Hi-C replicate is human 3D association; mouse zygotene timing is Fig. 4c in the current draft. |

## 2. Results/Methods consistency issues

1. **Control strength mismatch.** Results lines 221-225 say the unique-sequence
   control "confirms" the within-community contact is not a multi-mapping
   artefact. Methods lines 514-523 say random placement is an acknowledged
   limitation, the flanking control does not bound every bias from uniform MAPQ0
   distribution across paralogous arms, and PHR-window B/W should be read as an
   inflated upper bound on the artefact-controlled signal. These cannot both be
   true at the same claim strength. Align the Results to "argues against a
   simple multi-mapping artefact" or "supports the contact association after a
   unique-sequence control," leaving the upper-bound language in Methods.

2. **"Refutes" and "falsifying" are too absolute in Methods.** Methods line
   517 says the flanking control "refutes the multi-mapping artefact in
   PHR-internal regions," while lines 518-519 immediately say it does not bound
   possible MAPQ0 bias within the same community. Lines 524-526 then say the
   flanking control "falsif[ies] multi-mapping-driven inflation." The internal
   caveat is more nuanced than those verbs. Suggested Methods strength:
   "argues against broad multi-mapping-driven inflation" and "bounds the
   direction of the concern but not the exact PHR-internal magnitude."

3. **Human versus mouse lane is mostly correct but could be sharper.** Results
   lines 205-210 correctly place the direct bouquet-stage timing in mouse, and
   lines 211-220 correctly place human evidence in HG002 Pore-C, CHM13 Hi-C,
   CiFi, Dip-C, and sperm scHi-C contact maps. The risk is the transition
   "In human ... likewise" at line 211, which may imply the same meiotic timing
   evidence exists in human. Prefer "In human contact-map data..." or "In human
   somatic and gamete-derived contact maps..." depending on whether sperm scHi-C
   remains in the Results sentence.

4. **Fig. 4 title overstates directionality unless balanced by prose.** The
   title "Sequence similarity predicts three-dimensional proximity" at line 361
   is not wrong statistically, but it is the strongest verb in the figure. If
   the Results keep "model predicts" at lines 202-204, the figure title should
   probably use "tracks" or "is associated with" to avoid a chain of predictive
   language that sounds mechanistic.

5. **"Dominates" and B/W direction need reader help.** Fig. 4b says
   "within-community blocks dominate" with `B/W = 0.056` at lines 373-376, while
   Methods lines 520-526 discuss PHR and flanking `B/W` values where smaller
   numbers indicate stronger within-community enrichment. This is technically
   tractable but not self-explanatory. If space allows, use "within-community
   contacts are enriched" and let Methods define the ratio direction.

6. **Stage-resolution caveat belongs near the mouse claim or in Methods.**
   Methods lines 586-591 say the per-pair zygotene peak is resolution-dependent
   and is clearest at 10-20 kb, while the Results and caption give the zygotene
   peak without the caveat (`submission/paper.tex:205-210`, `377-386`). This is
   acceptable if kept as a Methods caveat, but the Results should avoid
   "cleanest test" unless the resolution dependence remains visibly bounded
   somewhere nearby.

7. **ED1 identity is consistent in scoped text.** Lines 647-665 define ED1 as
   the CHM13 Hi-C replicate of Fig. 4a. This is consistent with line 215 and the
   Fig. 4 caption line 372. Do not move mouse meiotic timing into ED1 unless the
   figure plan is changed across the manuscript; in this draft, mouse belongs to
   Fig. 4c.

## 3. Reviewer-era controls that should remain Methods-only

Keep these controls Methods-only unless a current figure panel or Results
sentence explicitly requires them:

- **Multi-resolution exclusion controls** at lines 556-565. They are useful as a
  robustness audit against acrocentric/sex/strong-community confounds, but they
  are not part of the current Fig. 4 spine. They should not be promoted into
  Results unless a figure is restored. In Methods, describe them as controls
  that bound confounding rather than as additional confirmation of the model.

- **Single-cell 3D controls** at lines 566-576. Dip-C, sperm scHi-C, PBMC Dip-C,
  and the zero-signal-arm control are supporting contact-map checks. Results
  lines 218-220 already mention Dip-C and sperm scHi-C; if the body is being
  compressed to the house abstract voice, keep the extra cell/control details in
  Methods and let Results say only that the enrichment is reproduced in
  additional maps. Avoid a Results detour into PBMC or `S_all` unless a panel is
  shown.

- **Strict-MAPQ and flanking unique-sequence control details** at lines 524-531.
  The Results need only the interpretive caveat: multi-mapping was retained and
  unique-sequence controls argue against a simple artefact. The exact
  `.allValidPairs`, MAPQ >= 30, Poisson-noise, and observed-over-expected
  normalisation details belong in Methods.

- **Mouse arm-level Mantel cross-check** at lines 587-593. This is a useful
  robustness check for the zygotene peak and can remain in Methods. Results and
  Fig. 4c should stay with the per-PHR-pair test because it is the direct
  parallel to the human scatter and avoids adding another metric layer.

Do not add reviewer-era forest plots, RPE-1/CEPH-style validation, FST, gene
enrichment, or 14-test 3D control language to this Results section unless those
figures return. The house voice wants the 3D claim to be a compact evidence
ladder, not a catalogue of every robustness analysis.

## 4. Phrase-level risks and suggested alternatives

| Current phrase | Line(s) | Risk | Suggested alternative |
| --- | --- | --- | --- |
| "persist" | 196 | Fine as the motivating biological question, but could imply the section will prove maintenance. | Keep, or use "Why are pseudo-homolog regions found at 41 of the 48 arms?" if the rewrite wants less mechanism. |
| "would bring ... and provide the opportunity" | 200-201 | Acceptable model framing; not a measured claim. | "could bring ... and create recombination opportunity" if the paragraph is softened. |
| "predicts" | 202, 361 | Good when the subject is "the model"; stronger as a figure title. | Results: keep "The model predicts." Figure title: "tracks" or "is associated with." |
| "should make more three-dimensional contact" | 202-203 | "Make" can sound causal. | "should show more three-dimensional contact." |
| "most strongly at the bouquet stage" | 203 | Good for mouse; not available in human. | "with the strongest relationship expected at the bouquet stage." |
| "The cleanest test is in mouse" | 205 | Could imply mouse is the cleanest test of the whole human model, not only the timing prediction. | "The stage-resolved test is in mouse" or "The direct bouquet-stage test is in mouse." |
| "only direct bouquet-stage measurement" | 209 | Safe only if "direct" modifies 3D measurement, not recombination. | "only stage-resolved bouquet 3D measurement of this relationship..." |
| "In human ... likewise" | 211 | Blurs human somatic/contact-map association with mouse meiotic timing. | "In human contact-map data, three-dimensional contact also increases..." |
| "recovers the community block structure directly" | 216-217 | "Directly" can overstate because ordering by community is imposed before visual recovery. | "shows the same community block structure" or "reveals within-community blocks after ordering by sequence community." |
| "The same within-community enrichment holds" | 218 | Fine if the datasets stay as support; a little broad if controls are caveated. | "The same within-community enrichment is observed..." |
| "confirms ... not a multi-mapping artefact" | 223-225 | Too strong against Methods caveats. | "argues against a simple multi-mapping artefact" or "supports the signal after a unique-sequence control." |
| "refutes" | 517 | Contradicts the following limitation at lines 518-519. | "argues against" or "addresses the broad concern of." |
| "inflated upper bound" | 522 | Strong and useful caveat, but should be mirrored by softened Results language. | Keep in Methods; align Results to "upper-bound/controlled signal" rather than "confirmed true signal." |
| "falsifying multi-mapping-driven inflation" | 526 | Too absolute given acknowledged random-placement limitation. | "arguing against multi-mapping-driven inflation as the sole explanation." |
| "within-community blocks dominate" | 376 | Strong visual language; B/W metric may confuse. | "within-community contacts are enriched." |
| "the same no-averaging test" | 382 | Useful, but can over-equate mouse timing and human contact-map evidence. | "the analogous per-pair, no-averaging test." |
| "peaks at the zygotene bouquet" | 384-386, 585-590 | Correct and central; keep, with Methods caveat about resolution dependence. | Keep. |

## Bottom-line recommendation

The section already has the right architecture: model premise, mouse meiotic
timing test, human Pore-C/Hi-C association, contact-map replication, and mapping
caveat. The main revision need is calibration, not restructuring. Soften
"confirms/refutes/falsifying" to "argues against/supports/bounds," keep random
placement as an acknowledged Methods limitation, and prevent the phrase
"predicts three-dimensional proximity" from becoming a causal claim that human
somatic contact maps prove meiotic maintenance. The defensible house-voice
take-home is: human contact maps show that sequence-similar subtelomeres contact
more often, mouse meiotic Hi-C shows the relationship is strongest at zygotene,
and together these observations support a bouquet/recombination-opportunity
model for subtelomeric concerted evolution.
