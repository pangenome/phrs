# B5 3D Apparatus Essentiality Audit

Date: 2026-06-17
Task: `manuscript-revision-b5`
Scope: judgment-support record only. I did not edit `submission/paper.tex`,
figures, captions, bibliography, or analysis outputs.

## Question

Apply the prompt essentiality test to the remaining 3D-control apparatus:

- five-set x five-resolution exclusion Mantel walk;
- W/B bootstrap;
- Mann-Whitney global test;
- observed/expected normalization;
- per-bin-pair normalization;
- any named normalization schemes currently attached to the 3D-contact claim.

Essentiality test used here: does this apparatus protect a current abstract-level
claim that is not already protected by the pointwise human contact evidence plus
the flanking unique-sequence control? If yes, keep it visible. If it only answers
a likely reviewer challenge, demote it to a compact Methods/Supplement support
role. If it is redundant, source-incomplete, or makes the manuscript look more
fragile than the evidence requires, cut unless an author explicitly chooses to
retain it.

## Abstract Claims Under Review

To avoid confusion with Cluster A task IDs, this record defines the relevant
abstract claims as follows:

| Claim ID | Current abstract claim | What must be protected |
|---|---|---|
| A4 | "Human Pore-C and Hi-C show that sequence-similar subtelomeres contact one another more often in three dimensions" | A positive human sequence-to-contact association, not a proof that every PHR-internal read is uniquely placed. |
| A5 | "Together, these data ... support a model in which nuclear architecture and recombination opportunity contribute to their concerted evolution" | A conservative multi-evidence model claim, where 3D contact is one line of support alongside sequence communities, mouse meiosis, and pedigree evidence. |

The current lead evidence for A4 is the community-free, per-sequence-pair human
Spearman family: HG002 Pore-C rho = 0.381, CHM13 Hi-C rho = 0.716, HG002 Hi-C
rho = 0.662, and HG002 CiFi rho = 0.191. The current primary defense against
multi-mapping is the flanking unique-sequence control. This audit assumes those
two layers remain available, as inventoried in `B0_3d_inventory.md`.

## Executive Recommendation

Keep only three 3D apparatus layers as load-bearing for the active manuscript:

1. **Pointwise per-PHR-pair contact normalization** for Fig. 4A / ED1, because
   it defines the lead A4 measurement.
2. **A concise MAPQ0/random-placement caveat plus flanking unique-sequence
   control**, because this is the essential anti-artifact defense for A4.
3. **A compact observed/expected statement for the community-ordered Pore-C
   matrix**, because Fig. 4B is a matrix visualization and O/E answers the
   obvious marginal-contact/chromosome-size objection.

Demote the exclusion Mantel walk, W/B bootstrap, and broad multi-resolution
apparatus to Supplement/Methods support. They are useful if a reviewer asks
whether the result is driven by acrocentric/sex/PAR communities, bin size, or
community labels, but they are not necessary to state A4 or A5 once the
community-free pointwise evidence and flanking control are visible.

Cut the Mann-Whitney global test from the manuscript-facing statistical
apparatus unless an author wants it retained as legacy provenance for a
specific table. It does not protect an abstract claim better than the pointwise
Spearman, W/B effect size, Mantel matrix check, and flanking control.

## Classification Table

| Apparatus item | Recommendation | Abstract claim protected | Does the flanking control already close the main threat? | Rationale | Author decision required |
|---|---|---|---|---|---|
| Five-set x five-resolution exclusion Mantel walk | **Demote** to compact Supplement/Methods support; do not lead with it | A4 indirectly; A5 only as robustness context | **Mostly yes** for multi-mapping. **No** for "single dominant community" confounding. | The walk answers a different reviewer question: whether acrocentric p-arms, sex chromosomes/PARs, all acrocentric arms, or the strongest community dominate the Mantel signal. B0 records that exclusion can strengthen the rho trajectory, which is reassuring. But A4 is already carried by community-free pointwise human correlations plus flanking unique-sequence evidence; A5 is a conservative model claim, not an exclusion-walk claim. The full 5 x 5 apparatus is too bulky for the main narrative. | Decide whether final paper includes a short sentence like "exclusion analyses gave the same or stronger matrix-level correlations" and a supplement table, or cuts this entirely from main Methods. |
| W/B bootstrap on 10,000 permutations | **Demote/keep only as provenance for Fig. 4B/community matrix effect sizes** | A4 support; not essential for A5 | **Yes** for the multi-mapping threat if flanking W/B is shown; **no** for the exact uncertainty model | W/B is a useful effect-size shorthand for within-community versus between-community contact in the block matrix, but it is community-label dependent and less clean than the pointwise community-free Spearman for A4. The bootstrap p-values can overstate precision because arm pairs and PHR pairs are not independent. Keep the W/B point estimate and, if needed, a conservative permutation provenance note; do not make bootstrap significance a headline. | Decide whether W/B p-values appear in captions/tables, or whether captions report only W/B effect sizes with pointwise Spearman as the statistical lead. |
| Mann-Whitney global test | **Cut** from active manuscript-facing apparatus unless retained only in source provenance | None uniquely; at most redundant A4 support | **Yes** for the main artifact threat | The Mann-Whitney global test compares contact distributions but does not address the hardest concern: non-independence, multi-mapping, and community/arm structure. It adds a named test without protecting a distinct abstract claim. It also competes with clearer statistics: pointwise Spearman for A4, W/B for block-matrix effect size, Mantel for matrix correlation, and flanking controls for mapping artifact. | Author must decide whether any existing figure/table still requires a Mann-Whitney p-value. If not, downstream patch should remove the named test from Methods. |
| Observed/expected inter-chromosomal normalization | **Keep, but demote to Fig. 4B/matrix-specific Methods and caption language** | A4 support for the matrix visualization; not necessary for A5 | **No**. Flanking controls mapping artifact; O/E controls marginal-contact and chromosome-size/contact-opportunity bias. | O/E is not the lead A4 statistic, but it protects the community-ordered contact matrix from the obvious objection that block enrichment is driven by chromosome-wide contact marginals or rare inter-chromosomal regimes. B0 records O/E enrichment remaining 8.6x to 34.4x. Keep it close to Fig. 4B, not as a broad statistical edifice. | Decide whether final Fig. 4B caption explicitly says "observed/expected inter-chromosomal normalization" or moves that detail entirely to Methods. |
| Per-bin-pair normalization / length-normalized contact for pointwise scatter | **Keep** | **A4 directly**; A5 indirectly as one pillar of the model | **No**. Flanking controls mapping artifact; per-bin-pair normalization controls PHR length/bin-count opportunity. | This is essential because Fig. 4A/ED1 compare PHR-pair sequence similarity to contact at exact coordinates. Without dividing by the number of bin pairs spanned by the two PHRs, larger PHRs have more possible bin-bin contacts and the pointwise Spearman could partly become a length/opportunity statistic. This is the cleanest named normalization to keep in the lead Methods. | Author decision is only wording level: keep the technical phrase in Methods, but avoid over-explaining it in Results. |
| "MAPQ0/random primary placement" policy | **Keep as a caveat, not as validation** | A4 boundary; protects against overclaiming rather than proving contact | **Partly.** Flanking control argues against broad inflation but does not prove true origin for every PHR-internal read. | The current Methods correctly says random placement is an acknowledged limitation rather than a validated method. This is essential honesty because A4 says "contact more often" but the PHR-internal data include multi-mappers. Keep the caveat compact and tie it directly to flanking unique-sequence evidence. | Author must decide whether the abstract A4 verb stays "show" or softens to "support" if the manuscript wants stricter mapping-language alignment. |
| Strict-MAPQ re-binning / MAPQ>=30 comparison | **Cut or hold pending source recovery** | Would protect A4 if complete, but currently not source-complete in B0 | **Yes enough for the current claim**, because flanking unique-sequence control is already the cleaner MAPQ-strict proxy | B0 found the script and intended output schema but not committed result tables. A source-incomplete strict-MAPQ claim should not remain as quantified manuscript apparatus. If recovered, it can be a supplement note showing that strict mapped flanks agree while PHR interiors collapse; if not recovered, cut the sentence. | Required J-task decision: locate/regenerate `comparison_summary.tsv` before retaining the claim, or authorize removal. Do not silently keep the quantified strict-MAPQ sentence. |
| "Random-ligation inflation" named O/E explanation | **Demote/rephrase** | A4 matrix support only | **No**, because it is a different bias class than mapping artifact | The phrase is too compressed. O/E can control expected inter-chromosomal contact given marginals, but "without random-ligation inflation" sounds like a stronger validated correction than the current apparatus needs. The manuscript can simply state that O/E normalization was used for matrix visualization and that the pointwise/contact conclusions do not rest on raw matrix marginals. | Author must decide whether to keep "random ligation" language. Recommendation: remove the phrase unless a methods specialist confirms the exact model. |
| Multi-resolution 5/10/20/50/100 kb W/B/Mantel stability | **Demote** to one supplement row/table | A4 robustness; not A5 | **No** for bin-size sensitivity; **yes** for mapping artifact if flanking control remains | Multi-resolution stability answers "is this a bin-size artifact?" and should be retained as compact robustness if space allows. It is not necessary in the main Results because the lead pointwise analysis is already stated at 50 kb and replicated across assays. | Decide whether a supplement table is in scope. If no supplement table exists, one Methods sentence is enough. |
| Independent Hi-C community detection / ARI on O/E matrices | **Cut from main; optional supplement only** | None uniquely; weak A4 support | **No**, but this does not target the main threat | The ARI values are modest and the analysis adds another community-detection layer when the main A4 evidence is intentionally community-free. It can reassure that contact matrices contain related structure, but it is not essential and risks distracting from the cleaner result. | Author decision: retain only if the paper wants a technical supplement on contact-derived communities. Otherwise cut. |

## Apparatus By Threat

| Threat to A4/A5 | Best required defense | Apparatus that can be demoted/cut |
|---|---|---|
| Multi-mapping or MAPQ0 read-placement artifact | MAPQ0 caveat plus flanking unique-sequence control | Strict-MAPQ comparison unless source tables are recovered; Mann-Whitney; full W/B bootstrap significance |
| PHR length or different numbers of contacted bins | Per-bin-pair / length-normalized pointwise contact | W/B bootstrap and Mantel do not replace this for Fig. 4A |
| Chromosome-size or matrix marginal-contact bias in Fig. 4B | Observed/expected inter-chromosomal normalization | Independent contact-community ARI; global Mann-Whitney |
| Single dominant community, acrocentric, or PAR/sex confounding | Compact exclusion Mantel summary if retained | Full five-set x five-resolution walk in main text |
| Resolution dependence | One compact multi-resolution stability statement/table | Repeated p-values for every resolution/statistic |
| Non-independence of thousands of PHR/arm-pair observations | Conservative language and avoiding headline astronomical p-values | Mann-Whitney/global p-values as headline evidence |

## Claim-Boundary Consequences

For A4, the essential manuscript logic should be:

1. Community-free pointwise human analyses show a positive sequence-to-contact
   relationship across Pore-C, Hi-C, and CiFi.
2. The analysis is limited by PHR multi-mapping, so PHR-internal contacts are
   not read-origin-perfect measurements.
3. Flanking unique-sequence controls argue against the result being a simple
   multi-mapping artifact.
4. Per-bin-pair normalization keeps the pointwise statistic from being a PHR
   length/bin-count artifact.
5. O/E normalization is a secondary matrix-display control for Fig. 4B.

For A5, none of the audited apparatus is individually load-bearing. A5 is
protected by triangulation: sequence communities, human 3D contact, mouse
meiosis, and pedigree evidence. The apparatus should prevent overclaiming the
human 3D line; it should not become the model claim itself.

## Author Decisions Required

These are J-task decisions for downstream fan-in or the guarded manuscript patch
task. They must not be silently applied by B5.

1. **A4 verb strength**: keep "Human Pore-C and Hi-C show..." or soften to
   "support..." given the MAPQ0/random-placement caveat. Recommendation: keep
   "show" only if the Results/Methods keep the caveat and flanking control
   adjacent; otherwise soften to "support".
2. **Strict-MAPQ sentence**: retain only if the missing result table is located
   or regenerated. Recommendation: cut the quantified sentence unless source
   recovery happens before final patching.
3. **Mann-Whitney global test**: decide whether any author still wants this
   named test in Methods. Recommendation: cut.
4. **Exclusion Mantel walk size**: decide whether to include a supplement table
   for the 5 exclusion sets x 5 resolutions or reduce it to one Methods sentence.
   Recommendation: supplement-only if a supplement exists; otherwise one
   compact sentence.
5. **W/B bootstrap p-values**: decide whether W/B remains an effect size only
   or carries permutation p-values. Recommendation: report W/B as effect size;
   avoid headline bootstrap p-values.
6. **Observed/expected language**: decide whether "random-ligation inflation"
   stays. Recommendation: replace with plain O/E marginal-normalization language
   unless a methods owner endorses the stronger phrase.
7. **Independent contact-community ARI**: decide whether a technical supplement
   wants this. Recommendation: cut from main manuscript.

## Downstream Consumption

- `manuscript-revision-bf-fanin` should use this record to collapse the B/F
  3D-contact package around pointwise human contact, flanking control, and
  conservative mouse wording from F3.
- `manuscript-revision-g0` should treat A4/A5 as supportable only at the
  bounded claim levels above, especially if F3 removes the unqualified
  zygotene-peak language.
- `manuscript-revision-paper-patch` should wait for author decisions before
  cutting or demoting any J-task apparatus.

## Validation Checklist

- Artifact exists: `paper_prep/manuscript_revision/B5_3d_apparatus_essentiality.md`.
- Every listed apparatus item is classified: exclusion Mantel walk, W/B
  bootstrap, Mann-Whitney global test, observed/expected normalization,
  per-bin-pair normalization, and named schemes including MAPQ0/random
  placement, strict-MAPQ, random-ligation language, multi-resolution stability,
  and independent contact-community detection.
- Recommendations are tied to A4/A5 abstract claims as defined above.
- J-task decisions are surfaced in "Author Decisions Required" and not silently
  applied.
