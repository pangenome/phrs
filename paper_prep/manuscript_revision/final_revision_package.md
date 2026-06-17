# Final Manuscript Revision Package

Date: 2026-06-17

Task: `manuscript-revision-final-fanin`

This package is the handoff to `manuscript-revision-paper-patch`. It synthesizes
completed upstream records and records decisions that remain author-only. It does
not apply new manuscript edits and does not invent decisions absent from the
prior decision records.

## Source Records Read

- `paper_prep/manuscript_revision/A_mechanical_fixes.md`
- `paper_prep/manuscript_revision/BF_3d_contact_synthesis.md`
- `paper_prep/manuscript_revision/CD_continuum_community_synthesis.md`
- `paper_prep/manuscript_revision/E_pedigree_audit.md`
- `paper_prep/manuscript_revision/G_abstract_intro_title_package.md`
- Supporting records: `00_inventory.md`, `01_fanout_graph.md`,
  `02_operating_rules.md`, `C0c_D1_resolution_sampling.md`,
  `F1_F2_orphan_audit.md`, and
  `F3_mouse_shape/F3_mouse_shape_report.md`.

## Completed and Blocked Task Table

| Area | Upstream task/artifact | Status | Manuscript authority/result |
| --- | --- | --- | --- |
| A mechanical reconciliation | `manuscript-revision-a-fanin`; `A_mechanical_fixes.md` | Completed | Applied supported mechanical edits to `submission/paper.tex`, `submission/bibliography.bib`, `submission/compile.sh`, `submission/README.md`, and `submission/BUILD_LOG.md`. Build attempted but blocked by missing local TeX package `geometry.sty`. |
| B/F 3D contact, Mantel, mouse and orphan evidence | `manuscript-revision-bf-fanin`; `BF_3d_contact_synthesis.md` | Completed | Decision-support only. Recommends conservative 3D evidence hierarchy, flanking unique-sequence control framing, MAPQ0 wording, mouse-softening, and FST/cM-Mb handling. |
| C/D continuum, community and sampling synthesis | `manuscript-revision-cd-fanin`; `CD_continuum_community_synthesis.md` | Completed | Decision-support only. Recommends two-tier continuum language, UPGMA as ordering not phylogeny, removal of NJ/bootstrap load-bearing language, and honest resolution/sampling wording. |
| E pedigree circularity/claim audit | `manuscript-revision-e1`; `E_pedigree_audit.md` | Completed | Decision-support only. Recommends treating Leiden as a credibility constraint, not the definition of the tested set; right-sizing replication and artifact claims. |
| G abstract/introduction/title/denominator package | `manuscript-revision-g0`; `G_abstract_intro_title_package.md` | Completed | Decision package only. Provides accepted-decision placeholders and text variants; no direct manuscript edit. |
| Final package | `manuscript-revision-final-fanin`; this file | Completed by this task | Handoff package for guarded manuscript patch. |
| Local LaTeX compile | `submission/BUILD_LOG.md`; A fan-in validation | Blocked by environment | `bash compile.sh` fails before manuscript content is processed because `geometry.sty` is not installed locally. This is not a manuscript-content failure. |
| Formal Leiden sample-resampling stability | `C0c_D1_resolution_sampling.md`; `CD_continuum_community_synthesis.md` | Not completed; follow-up only | Slurm-ready wrapper exists, but no formal PHR/sample-resampling ARI/NMI/co-clustering result exists. Do not claim formal sampling stability. |
| Author-only decisions | See decision register below | Blocked pending named author choice | Downstream patch may apply only recommendations already selected by an author/J decision record or clearly marked as non-judgment mechanical cleanup. |

## Decisions Still Requiring the Author

These are the unresolved choices explicitly surfaced by prior records. The patch
agent should not silently resolve them.

| Decision | Options from prior records | Recommended default in records | Downstream rule |
| --- | --- | --- | --- |
| Abstract q-arm list (`G-0a`/`G-0b`) | Keep q-arm list as a named local-density example, or remove the list if C0 is judged insufficient/artifact-prone. | `G_abstract_intro_title_package.md` recommends keeping the list only if authors accept the C0 dense-neighborhood framing. | If no author decision exists, use the safer non-list abstract variant or leave a placeholder for author review. |
| Title strength (`G-1`/`G-2`) | Keep "unorthodox recombination" title as an opportunity/compatibility claim, or soften to "sequence exchange". Avoid FST-promoted or over-discrete titles. | Keep the current title if authors accept the evidence boundary; otherwise use the softer title. | Do not retitle without author choice. |
| Denominator presentation (`G-3`) | Reconcile `465 assemblies`, `464 HPRC haplotypes from 232 individuals + CHM13`, `18,827 flanks`, `15,668 PHRs`, `41/48 arms`. | Use the explicit denominator language from G and A records. | Mechanical denominator cleanup is allowed where it is purely factual; do not change biological interpretation. |
| 3D evidence hierarchy | Lead with flanking unique-sequence control and treat PHR-internal all-points as downstream aggregate evidence, or retain the current stronger ordering. | BF recommends flanks as primary anti-artifact support and PHR-internal scatter as descriptive. | Apply only if accepted for the patch; otherwise preserve current claim while adding caveats already supported by A/B records. |
| MAPQ0/random primary wording | Explain disabled MAPQ filters as necessary to retain paralogous signal with one random primary alignment retained; avoid implying read-origin certainty. | BF recommends this wording. | This is safe as a caveat, but the evidence hierarchy change remains author-sensitive. |
| Mouse meiosis claim size | Keep "peaks at zygotene" only as descriptive/per-PHR contrast; avoid broad "prophase-I association" or overclaiming zygotene specificity. | BF/F3 recommend softening: positive through prophase I, with strongest descriptive support at zygotene and limited arm-block uncertainty. | If Fig. 4C remains, revise caption/results to match F3/BF. |
| FST role | Retain only as background-level ancestry structure/control, move to Extended Data/Methods, or cut. | F1/BF/G recommend not promoting FST to title or central result. | Do not introduce a main FST result unless the matched-background caveat is also present. |
| cM/Mb line | Clarify as this manuscript's derived analysis using Lalli 2025 recombination-rate input, or cut. | F1/BF recommend cutting if word count is tight. | Do not write "Lalli reported the anti-correlation." |
| Pedigree claim size | Present WashU evidence as consistent with recent exchange preferentially within sequence communities, with family/meiotic-unit limits and artifact caveats. | E recommends this softer framing. | Do not claim independent 538-meiosis replication or a validated event catalog. |
| Author Contributions and Competing Interests | Fill with author-confirmed statements. | A notes these remain TODO-only/author-truth statements. | Must be supplied or confirmed by authors; agents should not invent them. |
| Source-data accessions | Add exact public accessions/source notes for HG002 Pore-C/Hi-C/CiFi, Dip-C, sperm scHi-C, and mouse meiotic Hi-C where available. | BF/A3 recommend source-data rows and public citations rather than internal paths. | Do not invent accession numbers. If unavailable, leave an explicit source-data TODO. |

## Files Changed So Far

Already changed by upstream A fan-in:

| File | Change class | Record |
| --- | --- | --- |
| `submission/paper.tex` | Mechanical count/citation cleanup; removal of internal report/script/path pointers; bibliography placement cleanup; multibib/cite cleanup; construction-artifact cleanup. | `A_mechanical_fixes.md` |
| `submission/bibliography.bib` | Corrected `acrocentric_Guarracino2025ape` metadata while preserving the key. | `A_mechanical_fixes.md` |
| `submission/compile.sh` | Cleans stale `Meth.*` and `Supp.*` auxiliaries after moving to one reference list. | `A_mechanical_fixes.md` |
| `submission/README.md` | Updated build/story instructions for the current five-main-figure manuscript and bibliography state. | `A_mechanical_fixes.md` |
| `submission/BUILD_LOG.md` | Recorded the 2026-06-17 compile attempt and `geometry.sty` blocker. | `A_mechanical_fixes.md` |

Analysis/decision artifacts added upstream under `paper_prep/manuscript_revision/`:

| Artifact group | Key files |
| --- | --- |
| Mechanical audits | `A1_sample_counts.md`, `A2_bibliography_audit.md`, `A3_artifacts_audit.md`, `A4_guarracino_doi.md`, `A5_A7_mechanical_audit.md`, `A_mechanical_fixes.md` |
| 3D/contact/FST/mouse records | `B0_3d_inventory.md`, `B1_B3_3d_decision_record.md`, `B4_pvalue_mantel_audit.md`, `B5_3d_apparatus_essentiality.md`, `BF_3d_contact_synthesis.md`, `F1_F2_orphan_audit.md`, `F3_mouse_shape/*` |
| Continuum/community records | `C0_continuum/*`, `C0c_D1_resolution_sampling.md`, `C1_tree_essentiality.md`, `C2_bootstrap_audit.md`, `C3_qarm_language.md`, `CD_continuum_community_synthesis.md` |
| Pedigree/title packages | `E_pedigree_audit.md`, `G_abstract_intro_title_package.md` |

This final task adds only `paper_prep/manuscript_revision/final_revision_package.md`.

## Proposed Figure and Table Changes

| Figure/table | Proposed change | Source decision record | Author sensitivity |
| --- | --- | --- | --- |
| Fig. 2B/C heatmaps | Relabel left panel as UPGMA average-linkage ordering, not "UPGMA tree" or "phylogeny"; describe both panels as locally dense blocks on a broader similarity continuum. | `CD_continuum_community_synthesis.md` | Low if only relabeling/caption caveat; higher if abstract q-arm list is changed. |
| Fig. 2 q-arm language | Use "q-arm neighborhood" or "local-density example"; avoid "clade", "closed class", and "q-arm sextet" as a biological noun. | `CD_continuum_community_synthesis.md`; `G_abstract_intro_title_package.md` | Requires author decision if keeping/removing named q-arm list in abstract. |
| NJ/character-bootstrap material | Remove or de-emphasize from active claims; if retained, keep only as exploratory comparison/caveat. | `CD_continuum_community_synthesis.md` | Low for removing load-bearing text; author choice if a short caveat is desired. |
| Fig. 4A | Caption/results should label PHR-internal pointwise scatter as aggregate/descriptive and not as read-origin proof at every paralogous locus. | `BF_3d_contact_synthesis.md` | Medium because it changes evidentiary hierarchy. |
| Fig. 4B | Keep HG002 Pore-C community matrix as main visual support; report B/W as effect size and avoid headlining global Mann-Whitney p-values. | `BF_3d_contact_synthesis.md` | Medium. |
| Fig. 4C mouse | Soften from a strong zygotene-specific peak to a descriptive stage series: per-PHR-pair contrast supports zygotene over leptotene, but zygotene is not clearly above pachytene/diplotene under all uncertainty summaries. | `F3_mouse_shape_report.md`; `BF_3d_contact_synthesis.md` | Medium. |
| Fig. 5 pedigree | Caption should say within-Leiden/community membership is a credibility constraint and null-test label, not the rule that defines the tested event set. Add family/unit and artifact caveats. | `E_pedigree_audit.md` | Medium/high because it changes claim size. |
| Source-data table | Add processed repo paths plus public raw-read dataset/citation/accession rows for human 3D and mouse datasets. | `BF_3d_contact_synthesis.md`; `A3_artifacts_audit.md` | Requires real accession/source lookup or author input; do not fabricate. |
| FST/cM-Mb table/legend | If retained, move to Methods/Extended Data or annotate as control/background. Clarify cM/Mb provenance as our analysis using Lalli 2025 input. | `F1_F2_orphan_audit.md`; `BF_3d_contact_synthesis.md` | Requires author decision if adding/removing panels. |

## Exact Manuscript Edit Plan

The downstream patch should edit `submission/paper.tex` only where the plan below
is either mechanical or backed by an accepted decision record. Current anchors
refer to the manuscript inspected on 2026-06-17 after A fan-in.

1. **Title and abstract (`paper.tex:48-49`)**
   - Replace the abstract with one of the G-package variants.
   - If the author accepts q-arm dense-neighborhood support, use G Variant A or
     B and keep the named q-arm list as an illustrative local-density example.
   - If not, use G Variant C and remove the q-arm list.
   - Keep the denominator language: 465 near-complete assemblies; 464 HPRC v2
     haplotypes from 232 individuals plus CHM13; 18,827 flanks; 15,668 PHRs;
     41 of 48 arms.

2. **Community Results (`paper.tex:156-185`)**
   - Replace "discrete blocks" with "locally dense blocks embedded in a broader
     similarity continuum".
   - Replace "UPGMA tree/phylogeny" language with "UPGMA average-linkage
     ordering".
   - Keep "15 arm-level sequence-similarity communities" only as Leiden
     operating-point language, not as closed biological classes.
   - Use CD wording for named recovered systems: PAR1, PAR2, acrocentric p-arms,
     the 10p/18p TUBB8B system, the 4q/10q DUX4 system, and conditional q-arm
     neighborhood language.

3. **Mouse/3D Results (`paper.tex:192-225`)**
   - Reorder or rewrite the paragraph so human flanking unique-sequence controls
     are the primary anti-artifact support and PHR-internal Fig. 4A is
     descriptive aggregate evidence.
   - Replace the strong mouse "peaks at zygotene" claim with BF/F3-compatible
     wording: mouse meiotic Hi-C supports a positive sequence-contact
     relationship during prophase I, with the clearest descriptive zygotene
     contrast over leptotene but not a uniformly significant zygotene-specific
     maximum under all uncertainty summaries.
   - Avoid global pointwise p-value emphasis; use Mantel/Spearman effect sizes
     and intervals where available.

4. **Pedigree Results (`paper.tex:233-258`)**
   - Replace any implication that community membership defines the tested event
     set. State that high-quality inter-chromosomal patches are identified first
     by untangle/geometry/parent-child consistency and that within-community
     status is an orthogonal credibility label/null-test category.
   - Right-size replication: one three-generation WashU T2T pedigree with nested
     patches/events, not hundreds of independent biological replicates.
   - Add artifact-boundary wording from E: assembly-fragmentation artifacts are
     bounded for CEPH1463 by expected within-community depletion, but the
     WashU per-event artifact fraction is not directly estimated.

5. **Figure legends (`paper.tex:324-333`, `373-398`, `406-414`)**
   - Fig. 2: apply CD legend replacement, especially "average-linkage UPGMA
     ordering" and "locally dense blocks on a broader similarity continuum".
   - Fig. 4A/B/C: apply BF captions; explicitly describe PHR-internal contact as
     aggregate and add MAPQ0/random-primary caveat if the caption has room.
   - Fig. 5: apply E caption block or equivalent softer wording.

6. **Methods: community detection and clustering (`paper.tex:491-519`)**
   - Keep Leiden method details but add the resolution-scan qualification from
     CD/C0c: arm-level 15-community state is stable across resolution 1.13-1.18
     with selected 1.16 by maximum silhouette; sequence-level 50-community
     partition is a constrained k-NN Leiden setting (`k=75`, resolution `0.8`),
     not the unconstrained modularity maximum.
   - Replace UPGMA/NJ language so UPGMA is an ordering/algorithm comparison and
     NJ/character bootstrap is not load-bearing. Delete the NJ/bootstrap
     subsection if it remains claim-bearing.
   - State that formal sample-resampling stability is not claimed.

7. **Methods: Hi-C/Pore-C/CiFi (`paper.tex:521-591`)**
   - Preserve the MAPQ0 explanation from A/BF but make it more direct:
     paralogous PHR tracts cannot be uniquely assigned by MAPQ filtering, so raw
     reads were realigned with MAPQ filters disabled and one random primary
     alignment retained; flanking unique-sequence controls provide the reported
     artifact-control statistic.
   - Remove any statement that a strict-MAPQ PHR-internal result proves the same
     biological effect unless a committed result table is named.
   - Add public source-data/citation rows where exact accessions are known.

8. **Methods: mouse (`paper.tex:595-607`)**
   - Replace "peak at the zygotene bouquet" as a statistical claim with F3's
     exact stage-series framing.
   - Keep the Zuo 2021 and B6/CAST T2T provenance; add source-data accession
     rows only when exact accessions are verified.

9. **Methods: pedigree (`paper.tex:612-624`)**
   - Use E's Methods block: candidate patches from `odgi untangle nth-best=1`,
     high-quality filters, geometry/parent-child classification, and
     within-Leiden/community status as credibility constraint/null-test category.
   - Preserve PRDM9-independent recombination citations but avoid overclaiming
     mechanism from one pedigree.

10. **Limitations and data availability (`paper.tex:636-657`)**
    - Keep uncertainty statements but remove orphan FST/cM-Mb wording unless the
      author elects to retain those as controlled background analyses.
    - Where FST remains, use F1 wording: AFR/non-AFR values are 0.10-0.15 but
      statistically indistinguishable from matched genome-wide autosomal
      background.
    - Replace internal `/moosefs` or internal-report provenance with public
      citations/source-data rows. Do not invent accession numbers.

11. **Declarations (`paper.tex:282-286`)**
    - Leave Author Contributions and Competing Interests for author-confirmed
      text. Do not synthesize these statements.

12. **Post-patch validation**
    - Re-run targeted searches from A fan-in: stale counts, `citeMeth`/multibib,
      private paths, internal report/script pointers, `XTR`, k-symbol collision
      patterns, "discrete blocks", "phylogeny", "clade", "q-arm sextet", and
      unsupported "peak at zygotene" wording.
    - Attempt `cd submission && bash compile.sh`. If `geometry.sty` remains
      missing, record the same environment blocker in `BUILD_LOG.md`.

## Compute Artifacts and Slurm Logs

| Area | Durable artifacts | Slurm/log status |
| --- | --- | --- |
| C0a arm-level continuum | `C0_continuum/compute_arm_level_continuum.R`; `arm_pair_similarity_long.tsv`; `community_similarity_summary.tsv`; `c6_*`; `named_system_peak_similarities.tsv`; `arm_level_similarity_diagnostic.{pdf,png}` | Head-node/lightweight record; no Slurm log required by the fan-in records. |
| C0b sequence-level continuum | `C0_continuum/scripts/analyze_sequence_continuum.R`; `results/sequence_similarity_distribution.tsv`; `results/sequence_similarity_peaks.tsv`; `results/high_similarity_summary.tsv`; `results/c6_neighborhood_density.tsv`; `results/run_metadata.tsv`; `results/input_file_inventory.tsv`; `C0b_sequence_level_report.md` | Slurm job `1703691` on `octopus02` was cancelled after scanning rows through 15000/15667. Slurm job `1703692` on `octopus02` completed at `2026-06-17T15:40:17Z` with 4 CPUs and 80G request. Logs: `logs/c0b_seq_continuum.1703691.{out,err}` and `logs/c0b_seq_continuum.1703692.{out,err}`. |
| C0c/D1 resolution scans | `scan_file_inventory.tsv`; `arm_leiden_best_resolution.tsv`; `arm_leiden_resolution_by_k.tsv`; `sequence_leiden_scan_summary.tsv`; `sequence_upgma_scan_summary.tsv`; `summarize_resolution_scans.py`; `run_leiden_sampling_stability.slurm.sh` | No sampling-stability Slurm job was completed. Wrapper is Slurm-ready but intentionally exits unless a concrete `run_leiden_sampling_stability.R` implementation is supplied. |
| F3 mouse shape | `F3_mouse_shape_analysis.py`; `F3_stage_series.tsv`; `F3_zygotene_contrasts.tsv`; `F3_input_inventory.tsv`; `F3_mouse_shape_report.md` | No Slurm job submitted; F3 report records "No Slurm job was submitted. There is no Slurm job ID." |
| F1 FST/cM-Mb audit | `F1_F2_orphan_audit.md`; repo-local `scripts/ci/fst_block_jackknife.tsv` and related outputs cited there | Local/repo audit; no new Slurm record in final fan-in inputs. |
| A mechanical build | `submission/BUILD_LOG.md` | Build attempted locally; blocked by missing `geometry.sty`. |

## Short Risk List

- **Author-decision risk:** high-impact title, abstract q-arm list, 3D evidence
  hierarchy, mouse-claim strength, FST/cM-Mb retention, and declaration sections
  still need author confirmation.
- **Compile-environment risk:** the current local TeX install lacks
  `geometry.sty`, so manuscript compile validation is blocked until the TeX
  environment is fixed.
- **Provenance/accession risk:** several 3D datasets have internal processed TSV
  paths and citations but lack exact public accession/source-data rows in the
  inspected artifacts.
- **Overclaim risk:** old wording around "discrete blocks", "phylogeny",
  "q-arm sextet", "peak at zygotene", and independent pedigree replication
  would overstate the support if not patched.
- **Sampling-stability risk:** resolution scans and UPGMA comparisons exist, but
  formal sample-resampling stability for Leiden assignments is not a completed
  result.
- **Figure-source drift risk:** active `submission/fig` assets and older
  `paper_prep/figures`/deck-derived placeholders may not all reflect the final
  softened claims unless captions and source-data manifests are updated
  together.

## Validation

- Required deliverable exists: `paper_prep/manuscript_revision/final_revision_package.md`.
- Required contents covered: completed/blocked task table; author-decision
  register; files changed so far; proposed figure/table changes; exact
  manuscript edit plan; compute artifacts and Slurm logs; short risk list.
- Decision discipline: all author-only choices are carried forward from prior
  records and left unresolved rather than silently selected here.
