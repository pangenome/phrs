# Fig. 5 FISH / Orthogonal Validation Plan

**Status:** planning note, 2026-07-04.

This note defines a concrete validation path for the putative subtelomeric
exchange candidates used in Fig. 5 planning. It deliberately separates what is
already supported by the computational analysis from work that remains future
experimental validation. Until the future-validation steps are completed, the
autosomal events should remain described as *putative*, *candidate*, or
*PHR-exchange-compatible* rather than as cytogenetically validated exchanges.

## Scope

Primary validation target:

- `PAN027_chr9q_chr3q_PHR_candidate`: PAN027 paternal hap2 chr9q terminal
  window, compared with the PAN011 father joint target.

Secondary validation target:

- `PAN028_chr9q_chr3q_PHR_candidate`: PAN028 maternal hap1 chr9q terminal
  window, compared with PAN027, used as a transmission-consistency follow-up.

Assay controls:

- `PAR1_XY_positive_control`: known male X/Y PAR1 recombination sanity check.
- computational non-carrier or same-chromosome-only subtelomeric windows from
  the same sample class, selected before probe synthesis.

Deprioritized for immediate FISH:

- chr5q/chr1p from earlier ribbon drafts. The SweepGA/FastGA signal was not
  reproduced by wfmash, so it is not a first-pass validation target.

## Current Computational Evidence

These points are current evidence, not future FISH validation.

| Candidate | Current support | Claim boundary |
| --- | --- | --- |
| PAR1 X/Y positive control | Strict-path synteny row recovers chrX/chrY PAR1 structure in the PAN027 paternal comparison. | Positive control only; not autosomal PHR evidence. |
| PAN027 chr9q/chr3q | Strict 500 kb native chr9q window contains 445,737 bp same-chromosome chr9q context, 45,290 bp primary chr3q donor sequence, 7,273 bp chr15q and 1,207 bp chr16q side fragments, and a 493 bp low-confidence chr20q tail. The chr3q donor is concordant across SweepGA/FastGA f32/f16 and wfmash at >99.8% identity for the wfmash core. Arm-level community/Jaccard checks place chr9q/chr3q in the expected C3 high-sharing context. | Compatible with terminal PHR exchange, but not a breakpoint-resolved event call. The chr20q tail should not be interpreted as part of the event model. |
| PAN028 chr9q/chr3q follow-up | Strict PAN028 maternal chr9q window contains 449,356 bp same-chromosome chr9q context, 34,172 bp chr3q primary donor sequence, 15,166 bp chr15q side fragment, and 1,207 bp chr16q side fragment. This preserves the implicated chromosome-end classes from PAN027 into PAN028. | Transmission consistency only. The chr3q support in PAN028 is split across PAN027 haplotypes, so this is not a de novo-event proof. |
| WashU recombination comparison | Public `Cechova2025` supports pedigree assembly provenance. Personal communication from M. Cechova and T. Marschall says the relevant PAN011 to PAN027 paternal lineage does not contain an ordinary chr3 crossover that explains the signal. | Rules out a simple ordinary chr3-crossover explanation but does not prove interchromosomal exchange. Keep the negative recombination-map comparison as personal communication unless a public WashU table appears. |

Primary evidence files:

- `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/event_manifest.tsv`
- `paper_prep/_brainstorming/fig5_synteny_recombination_schematic/selected_segments.tsv`
- `paper_prep/_brainstorming/fig5_followup_transmission_check/REPORT.md`
- `paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/CROSS_ALIGNER_VALIDATION.md`
- `paper_prep/_brainstorming/fig5_donor_recipient_ribbon_draft/COMMUNITY_LINKAGE_NOTE.md`
- `paper_prep/_brainstorming/fig5_washu_recombination_comparison.md`

## Future Validation Goals

The future validation should answer three separate questions, in this order:

1. **Cytogenetic placement:** does chr3q-derived probe signal appear on the
   chr9q terminal region in the candidate carrier haplotype/cell line?
2. **Pedigree segregation:** is the same chr9q-localized chr3q signal present in
   PAN027 and transmitted to PAN028, with the expected absence or different
   configuration in the relevant parental haplotypes?
3. **Molecular structure:** do long molecules or reads connect the predicted
   chr9q recipient context to chr3q-derived sequence, and do they account for
   the chr15q/chr16q side fragments?

FISH can address goals 1 and 2 directly. Goal 3 requires orthogonal molecular
assays because FISH cannot resolve a 1-50 kb patch boundary inside a repetitive
subtelomeric 500 kb window.

## FISH Probe Design

Use oligo-paint or BAC/fosmid-derived probes only after in-silico specificity
screening against T2T-CHM13, PAN011, PAN027, PAN028, and the HPRCv2 PHR catalog.
Avoid probes dominated by TTAGGG arrays, TAR1, ubiquitous subtelomeric duplicons,
or sequence that lights up most chromosome ends.

Minimum first-pass four-color design:

| Probe | Purpose | Desired behavior |
| --- | --- | --- |
| chr9q recipient-anchor probe | Marks the candidate recipient terminal window and identifies chr9q on metaphase spreads. | Strong terminal chr9q signal in all relevant samples; minimal non-9q cross-hybridization. |
| chr3q donor-enriched probe | Tests whether chr3q-like sequence is cytogenetically present at chr9q in the candidate haplotype. | Signal on native chr3q and an additional co-localized or adjacent signal at chr9q in candidate carriers. |
| chr15q/chr16q side-fragment probe pool | Tests whether side-fragment signal is detectable and co-segregates with the candidate chr9q haplotype. | Optional weaker signal near the chr9q candidate; absence is not fatal because the side fragments are small. |
| chromosome identity controls | Confirms chromosome identity and prevents mistaking nearby subtelomeric signals for chr9q/chr3q. | Centromere or arm-paint controls for chr3, chr9, and optional chr15/chr16. |

Probe construction checklist:

- Tile candidate-specific intervals from `selected_segments.tsv`, not permissive
  multimap rows.
- Design separate candidate and control oligo pools: chr9q anchor, chr3q donor,
  chr15q side, chr16q side, PAR1 X/Y positive control, and a same-chromosome
  negative-control subtelomere.
- Mask simple repeats, segmental duplicon cores that occur on many unrelated
  ends, and any oligos with many high-identity off-target hits.
- Require an in-silico report for each oligo pool with expected target arms,
  off-target arms, number of retained oligos, total tiled bp, and predicted
  signal balance.
- Keep chr15q and chr16q side-fragment probes separate at first. Pooling them is
  acceptable only after individual probe performance is known.

## Sample Plan

Preferred samples, if viable cells or high-quality nuclei are available:

| Sample | Use |
| --- | --- |
| PAN011 | Relevant father for PAN027 paternal comparison; establishes parental chr9q/chr3q probe configuration. |
| PAN027 | Primary candidate carrier for chr9q/chr3q. |
| PAN028 | Transmission-consistency sample; expected to retain the implicated chr9q/chr3q chromosome-end classes. |
| PAN010 | Additional pedigree background control for maternal-side signals. |

Additional controls:

- at least two unrelated HPRC samples predicted to carry chr9q/chr3q C3 signal;
- at least two unrelated HPRC samples predicted not to carry the chr9q-localized
  donor configuration;
- a male positive-control sample for PAR1 X/Y co-localization;
- a technical negative control using a subtelomeric probe pair expected to stay
  on homologous/same-arm targets.

Cell material preference:

1. fresh or cryopreserved lymphoblastoid/fibroblast metaphase preparations from
   the relevant pedigree cell lines;
2. interphase nuclei from the same cell lines for co-localization frequency;
3. extracted high-molecular-weight DNA for fiber-FISH or molecular follow-up if
   metaphase quality is poor.

## Experimental Readouts

### Metaphase FISH

Primary readout:

- chr3q donor-enriched signal co-localizes with or lies immediately adjacent to
  the chr9q recipient-anchor signal on one chr9 homolog in PAN027.

Expected PAN027 result if candidate validates:

- one chr9 homolog has ordinary chr9q anchor signal only;
- the candidate chr9 homolog has chr9q anchor plus an extra chr3q donor-enriched
  signal at/near the terminal region;
- chr3q native loci retain their expected donor-probe signal.

Expected PAN028 result if transmitted:

- the chr9q-localized chr3q signal is present on the transmitted chr9q haplotype,
  with side-fragment signal considered supportive but not required.

Expected PAN011 result:

- native chr3q and chr9q signals establish parental configuration. A clean
  absence of chr3q donor signal on the relevant PAN011 chr9q homolog would
  support a lineage-specific interpretation; a present signal would shift the
  interpretation toward inherited standing subtelomeric polymorphism.

### Interphase FISH

Use interphase as a higher-throughput secondary readout, not as the primary
call. Score the fraction of nuclei where chr9q anchor and chr3q donor-enriched
signals overlap or lie within a pre-registered distance threshold. Interpret
only after metaphase confirms chromosome identity, because subtelomeric ends
physically cluster and can generate proximity without sequence transfer.

### Fiber-FISH

Use fiber-FISH if metaphase FISH shows a candidate chr9q-localized donor signal.
The desired molecule-level pattern is:

`chr9q anchor -> chr3q donor-enriched tract -> chr9q or subtelomeric continuation`

Fiber-FISH is especially useful for checking whether chr15q/chr16q side
fragments are part of the same molecule or are independent cross-hybridizing
features.

## Scoring Criteria

Pre-register these outcomes before unblinding carrier/control status:

| Outcome | Interpretation |
| --- | --- |
| chr3q donor-enriched signal appears on candidate chr9q in PAN027 and PAN028, absent from non-carrier controls | Strong cytogenetic validation of a transmitted chr9q-localized chr3q-derived subtelomeric configuration. |
| chr3q donor-enriched signal appears on PAN027 chr9q but not PAN028 | Supports PAN027-specific configuration but weakens transmission claim; re-check PAN028 haplotype assignment and probe sensitivity. |
| chr3q donor-enriched signal appears broadly on chr9q in many unrelated samples | Validates population-level chr9q/chr3q subtelomeric sharing but not pedigree-specific exchange. |
| chr3q donor-enriched signal only appears on native chr3q | FISH does not validate the candidate; computational signal may be below FISH resolution, probe design may be too specific, or the candidate may be an alignment artefact. |
| chr15q/chr16q side probes fail while chr3q validates | Acceptable first-pass result because side fragments are short and may be below robust FISH resolution. |
| chr15q/chr16q side probes light up many chromosome ends | Treat side-fragment FISH as uninformative; do not use it to support the event. |

Minimum validation threshold for a manuscript-strength FISH claim:

- at least 20 well-spread metaphases per sample for PAN011, PAN027, PAN028, and
  two non-carrier controls;
- blinded scoring by two readers or automated spot calling with manual audit;
- concordant chr9q-localized chr3q donor signal in PAN027 and PAN028;
- no comparable chr9q-localized chr3q signal in non-carrier controls;
- image examples plus a per-cell scoring table deposited with the figure
  source files.

## Orthogonal Non-FISH Validation

FISH should be paired with at least one sequence-level or molecule-level assay
before upgrading the language from "candidate" to "validated exchange".

Recommended assays:

| Assay | What it tests | Acceptance criterion |
| --- | --- | --- |
| Targeted adaptive-sampling ONT or Cas9-enriched HiFi over chr9q/chr3q intervals | Whether reads traverse chr9q context into chr3q-derived sequence in PAN027/PAN028. | Multiple independent molecules span the predicted junction or patch boundary with haplotype-consistent phase. |
| Fiber-FISH or optical genome mapping | Whether the donor and recipient probe order exists on single molecules. | Molecule-level co-linearity matching the predicted chr9q/donor/side-fragment layout. |
| ddPCR/qPCR copy-number assay for chr3q donor tract on chr9q carriers | Whether candidate carriers have the expected donor-tract dosage. | Carrier/control dosage separates cleanly after normalizing native chr3q copy. |
| Single-gamete or sperm long-read assay from a carrier father, if available | Whether exchange-compatible products appear in gametes. | Direct recombinant molecule/haplotype classes are observed; absence is not dispositive unless power is high. |
| Re-analysis when WashU recombination annotations become public | Whether official pedigree recombination calls expose or refute the candidate. | Replace personal communication with table/figure citation and update Fig. 5 text. |

## Manuscript Language Until Validation

Use:

- "putative chr9q/chr3q PHR-exchange-compatible candidate";
- "concordant across independent alignment strategies";
- "consistent with transmission into PAN028";
- "not explained by an ordinary chr3 crossover in the current WashU comparison";
- "requires FISH or molecule-level validation."

Avoid:

- "FISH-validated";
- "confirmed interchromosomal crossover";
- "breakpoint-resolved exchange";
- "de novo recombination event";
- "validated chr15q/chr16q side-fragment insertion."

## Immediate Next Actions

1. Extract candidate-specific intervals from `selected_segments.tsv` into a probe
   design manifest with columns for event, arm, haplotype, native interval,
   role, bp, and source file.
2. Run oligo uniqueness screening against CHM13 plus PAN011/PAN027/PAN028
   assemblies and HPRCv2 PHR intervals.
3. Select two non-carrier and two computational carrier controls from HPRCv2
   using the same chr9q/chr3q support criteria.
4. Confirm whether PAN011, PAN027, PAN028, and PAN010 metaphase-quality material
   or nuclei are available.
5. Perform PAR1 positive-control hybridization before interpreting autosomal
   candidate probes.
6. Score metaphase FISH first; use interphase and fiber-FISH as supporting
   assays after chromosome identity is established.

## Dependency Notes

This plan inherits the claim boundary from
`paper_prep/_brainstorming/fig5_washu_recombination_comparison.md`: cite
`Cechova2025` for assembly provenance, keep the chr3-negative recombination-map
comparison as personal communication unless the WashU recombination table is
public, and do not represent the chr9q/chr3q candidate as event-level validation
without future cytogenetic or molecule-level evidence.
