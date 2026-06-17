# B/F 3D Contact Evidence Synthesis

Date: 2026-06-17
Task: `manuscript-revision-bf-fanin`

## Purpose

This document synthesizes the completed B/F revision inputs for 3D contact,
Mantel/p-value handling, mouse meiosis, and orphan F_ST/cM-Mb evidence. It is a
decision-support artifact for downstream manuscript patching, not an edit to
`submission/paper.tex`.

The synthesis separates:

- data-supported results that can be stated now;
- proposed manuscript framing that still needs author/J-task approval;
- provenance and source-data notes that must be converted from internal paths to
public citations/accessions before submission.

## Sources Read

Primary dependency artifacts:

- `paper_prep/manuscript_revision/B0_3d_inventory.md`
- `paper_prep/manuscript_revision/B1_B3_3d_decision_record.md`
- `paper_prep/manuscript_revision/B4_pvalue_mantel_audit.md`
- `paper_prep/manuscript_revision/B5_3d_apparatus_essentiality.md`
- `paper_prep/manuscript_revision/F1_F2_orphan_audit.md`
- `paper_prep/manuscript_revision/F3_mouse_shape/F3_mouse_shape_report.md`
- `paper_prep/manuscript_revision/F3_mouse_shape/F3_stage_series.tsv`
- `paper_prep/manuscript_revision/F3_mouse_shape/F3_zygotene_contrasts.tsv`
- `paper_prep/manuscript_revision/A3_artifacts_audit.md`

Active manuscript and source manifests:

- `submission/paper.tex`
- `paper_prep/figures/fig3/sources.tsv`
- `paper_prep/figures/fig4/sources.tsv`
- `paper_prep/figures/ed8/sources.tsv`
- `end-to-end-report/report/05_hic_validation.md`
- `end-to-end-report/report/06_dipc_validation.md`
- `end-to-end-report/report/08_mouse.md`
- `end-to-end-report/report/07_integrated.md`

## Executive Synthesis

The strongest B/F revision path is:

1. Keep human PHR-internal pointwise contact as the primary human 3D result,
   because it is within-genome, per-sequence-pair, community-free, and directly
   tests contact at PHR coordinates.
2. Make flanking unique-sequence analysis the primary anti-mapping-artifact
   control, placed adjacent to the first human 3D claim.
3. Treat pointwise Spearman rho and sample size as descriptive effect sizes; do
   not use nominal pointwise p-values as inference.
4. Use row+column Mantel/permutation statistics, with finite-permutation floor
   language, as the matrix-level inferential family when an inferential p-value
   is needed.
5. Demote the broad 3D-control apparatus to Methods/Supplement unless a final
   figure needs it: W/B bootstrap p-values, Mann-Whitney global tests,
   exclusion walks, multi-resolution grids, independent contact-community ARI,
   Dip-C/sperm forest-plot p-value rows, and RPE-1 breadth.
6. Reframe mouse meiosis as broad prophase-I support for the
   sequence-to-contact relationship, not as an inferentially resolved zygotene
   peak.
7. Keep F_ST only as demoted population-structure context if retained at all;
   keep cM/Mb only as an honest null/limitation with clear Lalli-2025 input
   provenance.

## Evidence Hierarchy

| Evidence family | Role for active manuscript | Data-supported statement | Framing boundary |
|---|---|---|---|
| Human PHR-internal pointwise contact, Fig. 4A/ED1 | Primary human 3D evidence | Contact increases with PHR sequence similarity in HG002 Pore-C and replicates directionally in CHM13 Hi-C, HG002 Hi-C, and HG002 CiFi. | Descriptive rho/n only; contacts are aggregate PHR-window measurements, not read-origin truth at each paralogous locus. |
| Human flanking unique-sequence control | Primary anti-artifact control | Adjacent centromere-ward flanks show weaker but positive sequence/contact concordance in key datasets and positive arm-level Mantel in most datasets. | Flanks support broader subtelomeric-domain proximity and defend against broad MAPQ0 artifact; they do not prove exact PHR-internal read origins. |
| Row+column Mantel/permutation | Primary inferential statistic when needed | Arm-level similarity/contact matrices are positively concordant across most human datasets and mouse stages. | Use finite-permutation floor language; avoid nominal pointwise p-values and p = 0.0. |
| HG002 Pore-C community matrix, Fig. 4B | Main visual support | Community-ordered contact matrix shows within-community block structure; O/E normalization addresses marginal/chromosome-size bias. | Keep B/W as effect size; do not headline astronomical W/B/Mann-Whitney p-values. |
| Mouse meiotic Hi-C | Supportive cross-species and prophase-I evidence | Sequence/contact association is positive across leptotene, zygotene, pachytene, and diplotene; arm-collapsed point estimate is largest at zygotene. | Direct clustered contrasts do not resolve a zygotene-specific peak. |
| Dip-C GM12878 and sperm scHi-C | Supplementary/corroborative 3D support | Single-cell 3D distances support within-community proximity, with GM12878 positive and sperm community-based enrichment strong. | Keep caveats: PBMC negative control, sperm arm-level community-free aggregate weak/wrong-direction, MAPQ0/random placement retained. |
| RPE-1 Pore-C/CiFi | Optional supplement / likely cut from active main story | RPE-1 datasets support cell-type breadth and self-discovered community validation. | Reviewer-era breadth; not needed for the active main claim unless a supplement wants it. |
| F_ST | Demoted context only | Cross-arm/self-arm structural haplotypes preserve ordinary ancestry structure at approximately genome-wide background levels. | Not evidence for subtelomere-specific differentiation or ongoing exchange. |
| cM/Mb anti-correlation | Honest null/limitation only | Apparent anti-correlation using Lalli 2025 T2T-CHM13 short-read recombination rates collapses after low-callability arms are excluded. | This is our reanalysis using Lalli's map, not a Lalli-reported PHR result. |

## Primary Data-Supported Results

### Human Pointwise Contact

Supported:

- HG002 Pore-C pointwise Spearman: rho = 0.381, n = 2,830 for PHR
  Jaccard versus length-normalized contact at PHR coordinates. B0 records this
  as the lead Fig. 4A source from `data/human_HG002_porec_50000bp_seqlevel.tsv`
  and upstream `/moosefs/.../analysis/human/community_free/human_HG002_porec_50000bp_seqlevel.tsv`.
- Replicate pointwise effect sizes: CHM13 Hi-C rho = 0.716, n = 652;
  HG002 Hi-C rho = 0.662, n = 2,544; HG002 CiFi rho = 0.191, n = 2,757.
- The measured contact is length-normalized by the number of 50 kbp bin pairs
  spanned by the two PHRs, so longer PHRs are not given more contact
  opportunities by construction.

Do not overclaim:

- Do not say pointwise p-values prove independent read-level contact. The same
  PHRs, arms, haplotypes, contact-map rows/columns, and random-placement
  histories recur across dots.
- Do not say MAPQ0/random primary placement identifies the true chromosome of
  origin for every paralogous read.

Exact recommended reporting class:

> Report pointwise rho and n as descriptive effect sizes. Use "descriptive
> pointwise Spearman" once in the Results/caption/Methods, then use rho/n
> consistently. Remove nominal pointwise p-values from Fig. 4A, Fig. 4C, ED1,
> and the matching figure-label generators.

### Flanking Unique-Sequence Control

Supported:

- Human flanking community-free Spearman is weaker than the PHR-internal signal
  but positive in key datasets: CHM13 Hi-C rho = 0.136, HG002 Hi-C rho = 0.131,
  HG002 Pore-C rho = 0.038. HG002 CiFi is weak/non-significant, and HG02148 /
  NA19036 can become NaN because assembly fragmentation concentrates nonzero
  flanking Jaccard pairs on fragmented chromosomes.
- Human 100 kb flanking Mantel remains positive in most datasets: CHM13
  rho = 0.522, HG002 Hi-C rho = 0.520, HG002 Pore-C rho = 0.314, with HG02148
  weak/non-significant.
- The end-to-end report explicitly treats flanking regions as the clean-mapping
  control for multi-mapping artifact.

Do not overclaim:

- Flanking evidence is not the same biological measurement as contact inside
  homologous PHR tracts. It supports broader subtelomeric-domain proximity and
  argues against a broad read-placement artifact.

Exact recommended framing:

> Because paralogous PHR tracts cannot be uniquely assigned by MAPQ filtering,
> we treated PHR-window contacts as aggregate coordinate-level measurements and
> tested adjacent centromere-ward flanks as a clean-mapping control. The
> flanking analyses show the same positive direction in key datasets, arguing
> against broad inflation from retaining multi-mappers, while not establishing
> read-level origin inside the PHRs.

### Mantel and P-Value Handling

Supported:

- Arm-level row+column Mantel/permutation tests are the correct inference class
  for matrix-level sequence/contact concordance when a p-value is needed.
- Human 50 kb PHR Mantel values in the end-to-end report are positive for 7 of
  8 datasets: CHM13 0.656, HG002 Hi-C 0.657, HG02559 0.397, HG00658 0.276,
  NA19036 0.266, HG002 Pore-C 0.486, HG002 CiFi 0.308; HG02148 is weak
  (rho = 0.152, p = 0.085).
- Mouse arm-level Mantel is positive across stages, but should use
  finite-permutation floor wording rather than astronomical pointwise p-values.

Exact recommended language:

> Pointwise Spearman correlations are reported as descriptive effect sizes
> because sequence-pair observations share PHRs, chromosome arms, contact-map
> rows/columns and read-placement histories. Matrix-level inference used
> row+column Mantel/permutation tests. When no permutation exceeded the observed
> statistic in 10,000 permutations, report the finite floor as p < 1/10,001
> rather than p = 0 or an astronomical nominal pointwise p-value.

### MAPQ0 / Random Primary Placement

Supported by end-to-end documentation:

- `end-to-end-report/report/05_hic_validation.md` states that Hi-C uses
  HiC-Pro/Bowtie2 with `MIN_MAPQ = 0`, `RM_MULTI = 0`; Pore-C and CiFi use
  minimap2/pairtools with one primary alignment and no MAPQ filter.
- The same section states that each multi-mapped read keeps exactly one
  randomly chosen alignment, which adds symmetric noise to contact matrices and
  makes individual pair-level contacts in repetitive regions unreliable.
- `end-to-end-report/report/06_dipc_validation.md` states that Dip-C and sperm
  scHi-C use MAPQ=0 (`sam2seg -q 0`, `hickit --min-mapq=0`) with one primary
  alignment per read, retaining subtelomeric reads that default filters remove.

Exact recommended language:

> Hi-C, Pore-C, CiFi, Dip-C and sperm scHi-C were reprocessed with multi-mappers
> retained because default MAPQ-filtered deposited maps remove the paralogous
> subtelomeric signal-bearing sequence. Under this policy, each multi-mapped
> molecule contributes one primary placement rather than all possible
> placements. PHR-internal contacts are therefore aggregate coordinate-level
> measurements and not read-level truth sets. The adjacent flanking
> unique-sequence analyses provide the principal control against broad
> MAPQ0-driven inflation.

Exact language to avoid:

- "MAPQ0 validates true PHR contacts."
- "Random placement proves the contact is real for each PHR pair."
- "Strict MAPQ confirms the PHR-internal result" unless a source table is
  located or regenerated.

### Strict-MAPQ Status

Supported:

- `scripts/hic/mapq_strict_d_peerq1.py` exists and documents a planned
  workflow from upstream `.allValidPairs` to strict-MAPQ pairs, `.cool`,
  `.mcool`, `global_test.tsv`, `comparison_row.tsv`, and `comparison_summary.tsv`.
- B0/B1/B5 did not find committed strict-MAPQ result tables sufficient to keep
  a quantified strict-MAPQ claim as a pillar.
- A3 already provides safer public-facing wording if the strict-MAPQ table is
  absent.

Recommended status:

> Treat strict-MAPQ as source-incomplete. If final source-data tables are not
> available, replace any quantified strict-MAPQ claim with: "A strict-MAPQ
> re-binning workflow for upstream valid-pair files (both mates MAPQ >= 30) is
> included in the public analysis-code repository; the unique-sequence flanking
> analysis provides the reported artifact-control statistic."

### Mouse Meiosis

Supported:

- F3 recomputed cached 1 Mb / 50 kb mouse stage series:
  sequence-level exact-PHR rho = 0.372, 0.425, 0.428, 0.416 for leptotene,
  zygotene, pachytene, diplotene.
- Arm-pair-collapsed rho = 0.680, 0.715, 0.677, 0.574 for the same stages,
  with the largest point estimate at zygotene.
- Direct clustered contrasts do not establish a resolved zygotene-specific
  peak. Sequence-level zygotene is higher than leptotene only
  (delta rho = +0.062, 95% bootstrap CI 0.004 to 0.120, p = 0.033), but not
  higher than pachytene or diplotene. All arm-pair-collapsed zygotene contrasts
  have bootstrap intervals crossing zero.

Do not overclaim:

- Do not write that mouse sequence/contact correlation "peaks at zygotene" as
  an inferential result.
- Do not use mouse pointwise p-values as evidence of a stage-specific peak.

Exact recommended Results language:

> In mouse meiotic Hi-C, subtelomeric sequence similarity was positively
> correlated with inter-chromosomal contact across leptotene, zygotene,
> pachytene and diplotene. At 1 Mb and 50 kb resolution, exact PHR-pair
> correlations were similar across stages (Spearman rho 0.37-0.43), while
> arm-pair-collapsed correlations were higher (rho 0.57-0.71) with the largest
> point estimate at zygotene. Direct clustered contrasts did not resolve a
> zygotene-specific peak, so we treat the mouse data as evidence for a broad
> prophase-I sequence-to-contact association rather than a stage-specific
> maximum.

### F_ST and cM/Mb

F_ST supported status:

- F_ST is present, not absent. The upstream/source trail includes the
  `/moosefs/.../heterogeneity/fst_superpop_matrix.tsv` output, Fig. 2 source
  paths, `scripts/ci/fst_block_jackknife.tsv`,
  `scripts/ci/fst_per_arm_per_pair.tsv`, and
  `scripts/popgen/matched_fst_d_m6.py`.
- The safe interpretation is background ancestry preservation: AFR/non-AFR
  subtelomeric F_ST values around 0.10-0.15 are statistically indistinguishable
  from matched 1000G/HGDP genome-wide autosomal Hudson F_ST; no pair is elevated
  above matched background.

F_ST recommended role:

> Keep only as demoted Methods/Extended Data context if Fig. 2 population panels
> remain. It should not be promoted as a main causal result or as evidence for
> ongoing exchange.

cM/Mb supported status:

- The recombination-rate input is Lalli 2025 T2T-CHM13 short-read cM/Mb.
- The join to PHR cross-arm affinity and sequence-community assignments is our
  analysis.
- The apparent negative correlation collapses after excluding low-callability
  acrocentric/PAR arms, so current short-read recombination maps cannot test
  the recombination-rate relationship at PHRs.

cM/Mb recommended role:

> Clarify if retained; otherwise cut. Do not phrase the anti-correlation as
> "reported by Lalli 2025" unless the sentence also states that the correlation
> and low-callability filter are our reanalysis using Lalli's recombination map.

## Data Provenance Status

This table records what the inspected artifacts support now. "Internal source"
means the repo points to `/moosefs` or cached repo TSVs; "public source" means a
citation is present in the manuscript/bibliography trail. Exact accession
numbers were not consistently present in the inspected artifacts, so downstream
source-data export must not invent them.

| Dataset / analysis | Current internal source status | Public citation/source note status | Provenance action before submission |
|---|---|---|---|
| HG002 Pore-C pointwise Fig. 4A | Repo snapshot `data/human_HG002_porec_50000bp_seqlevel.tsv`; upstream `/moosefs/.../analysis/human/community_free/human_HG002_porec_50000bp_seqlevel.tsv`; Fig. 4B matrix sources in `paper_prep/figures/fig3/sources.tsv` and B0. | `submission/paper.tex` cites `Ulahannan2019` for Pore-C method/use; A2 bibliography audit marks `Ulahannan2019` as present. | Source-data table should map internal processed TSVs to public raw-read dataset/citation and repo-relative processed file. Exact accession not found in inspected artifacts. |
| HG002 Hi-C pointwise ED1 | Repo snapshot `data/human_HG002_hic_50000bp_seqlevel.tsv`; upstream sequence-level sweep under `/moosefs/.../analysis/human/community_free/`. | Manuscript methods cite Hi-C processing references (`hic3d_dixon2012`, `hic3d_imakaev2012`) but inspected artifacts did not expose a specific HG002 Hi-C accession. | Add public raw-data accession/source note in source-data table; keep repo processed TSV path for reproducibility. |
| HG002 CiFi ED1 / support | Repo snapshot `data/human_HG002_cifi_50000bp_seqlevel.tsv`; upstream `/moosefs/.../analysis/human/community_free/`; community-based sources in `paper_prep/figures/fig3/sources.tsv`. | `submission/paper.tex` cites `hic3d_cifi2025`; A2 confirms bibliography key. | Pair citation with exact raw-read accession/source note in source-data table if available; not found in inspected artifacts. |
| GM12878 Dip-C | Internal outputs under `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/`; 16-cell summary in `end-to-end-report/report/06_dipc_validation.md`. | `submission/paper.tex` cites `Tan2018`; A3 recommends using the public Dip-C source instead of internal report line pointers. | Replace internal report pointers with `Tan2018` and public accessions/source-data rows. Exact accession not found in inspected artifacts. |
| PBMC Dip-C negative control | Internal outputs under `/moosefs/.../dipc_t2t/pbmc_hg19/enrichment_corrected/`; hg19 projection caveat recorded in reports. | A3 provides public-facing wording: "18 cells from the public Dip-C dataset `Tan2018`". | Keep only if single-cell controls remain; cite `Tan2018`, not internal report line numbers. |
| Sperm scHi-C | Internal outputs under `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/`; 20-cell results in `end-to-end-report/report/06_dipc_validation.md`. | `submission/paper.tex` cites `Xu2025` plus scHi-C methods references; A2 confirms keys. | Add exact public accession/source note in source-data table; not found in inspected artifacts. Include caveat that arm-level community-free aggregate is weak/wrong-direction. |
| Mouse zygotene / stage-resolved meiotic Hi-C | Cached F3 TSVs under `paper_prep/manuscript_revision/F3_mouse_shape/`; repo data under `data/mouse_meiosis_sweep/`; upstream `/moosefs/.../mouse_T2T/...`; Fig. 4 legacy manifest points to Zuo 2021 zygotene source tables. | `submission/paper.tex` cites `Francis2025` for B6/CAST T2T assemblies and `Zuo2021` for mouse meiotic Hi-C; A3 notes GEO accessions should be listed in a source-data table but they are not present in the inspected artifact text. | Source-data table should name B6/CAST T2T assemblies, Zuo 2021 mouse meiotic Hi-C, GEO accessions if available from source-data table, and repo processed stage tables. |
| RPE-1 Pore-C/CiFi | Internal outputs under `/moosefs/.../RPE1_subtelo/` and `/moosefs/.../analysis/human/community_based/RPE1/`; B0 and `05_hic_validation.md` summarize. | Public source/citation not fully resolved in B/F artifacts. | Keep as optional supplement only if source-data provenance is resolved; otherwise cut from active story. |
| Lalli cM/Mb map | Internal table `/moosefs/.../recombination_maps/subtelomeric_recomb_rates.tsv`; `paper_prep/figures/ed8/sources.tsv` records "Lalli 2025 T2T-CHM13 subtelomeric recombination rate". | `submission/bibliography.bib` contains `Lalli2025`. | Clarify that Lalli supplies the cM/Mb input; our analysis performs the PHR-affinity join and low-callability filter. |

## Recommended Figure and Table Moves

Main figure / active manuscript:

- Keep Fig. 4A as the primary human pointwise scatter, but remove pointwise
  p-values from caption and generated labels.
- Keep Fig. 4B as matrix visualization support, reporting B/W and O/E
  normalization plainly. Do not headline the astronomical W/B/Mann-Whitney
  p-value.
- Keep Fig. 4C mouse only if caption/text are softened to broad prophase-I
  association; remove pointwise p-values and "resolved zygotene peak" language.
- Keep ED1 human replicate scatters if active, but remove p-values from captions
  and labels.

Supplement / Methods:

- Move or summarize Mantel row+column permutation tables as the inferential
  support family; report permutation floors.
- Move multi-resolution W/B/Mantel grids, exclusion walks, W/B bootstrap
  details, and O/E details to compact Methods/Supplement source-data tables.
- Demote Dip-C/sperm/RPE-1 breadth to supplementary robustness unless the final
  author direction explicitly wants a broad 3D-validation figure.
- If F_ST panels remain in Fig. 2 or Extended Data, add matched-background
  caveat in legend/Methods. Do not add a new main Results claim.
- If cM/Mb remains, make it a Discussion/Limitations honest-null sentence, not a
  positive evidence panel unless the figure is explicitly an "honest null".

Source-data tables:

- Create submission-facing source-data tables with columns like
  `panel`, `public_source`, `public_accession_or_citation`,
  `repo_relative_processed_file`, and `note`.
- Do not submit internal `/moosefs` roots, `SURVEY_*` paths, or
  `end-to-end-report/...` line pointers as public source citations.

## Exact Wording Blocks For Downstream Patch

### Human 3D Results

> In human contact-map data, aggregate contact at PHR coordinates increased
> with PHR sequence similarity. We measured this per inter-chromosomal PHR pair
> within a single genome and without community labels in HG002 Pore-C
> (descriptive pointwise Spearman rho = 0.381, n = 2,830), with the same
> direction in CHM13 Hi-C, HG002 Hi-C and HG002 CiFi. Because highly similar
> non-homologous PHR tracts cannot be assigned uniquely by a MAPQ threshold, we
> interpret PHR-window contacts as aggregate coordinate-level measurements.
> Adjacent centromere-ward flanking regions, which are more uniquely mappable,
> show the same positive direction in key datasets, arguing against broad
> inflation from retaining multi-mappers.

### Fig. 4A Caption

> Contact rises with similarity (descriptive pointwise Spearman rho = 0.381,
> n = 2,830; line, linear fit).

### Fig. 4B Caption

> HG002 Pore-C contact matrix ordered by sequence community after
> observed/expected inter-chromosomal normalization; within-community blocks
> dominate (B/W = 0.056). Matrix-level inference used row+column
> Mantel/permutation tests in Methods.

### Mouse Results / Caption

> Mouse meiotic Hi-C supports a broad prophase-I sequence-to-contact
> association; the arm-pair-collapsed point estimate is largest at zygotene, but
> clustered zygotene-versus-flanking-stage contrasts do not establish a
> statistically resolved zygotene-specific peak.

### Methods Statistic Caveat

> Pointwise Spearman correlations are reported as descriptive effect sizes
> because PHR-pair observations share chromosome arms, PHRs, matrix
> rows/columns and read-placement histories. Matrix-level inference used
> row+column Mantel/permutation tests; finite permutation floors are reported
> when no permutation exceeds the observed statistic.

### F_ST If Retained

> Hudson F_ST was computed on a binary cross-arm/self-arm structural-haplotype
> state across the significant arm/community blocks. AFR/non-AFR values are
> approximately 0.10-0.15 but are statistically indistinguishable from matched
> 1000G/HGDP genome-wide autosomal Hudson F_ST, so this result is best read as
> ancestry-preservation/background-population-structure context rather than a
> subtelomere-specific differentiation signal.

### cM/Mb If Retained

> A PHR-level comparison of Lalli 2025 T2T-CHM13 short-read cM/Mb estimates
> with our cross-arm affinity metric shows an apparent anti-correlation only
> before filtering; it collapses after seven low-callability acrocentric/PAR
> arms are excluded, so current short-read recombination maps cannot test the
> recombination-rate relationship at PHRs.

## Unresolved J Decisions

These remain author/J-task decisions and should not be silently applied by an
integrator:

1. Evidence ordering: direct PHR-internal pointwise scatter first with flanking
   control immediately after (recommended), or flanking control first as the
   opening sentence.
2. Abstract verb strength for human 3D: keep "show" only if the MAPQ caveat and
   flanking control are adjacent; otherwise soften to "support".
3. Strict-MAPQ claim: locate/regenerate `comparison_summary.tsv` and source-data
   rows before retaining a quantified strict-MAPQ sentence, or replace with the
   workflow-only/flanking-control wording above.
4. W/B p-values: decide whether any W/B bootstrap/Mann-Whitney p-values remain
   in captions/tables. Recommendation: B/W as effect size only in captions.
5. Mann-Whitney global test: decide whether any author wants this retained as
   legacy provenance. Recommendation: cut from manuscript-facing apparatus.
6. Exclusion Mantel walk: decide whether a compact supplement table exists; if
   not, reduce to one Methods sentence or cut.
7. O/E language: keep plain "observed/expected inter-chromosomal normalization";
   avoid "random-ligation inflation" unless a methods owner endorses that exact
   model.
8. Mouse claim strength: replace "peaks at zygotene" in abstract/Results/Methods
   with broad prophase-I association unless authors explicitly accept the
   descriptive-only caveat.
9. Single-cell breadth: decide whether Dip-C/sperm/RPE-1 are retained as a
   compact robustness package or demoted/cut from the active narrative.
10. F_ST placement: if Fig. 2 population panels remain, add the matched-control
    caveat; otherwise cut the F_ST confidence-interval promise from the
    limitations/statistics sentence.
11. cM/Mb limitation: clarify as our analysis using Lalli 2025 input or remove
    the clause.
12. Source-data/accession table: complete public accession/source notes for
    HG002 Hi-C/Pore-C/CiFi, Dip-C, sperm scHi-C, and mouse meiotic Hi-C before
    submission-facing data availability is finalized.

## Validation Against Task Requirements

- Source-data citations/source notes are separated from internal paths in the
  provenance table; accession gaps are explicitly marked rather than invented.
- Primary versus supplementary evidence is classified in the evidence hierarchy.
- Recommended figure/table moves are listed for Fig. 4A/B/C, ED1, Mantel,
  single-cell/RPE-1, F_ST and cM/Mb.
- MAPQ0/random-placement language is provided only where supported by
  `end-to-end-report/report/05_hic_validation.md` and
  `end-to-end-report/report/06_dipc_validation.md`.
- Provenance status is recorded for HG002 Hi-C, HG002 Pore-C, HG002 CiFi,
  GM12878/PBMC Dip-C, sperm scHi-C, mouse zygotene/stage-resolved Hi-C, RPE-1,
  and Lalli cM/Mb.
- Unresolved J decisions are surfaced explicitly.
- Data-supported results are separated from proposed framing and exact wording.
