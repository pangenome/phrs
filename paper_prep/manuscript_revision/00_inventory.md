# Manuscript Revision 00 Inventory

Date: 2026-06-17

Scope: inventory only. I used `rg`, `find`, `sed`, and file metadata. I did not
run figure scripts, LaTeX builds, R data loads, Slurm jobs, or stream the
sequence-level `similarity.tsv.gz`.

## Active Manuscript And Build

- Active manuscript: `submission/paper.tex`.
- Active bibliography: `submission/bibliography.bib`.
- Build directory: `submission/`.
- Build commands: `cd submission && make` or `cd submission && bash compile.sh`.
- Build script sequence: `pdflatex paper`, `bibtex paper`, `bibtex Meth`,
  `bibtex Supp || true`, then two more `pdflatex paper` passes.
- Bundled template files: `submission/jnl.cls`, `submission/mathphys.bst`.
- Build log: `submission/BUILD_LOG.md` records a successful 2026-05-27 build,
  but `submission/README.md` notes that the build log predates the 5-figure
  restructure and is stale.
- Bibliography structure in the active TeX: one `.bib` file is used with
  normal `\cite{...}` for the main bibliography and `\citeMeth{...}` for the
  Methods References multibib section. `compile.sh` also runs `bibtex Supp`,
  tolerating an empty Supplementary bibliography.

## Prompt Status And Gate Rules

- Revision prompt: `wg_manuscript_revision_prompt.md`.
- Git status: the prompt is currently untracked (`?? wg_manuscript_revision_prompt.md`).
  That is a prompt provenance risk for downstream workers unless it is expected.
- The prompt states that task clusters contain two kinds of work:
  mechanical tasks `(M)` and judgment tasks `(J)`.
- J tasks must not be auto-resolved. They require a named author decision,
  recorded in writing, before an agent edits, cuts, moves, or rewrites content.
- This inventory does not resolve any J task. It only maps source files,
  candidate dependencies, and likely compute requirements.

## Concrete Current Manuscript Anchors

- Abstract count/framing issue: `submission/paper.tex` currently says
  "466 near-complete assemblies (465 HPRC v2 haplotypes together with CHM13)"
  and frames the q-arm sextet as an identified linked group.
- Intro/results count issue: `submission/paper.tex` currently says 233
  individuals and 465 HPRC v2 haplotypes plus CHM13.
- Methods count issue: `submission/paper.tex` currently says
  "465 haplotypes (233 samples x 2, minus 1 because CHM13 contributes a
  single haploid)", which is arithmetically inconsistent.
- Sampling language appears in the Results and Methods and again in
  Limitations.
- Fig. 2 caption currently labels a UPGMA-ordered heatmap as the
  "`phylogeny' of subtelomeric arms" and reports Leiden k=15 ordering.
- Methods currently include a `Neighbour-joining tree and character-level
  bootstrap` subsection.
- Fig. 4 caption and figure scripts currently report pointwise Spearman
  statistics for human PHR-pair contacts and mouse zygotene contacts.
- Methods currently include multi-mapping/MAPQ, flanking-control,
  strict-MAPQ, exclusion-control, Mantel, Mann-Whitney, bootstrap, and
  observed-over-expected language.
- Limitations currently mention an FST block-jackknife CI even though no active
  main result or figure appears to carry an FST analysis.
- Data availability currently contains internal `/moosefs/...` paths.

## Figure And Asset Inventory

- Active manuscript figure assets live in `submission/fig/MainFigures/` and
  `submission/fig/ExtendedDataFigures/`.
- Active main-figure generators that use repo-local `data/` snapshots live in
  `submission/scripts/figures/`:
  - `make_fig2bc_jaccard_heatmaps.R`
  - `make_fig4a_human_scatter.R`
  - `make_fig4b_porec_community.R`
  - `make_fig4c_mouse_zygotene.R`
  - `make_ed1_human_contacts.R`
- `submission/scripts/figures/README.md` says 9 figure panels are computed
  from `data/` and 7 images are vendored. Vendored items include the Fig. 3
  browser panels and Fig. 5 pedigree untangle PDF.
- Paper-prep figure directories still exist under `paper_prep/figures/` with
  `sources.tsv`, captions, scripts, and rendered PDFs/PNGs. Some of these are
  older/reviewer-era figure packages and do not all match the active
  5-main-figure manuscript structure.
- Notable current assets:
  - Fig. 1: `submission/fig/MainFigures/Fig1a_genomewide.pdf`,
    `submission/fig/MainFigures/Fig1b_lengths.pdf`.
  - Fig. 2: `submission/fig/MainFigures/Fig2a_pggb_layout.png`,
    `submission/fig/MainFigures/Fig2bc_jaccard.{pdf,png}`.
  - Fig. 3: six UCSC browser PNGs in `submission/fig/MainFigures/`.
  - Fig. 4: `Fig4a_human_scatter.{pdf,png}`,
    `Fig4b_porec_community.{pdf,png}`,
    `Fig4c_mouse_zygotene.{pdf,png}`.
  - Fig. 5: `submission/fig/MainFigures/Fig5_pedigree_untangle.pdf`.
  - ED1: `submission/fig/ExtendedDataFigures/ED_Fig1_human_contacts.{pdf,png}`.

## Cluster Mapping: Sources And Dependencies

### Cluster A: Mechanical Reconciliation

Likely no-compute or light-text/data tasks.

- A-1 sample counts:
  - Active text: `submission/paper.tex`.
  - Citation entry: `submission/bibliography.bib`, HPRC Data Release 2 entry
    states 232 individuals / 464 haplotypes + CHM13 in its note.
  - Cached input tables that expose counts: `data/all-vs-all.1Mb.p95.id95.len.tsv`
    and `paper_prep/figures/ed1/sources.tsv`.
  - Do not assume 232 vs 233; prompt requires confirmation of the release
    actually run.
- A-2 duplicated references:
  - Active citation source: `submission/paper.tex`.
  - Active bib database: `submission/bibliography.bib`.
  - Build mechanics: `submission/compile.sh`, `submission/README.md`.
- A-3 construction artifacts:
  - Active text: `submission/paper.tex`.
  - Known leaked paths: `/moosefs/guarracino/HPRCv2/PHR_III/`,
    `/moosefs/erikg/phrs/`.
  - Public repo target: `https://github.com/ekg/phrs`.
- A-4 ref [10] DOI:
  - Active bib: `submission/bibliography.bib`.
  - This likely requires web verification from official/preprint sources before editing.
- A-5 XTR acronym:
  - Active text and Fig. 1 caption: `submission/paper.tex`.
  - Fig. 1 sources: `paper_prep/figures/fig1/sources.tsv`,
    `scripts/plot-impg-coverage.R`.
- A-6 empty required sections:
  - Active sections: `submission/paper.tex`.
- A-7 symbol collisions:
  - Active Methods: `submission/paper.tex`.
  - Community detection inputs: `data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`,
    `data/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`,
    external `/moosefs/.../similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`.

### Cluster B: 3D-Contact Control Apparatus

Mixed M and J. B-1, B-2, B-3, B-5 are J and must not be auto-resolved. B-4 has
a mechanical detection phase, but final reporting is J.

- Active text/captions: `submission/paper.tex`.
- Current Fig. 4A generator: `submission/scripts/figures/make_fig4a_human_scatter.R`.
  It reads `data/human_HG002_porec_50000bp_seqlevel.tsv` and reports pointwise
  Spearman p-values.
- Current ED1 generator: `submission/scripts/figures/make_ed1_human_contacts.R`.
  It reads the 50 kb human contact seq-level snapshots.
- Current Fig. 4B generator: `submission/scripts/figures/make_fig4b_porec_community.R`.
  It reads `data/hg002_porec_contact_matrix.tsv`,
  `data/hg002_porec_hic.arm-leiden.communities.tsv`,
  `data/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`, and
  `data/hg002_porec_global_test.tsv`.
- Human seq-level sweep: `data/human_seqlevel_sweep/`.
  `data/human_seqlevel_sweep/README.md` documents the schema and upstream
  `/moosefs/.../analysis/human/community_free/` source.
- Human resolution/bootstrap script: `scripts/human/human_seqlevel_resolution.R`.
  It uses repo-local `data/human_seqlevel_sweep/` and PHR-node bootstrap; this
  is light enough in principle but should be scheduled deliberately, not run in
  this gate.
- MAPQ-strict rerun script: `scripts/hic/mapq_strict_d_peerq1.py`.
  This requires upstream pairs/cooler data under `/moosefs/.../3d/` and is a
  system/data task, not a head-node quick check.
- Cached 3D control tables:
  - `data/hg002_porec_global_test.tsv`
  - `data/hg002_porec_contact_matrix.tsv`
  - `data/chm13_phr_pair_correlation.tsv`
  - `data/hg002_porec_phr_pair_correlation.tsv`
  - `data/human_HG002_porec_50000bp_seqlevel.tsv`
  - `data/human_HG002_hic_50000bp_seqlevel.tsv`
  - `data/human_HG002_cifi_50000bp_seqlevel.tsv`
  - `data/human_CHM13_hic_50000bp_seqlevel.tsv`
- External 3D apparatus sources in manifests:
  - `paper_prep/figures/ed5/sources.tsv` points to
    `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/{5000,10000,20000,50000,100000}bp/<sample>_global_test.tsv`.
  - It also points to no-acrocentric controls, O/E matrices, and community
    bootstrap tests under `/moosefs/.../analysis/human/`.

Concrete validation paths requested by task:

- B-4 pointwise p-values:
  - `submission/paper.tex`
  - `submission/scripts/figures/make_fig4a_human_scatter.R`
  - `submission/scripts/figures/make_ed1_human_contacts.R`
  - `submission/scripts/figures/make_fig4c_mouse_zygotene.R`
  - `data/human_seqlevel_sweep/*.tsv`
  - `data/mouse_meiosis_sweep/seqlevel/**/*.tsv`
- B-5 remaining apparatus:
  - `submission/paper.tex`
  - `paper_prep/figures/ed5/sources.tsv`
  - `scripts/hic/mapq_strict_d_peerq1.py`
  - `scripts/ci/bootstrap_ci_d_m12.py`
  - `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/`
  - `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/no_acrocentric/`

### Cluster C: Communities, NJ Tree, Bootstrap, Continuum

All C tasks are J and must not be auto-resolved. C-0 requires author choice of
the continuum statistic and likely cluster execution if the full sequence-level
matrix is used. C-1 requires a named author statement about whether the tree is
load-bearing. C-2 requires code/table inspection if the bootstrap survives.

- Active text/captions: `submission/paper.tex`.
- Active Fig. 2 heatmap generator: `submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R`.
- Active Fig. 2 assets: `submission/fig/MainFigures/Fig2bc_jaccard.{pdf,png}`.
- Arm-order/support/audit tables:
  - `submission/fig/MainFigures/arm_order_tree.tsv`
  - `submission/fig/MainFigures/arm_order_community.tsv`
  - `submission/fig/MainFigures/community_blocks.tsv`
  - `submission/fig/MainFigures/arm_inclusion_audit.tsv`
  - `submission/fig/MainFigures/source_audit.tsv`
- Arm-level similarity/community cached tables:
  - `data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`
  - `data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
  - `data/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`
- External arm-level inputs:
  - `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-upgma-k14.assignments.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`
- NJ tree package:
  - `paper_prep/figures/nj_tree_arms/README.md`
  - `paper_prep/figures/nj_tree_arms/nj_tree.R`
  - `paper_prep/figures/nj_tree_arms/nj_tree.newick`
  - `paper_prep/figures/nj_tree_arms/nj_tree_annotated.{pdf,png}`
- Character bootstrap:
  - `scripts/cladistics/char_bootstrap_d_m9.R`
  - It reads a per-PHR involvement TSV and a reference arm distance matrix,
    resamples PHR rows, and rebuilds an arm-level distance.
  - Default paths in the script are stale local desktop paths, but the header
    names the canonical PHR TSV as
    `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv`.
  - A repo-local snapshot exists as `data/all-vs-all.1Mb.p95.id95.len.tsv`.
- Heavy sequence-level input:
  - `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.similarity.tsv.gz`
  - `paper_prep/figures/ed2/sources.tsv` reports this as about 12.5 GB.
  - Do not stream this on the head node.

Concrete validation paths requested by task:

- C-0:
  - `data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`
  - `data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
  - `submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R`
  - `submission/fig/MainFigures/Fig2bc_jaccard.pdf`
  - Heavy full matrix path above, requiring Slurm if used.
  - Existing 0.1-3.0 Leiden scan was referenced in the prompt but I did not
    find a clearly named repo-local scan table during this lightweight pass;
    likely upstream under `/moosefs/.../similarity/` or in slide revision assets.
- C-2:
  - `scripts/cladistics/char_bootstrap_d_m9.R`
  - `data/all-vs-all.1Mb.p95.id95.len.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv`
  - `data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`

### Cluster D: Sampling Justification

D-1 is J and must not be auto-resolved. It may have a mechanical source-audit
phase, but the final wording and whether to add a stability check require an
author decision.

- Active text: `submission/paper.tex`.
- Current statements mention 11.6%/12% wfmash sampling, Erdős-Rényi
  connectivity, 230x threshold, and a Limitations caveat.
- Relevant Methods subsection: pangenome graph and Jaccard similarity in
  `submission/paper.tex`.
- Cached source data:
  - `data/all-vs-all.1Mb.p95.id95.len.tsv`
  - `data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`
  - `data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
- Potential stability-check dependencies:
  - full similarity matrix under `/moosefs/.../pggb/...similarity.tsv.gz`
  - sequence-level assignment table under `/moosefs/.../similarity/hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv`
  - any existing Leiden resolution/subsampling scans under `/moosefs/.../similarity/`.

Concrete validation path requested by task:

- D-1:
  - `submission/paper.tex`
  - `data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`
  - `data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/pggb/...similarity.tsv.gz` if new
    sampling/subsampling checks are requested.

### Cluster E: Pedigree

All E tasks are J and must not be auto-resolved.

- Active text/caption: `submission/paper.tex`.
- Active Fig. 5 asset: `submission/fig/MainFigures/Fig5_pedigree_untangle.pdf`.
- Active vendored-generator note: `submission/scripts/figures/README.md` says
  Fig. 5 is copied from `end-to-end-report/pedigree-plots/washu/...untangle.pdf`
  and its generator needs off-repo moosefs data.
- Pedigree source survey: `paper_prep/surveys/SURVEY_14_pedigree_recombination.md`
  if present; `paper_prep/figures/fig4/sources.tsv` points to it.
- Pedigree script: `scripts/pedigree/monte_carlo_null_d_m4.py`.
- Upstream patch tables in manifests:
  - `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/all_pedigrees_patches.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle/recombination/patches.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/ceph1463/hifiasm/untangle/recombination/patches.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/ceph1463/verkko/untangle/recombination/patches.tsv`

### Cluster F: Orphan Defenses And Vestigial Elements

F-1 and F-3 require J resolution after mechanical detection. Do not auto-resolve
whether to delete, restore, or reinterpret these analyses.

- F-1 FST:
  - Active orphan limitation: `submission/paper.tex`.
  - Cached outputs: `scripts/ci/fst_block_jackknife.tsv`,
    `scripts/ci/fst_per_arm_per_pair.tsv`, `scripts/ci/results_d_m12.json`.
  - Scripts: `scripts/ci/bootstrap_ci_d_m12.py`,
    `scripts/popgen/matched_fst_d_m6.py`.
  - Older/source figure manifest: `paper_prep/figures/fig2/sources.tsv` points
    to `/moosefs/.../heterogeneity/fst_superpop_matrix.tsv` and
    `/moosefs/.../heterogeneity/cross_arm_superpop_enrichment.tsv`.
- F-2 cM/Mb:
  - Active limitation line: `submission/paper.tex`.
  - I did not find a concrete repo-local script in this pass beyond the text;
    this needs follow-up search in reports/history or author clarification.
- F-3 mouse Mantel double-reporting:
  - Active text/caption: `submission/paper.tex`.
  - Current Fig. 4C generator: `submission/scripts/figures/make_fig4c_mouse_zygotene.R`.
  - Mouse Mantel scripts: `scripts/mouse/mantel_d_m5.py`,
    `scripts/mouse/mantel_d_m5.R`.
  - Mouse shape/contrast scripts: `scripts/mouse/mouse_significance.R`,
    `scripts/mouse/mouse_stage_resolution_grid.R`.
  - Cached stage tables:
    - `data/mouse_leptotene_phr_20000bp_seqlevel.tsv`
    - `data/mouse_zygotene_phr_20000bp_seqlevel.tsv`
    - `data/mouse_pachytene_phr_20000bp_seqlevel.tsv`
    - `data/mouse_diplotene_phr_20000bp_seqlevel.tsv`
    - 50 kb siblings in `data/`.
  - Sweep tables: `data/mouse_meiosis_sweep/seqlevel/{1Mb,2Mb,4Mb}/`.
  - Upstream arm-level Mantel input path in script header:
    `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/subtelo_1Mb/similarity/mouse.dist_matrix.tsv`.
  - Upstream Hi-C matrix directory in script header:
    `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp`.

Concrete validation paths requested by task:

- F-1:
  - `submission/paper.tex`
  - `scripts/ci/bootstrap_ci_d_m12.py`
  - `scripts/ci/fst_block_jackknife.tsv`
  - `scripts/ci/fst_per_arm_per_pair.tsv`
  - `scripts/popgen/matched_fst_d_m6.py`
  - `paper_prep/figures/fig2/sources.tsv`
- F-3:
  - `submission/paper.tex`
  - `submission/scripts/figures/make_fig4c_mouse_zygotene.R`
  - `scripts/mouse/mantel_d_m5.py`
  - `scripts/mouse/mouse_significance.R`
  - `scripts/mouse/mouse_stage_resolution_grid.R`
  - `data/mouse_meiosis_sweep/seqlevel/`

### Cluster G: Title, Abstract, Introduction Alignment

All G tasks are J and must be last. They must not be auto-resolved until B-F
decisions are recorded.

- Active title/abstract/introduction: `submission/paper.tex`.
- Framing source prompt: `wg_manuscript_revision_prompt.md`, especially the
  two-tier continuum reframing and the heatmap-density question.
- Depends on outcomes from:
  - C-0/C-1/C-2 for the continuum/q-arm/tree/bootstrap claims.
  - B-1 through B-5 for 3D evidence ordering.
  - F-1 for whether "concerted evolution" has a live homogenization/FST result.
  - G-3 denominator clarity, also tied to A-1.

## No-Compute Judgment/Text Tasks Vs System/Data/Cluster Tasks

No-compute judgment/text tasks:

- Record author decision on the q-arm heatmap-density question before changing
  the abstract.
- Record author decision for B-1/B-2/B-3/B-5 before rearranging 3D evidence.
- Record author decision for C-1 on whether the NJ tree is load-bearing.
- Record author decision for C-3 on whether the q-arm list can be named.
- Record author decision for D-1 on whether to state connectivity-only or add a
  stability check.
- Record author decision for E-1/E-2/E-3 before changing pedigree claims.
- Record author decision for F-1/F-3 before deleting/restoring/reframing
  orphan analyses.
- Record author decision for G-0/G-1/G-2/G-3 before title, abstract, and
  introduction rewriting.

Mechanical text/source tasks that are light:

- Find and list all pointwise p-values in `submission/paper.tex` and figure
  scripts.
- Find internal mount paths and lab-notebook citations in `submission/paper.tex`.
- Find all count triples in `submission/paper.tex`, captions, and `bibliography.bib`.
- Check active citation keys and duplicate bibliography entries.
- Inspect small `sources.tsv`, README, script headers, and cached summary TSVs.

System/data/cluster tasks:

- Any analysis that streams the full sequence-level Jaccard file:
  `/moosefs/.../smooth.final.similarity.tsv.gz` (~12.5 GB).
- Any MAPQ-strict re-binning from `.allValidPairs`, `cooler cload`, or
  `cooler zoomify`.
- Any full re-run of community detection or subsampling stability on the
  sequence-level matrix.
- Any figure rebuild that loads large moosefs-only inputs through R/guix.
- Any Slurm-scale permutation/bootstrap beyond small cached repo-local tables.

## Small Cached Tables Vs Heavy Inputs

Small cached/repo-local tables suitable for lightweight inspection:

- `data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`
- `data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
- `data/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`
- `data/hg002_porec_contact_matrix.tsv`
- `data/hg002_porec_global_test.tsv`
- `data/hg002_porec_hic.arm-leiden.communities.tsv`
- `data/human_*_50000bp_seqlevel.tsv`
- `data/human_seqlevel_sweep/*.tsv`
- `data/mouse_*_phr_20000bp_seqlevel.tsv`
- `data/mouse_*_phr_50000bp_seqlevel.tsv`
- `data/mouse_meiosis_sweep/seqlevel/**/*.tsv`
- `scripts/ci/*.tsv`, `scripts/ci/*.json`
- `submission/fig/MainFigures/*audit.tsv`, `*order*.tsv`, `community_blocks.tsv`
- `paper_prep/figures/*/sources.tsv`
- `data/all-vs-all.1Mb.p95.id95.len.tsv` is 4.7 MB and appears safe for
  targeted header/line-count inspection, but avoid broad transformations in
  this gate.

Large or external inputs requiring Slurm/moosefs-aware handling:

- Full sequence-level similarity:
  `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.similarity.tsv.gz`
  (~12.5 GB in `paper_prep/figures/ed2/sources.tsv`).
- Upstream 3D pairs/mcool/cooler files under `/moosefs/.../3d/`.
- Human community-based multi-resolution controls under
  `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/`.
- Pedigree upstream patch and untangle-generation data under
  `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/`.
- Mouse upstream arm matrices and O/E matrices under
  `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/`.
- Figure scripts under `paper_prep/figures/` that declare guix manifests or
  moosefs-only inputs.

## Proposed Revision Artifact Convention

Use a single revision workspace:

- Root: `paper_prep/manuscript_revision/`
- Gate inventory: `00_inventory.md`.
- Author decisions: `paper_prep/manuscript_revision/decisions/`.
  Suggested filenames: `C0_continuum_decision.md`,
  `B3_3d_ordering_decision.md`, etc.
- Mechanical audits: `paper_prep/manuscript_revision/audits/`.
  Suggested filenames: `A1_counts.tsv`, `B4_pointwise_pvalues.tsv`,
  `A3_construction_artifacts.tsv`, `F1_fst_presence.md`.
- Lightweight derived tables from repo-local cached data:
  `paper_prep/manuscript_revision/tables/`.
- Slurm or moosefs-heavy run records:
  `paper_prep/manuscript_revision/slurm/<cluster-id>/` with command,
  input paths, output paths, git commit, and logs.
- Draft patches or prose proposals:
  `paper_prep/manuscript_revision/drafts/`, separate from `submission/paper.tex`
  until the relevant J decisions are recorded.
- Do not store large regenerated matrices or raw heavy outputs in git; store
  manifests and small summarized TSVs/MD instead.

## Risks And Blockers

- J-task blocker: most high-impact edits require named author decisions.
  The explicit J tasks that must not be auto-resolved are B-1, B-2, B-3, B-5,
  C-0, C-1, C-2, C-3, D-1, E-1, E-2, E-3, F-1 final disposition, F-3 final
  reconciliation, and G-0 through G-3. B-4 has a mechanical detection phase,
  but final reporting is J.
- Prompt provenance: `wg_manuscript_revision_prompt.md` is untracked in this
  worktree. Downstream tasks may miss the prompt if they run from a clean clone.
- Missing/discoverability issue: I did not find a clearly named repo-local
  Leiden resolution scan table for the prompt's 0.1-3.0 resolution evidence.
  It may be in moosefs, slide assets, or git history.
- Potential manuscript/source drift: `paper_prep/figures/fig4/sources.tsv`
  describes older pedigree/RPE-1/mouse panel content, while active
  `submission/` Fig. 4 is human scatter + Pore-C community + mouse zygotene.
  Downstream workers should treat `submission/` as active and `paper_prep/figures/`
  as provenance/older figure-package context unless an author says otherwise.
- README inconsistency: `submission/README.md` contains statements both that
  reviewer-era analyses were retained as text and that they were cut; active
  `submission/paper.tex` still contains some reviewer-era apparatus. Confirm
  against the prompt before using README prose as policy.
- Data availability risk: active manuscript still names internal `/moosefs/`
  paths.
- Count risk: active manuscript counts conflict with the HPRC release note in
  `bibliography.bib`. A-1 needs release/version confirmation, not just a typo fix.
- Head-node safety: do not stream the 12.5 GB sequence-level similarity file,
  do not run MAPQ strict re-binning, and do not run guix/R figure rebuilds with
  moosefs-heavy inputs on the head node.
