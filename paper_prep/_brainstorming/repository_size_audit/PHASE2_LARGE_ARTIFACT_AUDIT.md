# Phase 2 large-artifact audit

**Audit date:** 2026-07-17 UTC

**Audit target:** live `pangenome/phrs` refs at `2026-07-17T14:45:26Z`

**Status:** reviewable proposal only. **No cleanup action in this document has
been performed or authorized.**

## Executive finding

The July rewrite achieved its intended result. With one deterministic packing
method, the commit immediately before the copy-enrichment series packs to
505,421,184 bytes, current rewritten `origin/main` packs to 510,205,212 bytes,
and the unreduced surgery-time backup lineage packs to 1,005,907,764 bytes. The
unreduced July series therefore added 500,486,580 packed bytes; the rewritten
lineage adds only 4,784,028 bytes over the pre-series baseline. The phase-1
reduction relative to the unreduced lineage is 495,702,552 bytes.

A conservative second pass could remove generated proof/raster history,
reproducible Fig. 5 scratch renders, verified raw Fig. 2 duplicates, obsolete
LaTeX intermediates, and obsolete WG refs. It must **not** remove the canonical
CHM13 GFF3, compressed manuscript inputs, retained V7 outputs, directly used
submission figures, paper figure packages, or the final/frozen slide decks.
The candidate groups overlap, so their isolated costs must not be added to
promise a total saving.

## Methodology

1. I read live refs with `git ls-remote --heads --tags`; I did not fetch or
   change local or remote refs. The final stable origin snapshot had 30 heads:
   `main`, the archive head, and 28 WG heads. `origin/main` was
   `2d95f4d181760bab6547ef07c1da56b894aa823c`.
2. I identified the July copy-enrichment series from first-parent history. Its
   first rewritten commit is `726bc141` and its immediate parent is
   `0d9ce461105cb31a388a77afd63a47eb453b1b17`; the rewritten range through
   `d572d11e` has 19 commits. This is the pre-series comparison point, not an
   estimate based on dates.
3. I cloned `git@github.com:ekg/phrs-backup.git` read-only to a temporary bare
   directory. The backup's surgery-time `main` is `73087d341e53` (parent
   `6dbbd35e`); its current `main`, including the three archive-note commits, is
   `41f37be47c4b`. The backup itself was not modified.
4. Every comparable history baseline used Git 2.20.1 and the same command:

   ```text
   git pack-objects --stdout --revs --window=10 --depth=50 --threads=1 \
     --delta-base-offset <revisions | wc -c
   ```

   Values are exact bytes. “MiB” means bytes / 1,048,576. Path-group costs use
   the same pack options, but feed the sorted unique blob IDs for that exact
   path group without `--revs`. Those are isolated blob-pack costs: useful for
   ranking, but non-additive because groups overlap and a larger set can find
   better deltas.
5. HEAD presence came from `git ls-tree`; historical versions came from
   `git rev-list --objects origin/main -- <paths>`. Dependencies were searched
   with `rg`, and regeneration claims were checked against committed scripts,
   READMEs, manifests, and the actual `submission/paper.tex` image paths.
6. Exact duplicates were established by Git object ID and SHA-256. For gzip
   pairs I additionally compared `gzip -dc ... | sha256sum`; equal filenames,
   sizes, or biological meaning alone were not treated as duplicate evidence.
7. No `gc`, prune, repack-in-place, history rewrite, ref edit, deletion, push,
   or backup mutation was run. Temporary pack streams and the temporary backup
   clone were outside the repository.

## Comparable size baseline

| Reachable set | Tip/snapshot | Exact packed bytes | MiB | Interpretation |
|---|---|---:|---:|---|
| Commit before July series | `0d9ce461` | 505,421,184 | 482.007 | Common pre-enrichment baseline |
| Current rewritten origin main | `2d95f4d` | 510,205,212 | 486.570 | Public main history after phase 1 |
| All 30 live origin heads | live at timestamp above | 603,025,511 | 575.090 | Main plus archive/WG reachability |
| Backup surgery-time main | `73087d34` | 1,005,907,764 | 959.308 | Unreduced series plus archive note |
| Backup current main | `41f37be4` | 1,005,919,858 | 959.320 | Same, plus two small note fixes |
| All 30 backup heads | live at audit time | 1,102,290,072 | 1,051.226 | Unreduced main plus WG objects |

The all-origin full-pack delta over current main is 92,820,299 bytes. A second
calculation packing only the 194 non-main reachable blob IDs in isolation gave
100,885,869 bytes (241,817,714 uncompressed bytes). The numbers differ because
the full union can delta branch blobs against main blobs; the isolated branch
set cannot. The full-union delta is the appropriate whole-repository baseline,
while the isolated number describes the branch-only content itself.

### Why GitHub, Git, and `du` disagree

At audit time the GitHub repository API reported `size=594753`, conventionally
KiB: 609,027,072 bytes (580.813 MiB). That is close to, but 6,001,561 bytes
larger than, the reproducible all-origin pack. GitHub's value is server storage,
not a promised output of a specified `pack-objects` invocation; pack indexes,
server packing choices, ref/log overhead, and update lag can all contribute.
It should be used as an operational trend, not as the scientific baseline.

The shared local Git directory reported 610.52 MiB in packs plus 441.24 MiB of
loose objects via `git count-objects -vH`, while `du -sb` reported
1,837,010,898 bytes. It is shared by many WG worktrees and contains local refs,
reflogs, loose objects, and MooseFS directory metadata not reachable from live
origin. It was deliberately not garbage-collected for this audit.

The checked-out worktree had 918,831,666 bytes of tracked file contents
(849,554,858 bytes after counting identical blobs once), but `du -sb` reported
1,374,273,112 bytes. On this MooseFS checkout, directory entries themselves
contributed 455,441,423 apparent bytes; the 3,114 regular tracked files summed
exactly to 918,831,666 bytes. A checkout is also uncompressed, whereas a Git
pack uses zlib and cross-version deltas. These values answer different
questions and should not be compared as though they used the same units or
object set.

## Ranked path groups

The table is ranked approximately by isolated historical packed cost. Parent
and child rows intentionally overlap (for example, `data/` includes the GFF3,
and `slides/` includes PDF and proof-PNG subgroups).

| Exact path group | Category / purpose | HEAD pack bytes | Main-history pack bytes | Recommendation (proposal only) |
|---|---|---:|---:|---|
| `slides/` | Decks, frozen revision sources, proofs, and review assets | 153,291,883 | 181,588,038 | **Retain decks/sources; split generated proofs below** |
| `data/` | Vendored canonical analysis and figure inputs | 135,142,395 | 135,142,395 | **Retain**, except independently verified raw duplicates below |
| Four Fig. 5 lead directories: `fig5_homolog_vs_interchrom_whole_genome_ribbon_draft/`, `fig5_extended_maternal_whole_genome_ribbons/`, `fig5_whole_genome_alignment_overview/`, `fig5_whole_genome_length_scaled_tracks/` under `paper_prep/_brainstorming/` | Generated scratch renders, tables, and generators | 26,102,827 | 92,414,553 | **Remove only reproducible output subsets; retain scripts/tables used downstream** |
| Enumerated `slides/**/*.pdf` paths | Final decks plus supporting/revision PDFs | 70,878,930 | 70,429,095 | **Retain final/frozen decks; investigate supporting PDFs**. The historical pack is slightly smaller because extra versions improve delta selection |
| `paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft/*.{pdf,png,svg}` | Fig. 5 scratch render versions | 12,562,796 | 55,427,660 | **Propose remove render history**, after reproducing current Fig. 5 |
| Enumerated `slides/*/_typst/page-NN.png` paths | Typst proof PNGs (116 at HEAD; 210 unique historical blobs) | 27,217,380 | 54,082,810 | **Propose remove**, after link updates and clean deck rebuilds |
| `data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz` | Canonical CHM13 RefSeq/Liftoff annotation | 52,996,471 | 52,996,471 | **Retain / explicit exclusion** |
| `submission/fig/` | Manuscript-facing and retained submission assets | 22,633,894 | 46,200,148 | **Retain directly used files**; two unused legacy assets are investigated below |
| `paper_prep/_brainstorming/fig5_extended_maternal_whole_genome_ribbons/**/*.{pdf,png,svg}` | Maternal source-render variants | 8,303,482 | 29,680,415 | **Propose remove renders only**, not run tables/scripts |
| `paper_prep/_brainstorming/fig2a_community_colored/` | Raw Fig. 2 layout/community copies plus scratch renderer | 19,047,375 | 19,047,375 | **Propose remove the two raw tables after reference update** |
| `data/fig2a_pggb_layout.og.lay.tsv.gz` + `data/fig2a_node_community.tsv.gz` | Canonical compressed Fig. 2 inputs | 18,923,212 | 18,923,212 | **Retain / explicit exclusion** |
| `paper_prep/_brainstorming/pedigree_direct_sweepga_concordance/` | Review-only SweepGA PAF package | 8,431,126 | 8,431,126 | **Investigate**; raw/filtered PAF subset is 8,385,571 bytes but regeneration is external and costly |
| Historical `paper_prep/submission/` | Superseded manuscript drafts, absent at HEAD | 0 | 7,281,938 | **Investigate history-only removal**; low saving versus provenance value |
| `inter-chr-plots/` + `identity_heatmaps/` | Vendored canonical pre-rendered figure inputs | 5,350,537 | 5,350,537 | **Retain / explicit exclusion** |
| `paper_prep/figures/` | Canonical paper figure packages and outputs | 4,721,030 | 4,721,329 | **Retain / explicit exclusion** |
| `paper_prep/lit_review/_render/` | 38 generated PPM page intermediates | 4,326,304 | 4,326,304 | **Propose remove**; retain Typst, Markdown, bibliography, and PDF |
| `paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v7/` | Retained final V7 results, validation, community attribution | 3,829,016 | 3,829,540 | **Retain / explicit exclusion** |
| Six overview outputs in `fig5_whole_genome_alignment_overview/` | Generated PDF/PNG/SVG and three TSVs | 3,572,280 | 3,572,280 | **Propose remove renders; investigate generated TSV retention** |
| Fourteen output paths in `fig5_whole_genome_length_scaled_tracks/` | Generated render and track-table outputs | 1,543,658 | 3,136,228 | **Propose remove after regeneration check** |
| `paper_prep/_brainstorming/comprehensive_copy_background.csv` | Exploratory 58,230-row copy background | 412,199 | 412,199 | **Investigate**, not a proven duplicate and not identified as “CQB” |
| Historical `paper_prep/submission/*.{aux,bbl,blg,log,out,...}` | LaTeX build intermediates, absent at HEAD | 0 | 63,397 | **Propose remove from history** (subset of historical submission row) |

## Dependency, reproducibility, and retention evidence

| Candidate or protected item | At HEAD? | Reproduction/dependency evidence | Directly used by current submission? | If later removed, where protected? | Recommendation (proposal only) |
|---|---|---|---|---|---|
| Canonical CHM13 GFF3 | Yes | Default input in `build_inputs.py:633` and `build_term_maps.py:23`; documented as the sole physical annotation in `STATISTICAL_SPEC.md:43`; checksum frozen in `ontology_v6/INPUT_MANIFEST.tsv:2` | Not an image included by LaTeX, but it is a canonical analysis input | Origin and backup | **Retain** |
| Canonical compressed Fig. 2 tables | Yes | `submission/scripts/figures/make_fig2a_pggb_layout.R:40-41` reads them directly and writes the manuscript panel; fully repo-local | Yes, through reproducible `Fig2a_pggb_layout.png` | Origin and backup | **Retain** |
| Raw Fig. 2 layout/community tables | Yes | Content-identical to canonical gzip payloads. The long raw layout name remains in old render provenance and revision scripts; current manuscript generator uses compressed paths | No | Canonical gzip payloads remain; old versions also in backup | **Propose removal only after changing/archiving those stale raw-path references** |
| Fig. 5 combined PNG at scratch and `submission/fig/MainFigures/Fig5_whole_genome_recombination.png` | Yes | Exact duplicate. Current `paper.tex:515-521` instead uses Fig5A and separate Fig5B/C/D PNGs; `compose_main_fig5.py:287-292` can regenerate the combined scratch image | No, despite stale scratch README wording | Backup; current panel generator and canonical `data/fig5_*` inputs remain | **Propose remove both copies** |
| Other main Fig. 5 scratch renders | Yes | Generated by scripts retained beside their run tables. Historical render pack is large because many PDF/SVG versions accumulated | No | Backup and retained generators/tables | **Propose remove render extensions only** |
| Extended maternal renders | Yes | Build command is committed at `fig5_extended_maternal_whole_genome_ribbons/README.md:16-27`, but it currently resolves some inputs through absolute MooseFS/worktree paths (`:29-32`). Exact run TSVs are read by `submission/scripts/figures/make_fig5_pedigree_supporting_material.py:25-30` | The scratch PDF/PNG/SVG files are not included by `paper.tex`; run TSVs support submission material | Backup; scripts and run TSVs retained | **Propose remove renders only after a successful clean rerun; retain every TSV, manifest, note, and script** |
| `submission/fig/ExtendedDataFigures/ED_Fig2a_*` and `ED_Fig2b_*` maternal PNGs | Yes | Retained source assets described by the maternal README, but the manuscript has only ED1 and `paper.tex:775` includes `ED_Fig1_human_contacts.pdf` | No | Backup and maternal renderer | **Investigate removal**; resolve the README/source-asset policy first |
| Overview outputs | Yes | Inputs and commands are enumerated in `fig5_whole_genome_alignment_overview/README.md:8-26,72-100`; package is not referenced by `paper.tex` | No | Backup; generators and upstream manifests remain | **Propose renders; investigate TSVs**, because upstream filtered PAF availability must be rechecked in a clean clone |
| Length-track outputs | Yes | `fig5_whole_genome_length_scaled_tracks/README.md` supplies the build command and exact output names; not referenced outside scratch | No | Backup; script and manifests remain | **Propose removal after clean rerun** |
| Direct SweepGA raw/filtered PAFs | Yes | Explicitly review-only (`pedigree_direct_sweepga_concordance/README.md:1-6`). Inputs include absolute MooseFS FASTA/lists (`:7-19`), and rerun needs SweepGA/FastGA 0.1.1 (`:23-37`) | No; the evidence review references its summaries/report, not raw PAFs | Backup and external biological inputs, but not self-contained | **Investigate, do not automatically remove**; 8.39 MB is modest relative to rerun cost |
| Literature-review PPM pages | Yes | No repository reference to any `page-NN.ppm`; authoritative `SYNTHESIS_v2.typ` and the 406,728-byte `SYNTHESIS_v2.pdf` remain. PPM is a deterministic PDF rasterization intermediate | No | Backup; source Typst/PDF retained | **Propose remove** |
| Slide proof PNGs | Yes | Generated by the documented Typst `--ppi 144 ... page-{0p}.png` build. Revision notes cite the proofs (for example `REVISION_NOTES_V9.md:8`), and a few revision READMEs use pages as visual baselines | No | Backup; frozen Typst sources/assets and final PDFs retained | **Propose remove only after updating proof links and validating every frozen deck selected for preservation** |
| Historical `paper_prep/submission/` | No | Last touched at `b86ee424` on 2026-06-14. It has 98 historical blobs; only two (55,932 raw bytes) are also current `submission/` blob IDs, so this is old draft history, not an exact duplicate tree | No | After rewrite, private backup only | **Investigate**; provenance loss may outweigh 7.28 MB |
| Historical LaTeX intermediates | No | Standard products of `compile.sh`; 14 versions across nine paths, all absent at HEAD | No | Private backup only | **Propose remove** |
| V7 outputs | Yes | Phase-1 policy explicitly preserves final V7 result tables, validation, community attribution, and heatmaps (`ARTIFACT_ARCHIVE.md:13-15,32-34`). The 36 HEAD blobs include `TERM_RESULTS.tsv.gz`, exact attribution, scripts, manifests, tests, and reports | Not currently included by `paper.tex`; intentionally retained paper-facing analysis record | Origin and backup | **Retain** |
| Current manuscript figures | Yes | `paper.tex:395-521,775` resolves 19 images totaling 4,118,515 raw bytes; all exist. Fig. 2 and Fig. 5 generators read vendored `data/` inputs (`make_fig2a...:40-41`; `make_fig5_panels.sh:45-52`) | **Yes** | Origin and backup | **Retain all 19 exact referenced paths** |
| Canonical paper figures and pre-renders | Yes | Project figure contracts identify `paper_prep/figures/`, `inter-chr-plots/`, and `identity_heatmaps/` as paper-facing or vendored upstream artifacts | Some are upstream sources rather than direct LaTeX paths | Origin and backup | **Retain** |

The directly included LaTeX files are the ORCID icon; Fig1a/1b; Fig2a/2bc;
six Fig. 3 browser PNGs; Fig4a/4b/4c; Fig5A/B/C/D; and
`ED_Fig1_human_contacts.pdf`. This allowlist, derived from the actual TeX rather
than directory names, should be a hard rewrite test.

## Verified duplicate claims

| Paths | Git/SHA-256 evidence | Meaning for cleanup |
|---|---|---|
| `data/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.tsv` and `paper_prep/_brainstorming/fig2a_community_colored/layout.tsv` | Same Git blob `952f0ef1`; raw SHA-256 `2d73f55f807f0b037b69cf909cde8739d1320577a504f321571fb992a3b6a32f` | Removing only one path saves checkout space but no packed blob. Remove both raw paths to save the blob |
| Above raw layout and decompressed `data/fig2a_pggb_layout.og.lay.tsv.gz` | Decompressed SHA-256 is the same `2d73...a32f`; gzip-file SHA-256 is `6563b469...bcbe2` | Canonical gzip is a content-equivalent replacement |
| `fig2a_community_colored/node_comm.tsv` and decompressed `data/fig2a_node_community.tsv.gz` | Decompressed/raw SHA-256 `a5d0787c2b6178585d9e2e07e30cdc80ed855972324895f82fef70ec2366aa8d`; gzip SHA-256 `c0ccdc90...d16f1f2` | Canonical gzip is a content-equivalent replacement |
| Scratch and submission `Fig5_whole_genome_recombination.png` | Same Git blob `77a380a2`; SHA-256 `f6c2536aabe58e2976daee16d6031995545203ad1502b88faa50c3d4df7b0ea1` | Both are unused by current TeX; removing both releases the blob |
| V8 and V9 `nearest_same_superpop_mds_distances.tsv` in their named revision-asset directories | Same Git blob `3d5cc182`, 4,180,581 bytes | No packed saving while either frozen copy remains; do not claim two-object savings |
| Three raw `*.sweepga_many_many_j0.paf.gz` files and corresponding `filtered_paf/*.many_many_noscaffold.paf.gz` in `pedigree_direct_sweepga_joint_parent/` | Same Git blobs per pair. Compressed/decompressed SHA-256 prefixes: PAN027mat `dfb63c55` / `7501f8d8`; PAN027pat `9e5e0be4` / `21b1d542`; PAN028mat `c3218682` / `f606e4cd` | Exact path duplication, but Git already stores one blob per pair. Removing one alias provides no packed-history saving if the other remains |

The two maternal source PDFs are **not** duplicates of the manuscript Fig5C/D
PDFs: their full SHA-256 values differ (`c18f9b92...` versus `e0f960ff...` for
the PAN027 maternal comparison, and `1d73e085...` versus `dc22b5b2...` for
PAN028). They are generated variants, not hash duplicates.

The GFF3 gzip SHA-256 is
`a1c8e61cb4e60a3af3a18599b7d5551a72a1b0317bdffad42ae7fa36e73da968`;
its decompressed SHA-256 is
`286a8a79a3d940c389727e45b7c9994d126144a159f4fc47f7da95fa7e03f171`.
No second matching GFF3 blob was found.

## CHM13 and “CQB” check

No literal, case-insensitive `CQB` path component was found in any current
origin history or in the backup's refs. The only path matches from a broader,
explicit `copy.*background|comprehensive.*background` search were:

- `paper_prep/_brainstorming/build_genome_wide_copy_background.py`
- `paper_prep/_brainstorming/comprehensive_copy_background.csv`

I do **not** infer that either name expands or represents the acronym “CQB.”
The CSV is 2,555,317 raw bytes, packs alone to 412,199 bytes, and has SHA-256
`645e9d0c3844018fdc161c5a3d93961fee837f53130f7a843efd30d20fc30e9b`.
`SURVEY_DATA_inventory.md:34,114` calls it a strict eight-column superset of
`genome_wide_gene_copies.csv` and says one should be selected before SI freeze;
it is also read by `improved_copy_weighted_enrichment.R:16` and catalogued as
parked/methodologically caveated by the slide inventory. This is semantic
redundancy, not a hash duplicate. Recommendation: investigate during SI freeze,
not phase-2-remove by acronym.

## Origin-ref contribution

The archive head `d572d11e`, WG heads `d0444703` and `8f23c3d6`, and `main`
have no blobs outside current main. The largest individual non-main sets are
below. Costs are isolated packs of blobs reachable from that tip and not main;
they overlap, so they are not additive.

| Live origin ref | Tip | Non-main blobs | Isolated pack bytes | Backup protection |
|---|---:|---:|---:|---|
| `wg/agent-2688/fig5-raw-fasta-sweepga-f16-no-chop-merge-panels` | `6f0b9e36` | 25 | 39,493,534 | Exact tip also present in private backup |
| `wg/agent-2907/fig5-pre-impg-best-period-similarity` | `2943d58e` | 18 | 26,079,975 | Exact tip also present in private backup |
| `wg/agent-2949/fig5-whole-genome-full-homologous-chain-ribbon-layer` | `4f20d20d` | 8 | 14,129,949 | Exact tip also present in private backup |
| `wg/agent-944/fix-bog-review` | `759b6732` | 27 | 11,759,371 | Exact tip also present in private backup |
| `wg/agent-941/render-bog-annotated` | `f5f94956` | 22 | 8,407,118 | Exact tip also present in private backup |
| `wg/agent-2762/fig5-raw-manymany-impg-similarity-fullbed` | `ef24d4ae` | 24 | 6,197,338 | Exact tip also present in private backup |
| `wg/agent-2741/fig5-length-scaled-whole-genome-tracks` | `992ed4dd` | 94 | 5,298,220 | Exact tip also present in private backup |
| `wg/agent-2733/fig5-whole-genome-alignment-overview` | `b5ddc6f1` | 73 | 3,638,392 | Exact tip also present in private backup |

The remaining 18 contributing WG refs each pack to 888,647 bytes or less
(range 1,352–888,647). Their content is mainly small Fig. 5 experiments,
manifests, and CHM13/WG reports. Branch deletion without path cleanup would
remove branch-only objects only after the server makes them unreachable and
repackages; it would not shrink main. **Proposal only:** confirm each task is
closed/merged and delete obsolete public WG refs in a separate, reviewed ref
hygiene operation. Do not rewrite active task refs.

## Conservative phase-2 path/ref list

This is a proposal for review, **not an instruction already executed**.

### Tier A: strong path candidates

- `paper_prep/lit_review/_render/`
- The 116 current proof paths matching these exact directory-level patterns:
  `slides/v2/_typst/page-*.png`, `slides/v2-zoom/_typst/page-*.png`,
  `slides/v2-review/_typst/page-*.png`, and
  `slides/v2-review-zoom/_typst/page-*.png`; rewrite their historical versions
  only after updating committed links.
- The two raw Fig. 2 tables in
  `paper_prep/_brainstorming/fig2a_community_colored/` plus the long-named raw
  layout TSV under `data/`; retain both `data/fig2a_*.tsv.gz` files.
- Generated `*.pdf`, `*.png`, and `*.svg` files in
  `fig5_homolog_vs_interchrom_whole_genome_ribbon_draft/`, while retaining
  scripts, README/notes, run tables, and manifests.
- Both exact copies of `Fig5_whole_genome_recombination.png` (scratch and
  submission), because the current TeX uses the separate Fig5A-D files.
- Historical LaTeX intermediates under `paper_prep/submission/`:
  `Meth.aux`, `Meth.bbl`, `Meth.blg`, `Supp.aux`, `paper.aux`, `paper.bbl`,
  `paper.blg`, `paper.log`, and `paper.out` (plus any independently enumerated
  historical build-only siblings).

### Tier B: conditional candidates

- Generated render extensions in
  `fig5_extended_maternal_whole_genome_ribbons/`, and unused
  `submission/fig/ExtendedDataFigures/ED_Fig2a_PAN027_maternal_recombination.png`
  / `ED_Fig2b_PAN028_maternal_recombination.png`, only after resolving their
  source-asset policy and proving regeneration. Retain all maternal run TSVs.
- The six named overview outputs listed in its README and the fourteen named
  length-track outputs, only after a clean-clone rerun. Retain generators,
  manifests, validation scripts, and any compact table chosen as review
  evidence.
- `pedigree_direct_sweepga_concordance/{raw_paf,filtered_paf}/`: investigate
  moving large rerunnable PAFs to external/checksummed storage while retaining
  compact reports, summaries, manifests, and scripts. Do not remove merely
  because they are generated; the run is external and expensive.
- Historical `paper_prep/submission/`: consider only if authors accept losing
  public draft history and the private backup has been independently verified.
- The eight high-contribution WG refs in the origin table, followed by the
  smaller closed WG refs. Ref removal is a separate proposal, not a path filter.

### Explicit exclusions from phase 2

- `data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz`.
- `data/fig2a_pggb_layout.og.lay.tsv.gz`,
  `data/fig2a_node_community.tsv.gz`, the three `data/fig5_*.class_winners...gz`
  tables and their FAI files, and other direct submission generator inputs.
- All 19 image paths actually included by current `submission/paper.tex`.
- `paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v7/`.
- `paper_prep/figures/`, `inter-chr-plots/`, and `identity_heatmaps/`.
- Final/frozen deck PDFs, the paper-story PPTX, Typst sources, revision notes,
  and unique revision assets. Only enumerated proof PNGs are candidates.
- The already removed phase-1 enrichment directories listed in
  `ARTIFACT_ARCHIVE.md:18-30`; do not repeat that rewrite or re-import them.
- The private backup, all `/moosefs/...` source-of-truth data, and any active WG
  task ref.
- `comprehensive_copy_background.csv` until SI/data owners decide between it
  and the smaller table; no action is justified by the unsubstantiated acronym
  “CQB.”

## Full safety protocol for any later rewrite

Every step below is required before a future cleanup. Nothing here authorizes
the rewrite.

1. **Obtain explicit owner approval and freeze writes.** Name the approved path
   and ref list, announce a maintenance window, stop merges/WG dispatch, and
   record protected-branch and collaborator constraints. Abort if origin changes
   after the freeze.
2. **Capture immutable before-state evidence.** Record `git ls-remote` for every
   head/tag, GitHub API size/time, default branch, open PR bases, release/tag
   refs, submodule/LFS state, and SHA-256 of the approved filter specification.
   Save the 30-ref mapping in the private audit record, not as an unreviewed
   repository blob.
3. **Verify two recovery copies.** Run `git fsck --full` on the existing private
   backup and on a fresh second mirror; verify that `73087d34`, `41f37be4`, all
   candidate WG tips, tags, and every removal-candidate blob are readable.
   Make the recovery mirrors read-only. Never filter, prune, or GC them.
4. **Use a disposable clone only.** Pin and record the history-rewrite tool and
   version. Run the exact path/ref filter in a newly created bare clone. Never
   experiment in `/moosefs/erikg/phrs/.git`, a working checkout, or the backup.
5. **Separate path cleanup from ref cleanup.** Produce an explicit mapping of
   old to new commits for every rewritten public ref. For obsolete WG refs,
   obtain task-owner confirmation and list exact old tips; do not delete active
   refs or rely on broad `refs/wg/*` globs.
6. **Enforce retain and remove allowlists.** Assert that every excluded path
   above remains. Assert that only approved candidate paths disappear. Compare
   HEAD trees and file counts, allowing only the reviewed HEAD deletions; fail
   on any unrelated mode, symlink, filename, or content change.
7. **Re-run hash checks.** Confirm canonical Fig. 2 gzip payloads still have the
   decompressed hashes recorded above, the GFF3 hashes match, the V7 manifest
   validates, and all 19 TeX image paths exist. Verify that each claimed
   duplicate removed has a retained hash-identical or decompressed-identical
   counterpart where promised.
8. **Validate reproducibility before publishing.** In a clean checkout, run the
   submission build (`submission/compile.sh`) and inspect the final LaTeX log;
   regenerate Fig2a and Fig5A-D from vendored inputs; compile each preserved
   Typst deck; run overview/length/maternal validation if their outputs are
   removed; check links/manifests and repository tests. Missing MooseFS/tool
   prerequisites are a blocker, not permission to skip.
9. **Measure with the same baseline.** Pack rewritten main and the full proposed
   public ref set using the exact options in this report. Report exact bytes,
   object counts, and the non-additive nature of group estimates. Investigate
   any saving materially different from the dry-run forecast.
10. **Require independent review.** A second person must inspect the old/new ref
    map, retained-tree diff, removal list, hashes, builds, and recovery commands.
    Sign off on scientific provenance, manuscript assets, and frozen-deck policy
    separately from Git mechanics.
11. **Push explicitly and lease every ref.** With branch-admin approval, push
    only reviewed refs using per-ref
    `--force-with-lease=<ref>:<captured-old-oid>`. Never use an unqualified
    `--force`, wildcard mirror push, or tag overwrite. If any lease fails, stop
    and restart the before-state audit.
12. **Verify the server from scratch.** Re-read `ls-remote`, make a fresh clone,
    repeat fsck, same-method packing, hash tests, submission build, figure/deck
    checks, and confirm GitHub's default branch/tags/releases. Server size may
    lag and is not the acceptance criterion.
13. **Coordinate recovery and collaborators.** Publish the old/new commit map,
    require fresh clones (not merges from old history), rebase/recreate open
    PRs and WG worktrees, and document how to recover any removed object from
    the private mirror without pushing the unreduced history back accidentally.
14. **Retain rollback capability.** Keep the maintenance freeze until fresh-clone
    validation passes. Preserve both recovery mirrors and the captured ref map
    for an agreed retention period. Do not run destructive GC on any source or
    backup as part of the rewrite; server reclamation can occur later under a
    separate, approved operation.

## Audit conclusion

Phase 1 reduced the intended enrichment history by about 495.7 MB under a
consistent packing method. Phase 2 has credible, reviewable candidates, but the
largest low-risk opportunities are generated proof/render history and obsolete
WG reachability—not canonical biology inputs. The CHM13 GFF3, retained V7 layer,
manuscript figures, and direct submission inputs are explicitly protected. Any
later cleanup must be a new authorized operation following the complete safety
protocol above.
