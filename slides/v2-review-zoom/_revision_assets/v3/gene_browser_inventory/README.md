# Gene Browser Inventory for Review Zoom v3

Task: `review-zoom-v3-gene-browser-inventory`

Scope: find the interval, community, gene, and repeat annotation inputs needed to
replace weak slide 14a-style summaries with focused genome-browser-style PHR
views. This is documentation and data inventory only; no deck source or rendered
slide assets are changed here.

## Deliverables

- `README.md`: source decisions, target-locus table, and renderer guidance.
- `target_loci.tsv`: machine-readable candidate browser-view loci.
- `track_schema.tsv`: minimal tracks the downstream renderer should draw.
- `source_inventory.tsv`: source paths, coordinate systems, and stale-path
  caveats.

## Coordinate Systems

Use these two coordinate systems explicitly. Do not mix them without conversion.

1. **CHM13 v2.0 chromosome coordinates, BED 0-based half-open.**
   - Primary for the first rendered panels because `chm13.phrs.bed` and
     `phrs.genes.gff3` give simple chromosome intervals and gene labels.
   - `phrs.genes.gff3` and `chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz` are GFF3,
     therefore their feature coordinates are 1-based closed. Convert them to
     BED if the renderer expects half-open intervals.
2. **HPRCv2 1Mb trimmed PanSN sequence coordinates.**
   - Primary for cohort/haplotype rows and for acrocentric p-arm rows that are
     present in the HPRCv2 projected table but not all present in the repo-root
     CHM13 BED.
   - In `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv`,
     `seq` names encode the sample/haplotype/chromosome/original interval, and
     `region_start`/`region_end` are 0-based half-open offsets within that
     sequence. Draw the PHR from `region_start` to `region_end`, not necessarily
     the full 1Mb sequence window.

## Canonical Source Decisions

Canonical for v3 browser panels:

- PHR intervals:
  - `chm13.phrs.bed`
  - `CHM13-HG002.sub-telo-phrs.bed`
  - `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv`
- Community labels:
  - `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden.communities.tsv`
- Gene coordinates:
  - `phrs.genes.gff3` for CHM13 PHR gene rows.
  - `chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz` if transcripts/exons outside the
    PHR-intersected gene-only file are needed.
  - `/moosefs/guarracino/HPRCv2/PHR_III/hprc_annotations/*.gff3.gz` for
    haplotype-specific panels.
  - No `.gtf` or `.gtf.gz` files were found in the repo or the searched
    HPRCv2 PHR_III tree; the canonical gene-annotation inputs are GFF3.
- Repeat/TAR1 markers:
  - `/moosefs/guarracino/HPRCv2/PHR_III/hprc_repeatmasker/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed.gz`
  - `/moosefs/guarracino/HPRCv2/PHR_III/hprc_repeatmasker/*.RepeatMasker.bed.gz`
  - `/moosefs/guarracino/HPRCv2/PHR_III/annotations/subtelomeric_annotations.1Mb.rds`
- Gene-family and marker summaries:
  - `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_summary_table.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_enrichment.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_gene_families.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_enrichment_fisher.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/plots/d4z4_dux4l_by_community.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/plots/d4z4_perwindow_signal.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_per_gene.csv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/plots/il9r_distribution.csv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_by_arm.tsv`
  - `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/tar1_positional/tar1_positional_per_arm.tsv`

Noncanonical or stale for this purpose:

- `slides/v2-review-zoom/_typst/assets/s14_dux4.png`,
  `s14_or4f.png`, and `s14_tar1.png` are current deck crop caches, not data
  sources. The git provenance audit shows they are copied agent-878 crops with
  no committed crop geometry.
- `slides/v2/slide_14_gene_biology.R` is the best committed slide-level source
  for the old DUX4/OR4F/TAR1 summary, but it produces arm-level charts rather
  than browser tracks. Use it for provenance and sanity checks only.
- `paper_prep/_brainstorming/gene_copy_summary.csv` and
  `paper_prep/_brainstorming/EXECUTIVE_SUMMARY.md` are parked historical
  copy-number-aware enrichment outputs after commit `c50b99c`. They are useful
  context but should not override the HPRCv2 PHR_III tables.
- `/moosefs/guarracino/HPRCv2/PHR_III/flanking/similarity/hprcv2.1Mb.flanking.communities.tsv`
  is a 100kb flanking-community product with different meaning. Do not use its
  C1/C2 labels as the PHR arm-level C1-C15 labels.
- Sequence-level community files such as
  `hprcv2.1Mb.subtelo.seq-leiden-k50.assignments.tsv` reuse `C*` labels with a
  different partition. Use them only if the renderer explicitly shows
  sequence-level subtyping and relabels them.

## Target Locus Table

The target list below is intentionally focused. It covers the minimum requested
biologies and leaves one optional spare example if the renderer has space.

| Target | Coordinates | Community | Draw | Why it belongs |
|---|---|---|---|---|
| DUX4/D4Z4 C1 chr4q | CHM13 v2.0 BED: `chr4:193392741-193572740` | C1, chr4_qarm/chr10_qarm | PHR block, DUX4, DUX4L/DUX4-like array, DBET, FRG2, RPL23AP84, D4Z4/DUX4L label | Direct replacement for weak slide 14a. It shows the clinically familiar D4Z4/DUX4 locus as a discovered C1 PHR. |
| DUX4/D4Z4 C1 chr10q | CHM13 v2.0 BED: `chr10:134574995-134754994` | C1, chr4_qarm/chr10_qarm | PHR block, FRG2B, many DUX4L family copies, RPL23AP60, D4Z4/DUX4L label | Paired homologous chr10q view proves this is a systematic chr4q/chr10q community, not a one-off gene anecdote. |
| OR4F C3 chr3q | CHM13 v2.0 BED: `chr3:200846363-201101362` | C3, f7501/major duplicon block | PHR block, OR4F5, OR4F8BP, WASH8P, DDX11L8, FAM138D, SEPTIN/GTF2IP labels | OR4F-rich subtelomere with the largest OR4F count in the HPRCv2 OR4F summary: chr3_qarm has 846 total OR4F entries. |
| OR4F decay C8 chr15q | CHM13 v2.0 BED: `chr15:99625359-99750358` | C8 singleton | PHR block, OR4F28P, OR4F4, FAM138E, WASH3P, DDX11L9, pseudogene-fraction label | Shows the extreme OR4F decay endpoint: chr15_qarm has 854 OR4F annotations, 852 pseudogene annotations, 99.8% pseudogene fraction. |
| TAR1-rich C2 chr18p | CHM13 v2.0 BED: `chr18:2017-267017` | C2, chr10_parm/chr18_parm | PHR block, TAR1 RepeatMasker blocks, TUBB8B, IL9RP4, TAR1 copy-density label | Best positive TAR1 browser panel: chr18_parm is 100% TAR1-positive with mean 4.0 copies per sequence in the HPRCv2 table. |
| Acrocentric C7 p-arm group | HPRCv2 PanSN relative intervals: `CHM13#0#chr13:2544-502543_chr13_parm:0-500000`; `CHM13#0#chr14:2075-502074_chr14_parm:0-500000`; `CHM13#0#chr15:3258-503257_chr15_parm:0-500000`; `CHM13#0#chr21:2505-502504_chr21_parm:0-500000`; `CHM13#0#chr22:4138-504137_chr22_parm:0-500000` | C7, acrocentric p-arms | Small multiples with PHR blocks, TAR1 markers, MTCO-family pseudogene labels, rDNA-adjacent C7 band | Required acrocentric p-arm example. C7 is the rDNA-adjacent homogenization community and a strong community-level story. |
| IL9R/f7501 C3 chr9q | CHM13 v2.0 BED: `chr9:150390025-150615024` | C3, f7501/major duplicon block | PHR block, IL9RP1, SEPTIN14P22/P13, GTF2IP10, RPL23AP47, C3/f7501 label | Top gene-enrichment/community example. It ties C3 to IL9R pseudogene biology and the variable/AFR-enriched f7501 arms. |
| DDX11L/WASH C5 chr12p | CHM13 v2.0 BED: `chr12:2783-72783` | C5, RPL23A/WASH/DDX11L duplicon block | PHR block, DDX11L8, WASH8P, FAM138D, partial IQSEC3, optional TAR1 | Optional second gene-enrichment/community panel. C5 is the clean DDX11L/WASH subtelomeric duplicon example. |

The complete version with source-file paths and renderer notes is in
`target_loci.tsv`.

## Minimal Track Schema

The downstream renderer should draw at least these tracks:

| Track | Required | What to draw | Source |
|---|---|---|---|
| PHR intervals | yes | Filled region block for each selected PHR | `chm13.phrs.bed`; `CHM13-HG002.sub-telo-phrs.bed`; HPRCv2 `all-vs-all.1Mb.p95.id95.len.tsv` |
| Projected homologous regions | yes for HPRCv2 panels | Relative `region_start` to `region_end` block inside each PanSN sequence, with involved arms as tooltip/label | `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv` |
| Community label | yes | Colored band and text label, using arm-level C1-C15 | `hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`; `hprcv2.1Mb.subtelo.arm-leiden.communities.tsv` |
| Gene models | yes | Directional gene glyphs, compact labels, biotype styling | `phrs.genes.gff3`; full CHM13 Liftoff GFF3 if transcript detail is needed |
| Gene-family labels | yes | Highlight DUX4L, OR4F, IL9R/IL9RP, RPL23A, SEPTIN14, DDX11L, WASH, MTCO families | HPRCv2 community gene enrichment/family tables plus OR4F/IL9R summaries |
| Repeat markers | yes | Separate repeat lane, especially TAR1 blocks; D4Z4 via DUX4L array proxy unless a D4Z4 BED is provided | CHM13/HPRC RepeatMasker BEDs and `subtelomeric_annotations.1Mb.rds` |
| Marker summaries | optional | Small side labels: DUX4L copies, OR4F pseudogene fraction, TAR1 prevalence/copy count | `d4z4_dux4l_by_community.tsv`; `or4f_pseudogene_fraction.csv`; `community_tar1_by_arm.tsv`; `tar1_positional_per_arm.tsv` |
| Copy/community context | optional | Subtitle or legend with top family and statistical caveat | `community_summary_table.tsv`; `community_enrichment_fisher.tsv` |

The machine-readable version is in `track_schema.tsv`.

## Renderer Guidance

1. Start with CHM13 v2.0 single-locus panels for DUX4/D4Z4, OR4F, TAR1, IL9R,
   and C5 examples. They have straightforward chromosome coordinates and local
   `phrs.genes.gff3` gene rows.
2. Use a small-multiple layout for C7 acrocentric p-arms. Use HPRCv2
   `all-vs-all.1Mb.p95.id95.len.tsv` for chr13p because the repo-root CHM13 PHR
   BED omits that row while the HPRCv2 projected interval table includes it.
3. Keep repeat/TAR1 blocks visually separate from gene models. TAR1 is a repeat,
   not a gene.
4. Put community labels in a narrow color band above or below the PHR block.
   Use the arm-level C1-C15 assignments only.
5. Use gene-family badges for recognizable biology rather than dense labels for
   every gene. The most useful badges here are DUX4L/D4Z4, OR4F, IL9R/IL9RP,
   DDX11L/WASH/RPL23A/SEPTIN14, TAR1, and MTCO.
6. Preserve statistical caveats in small text outside the browser track. For
   example, OR-family community rows are qualitative support in the reviewed
   HPRCv2 Fisher table, not BH-significant enrichment.

## Validation Notes

- Source paths are recorded in `source_inventory.tsv` and summarized above.
- The target-locus table names the coordinate system for every coordinate row.
- `target_loci.tsv` specifies what each downstream browser panel should draw:
  interval blocks, gene labels, repeat markers, community bands, and optional
  numeric callouts.
- Git history and prior notes were used to distinguish canonical data products
  from stale crop assets and parked enrichment intermediates:
  - `slides/v2-review-zoom/_revision_assets/git_provenance/README.md`
  - `slides/v2-review-zoom/_revision_assets/14_gene_background/README.md`
  - `slides/v2-review-zoom/_revision_assets/14_gene_enrichment_or4f/README.md`
  - `git log -- chm13.phrs.bed phrs.genes.gff3 chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz subtelomeric_analysis_report.md paper_prep/figures/ed3 paper_prep/figures/ed4 paper_prep/surveys/SURVEY_03_gene_enrichment.md slides/v2/slide_14_gene_biology.R`
