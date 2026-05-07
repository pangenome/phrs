# Gene Browser Panels for Review Zoom v3

Task: `review-zoom-v3-gene-browser-render-assets`

This directory contains five focused, slide-ready genome-browser-style panels
rendered from the upstream inventory task:

`slides/v2-review-zoom/_revision_assets/v3/gene_browser_inventory/`

The panels are candidates for the v3 deck fan-in. They intentionally cover a
small set of recognizable PHR biologies instead of dumping every subtelomeric
region.

## Generated Panels

| Panel | Files | What it shows |
| --- | --- | --- |
| `panel_01_dux4_d4z4_c1_chr4_chr10` | `.png`, `.pdf` | Paired C1 chr4q/chr10q PHR views. DUX4, DUX4L/D4Z4 proxy array, FRG/DBET, RPL23A pseudogene markers, and TAR1 repeats are shown with the same grammar in both rows. |
| `panel_02_or4f_c3_chr3q` | `.png`, `.pdf` | C3 chr3q OR4F-rich PHR with OR4F5/OR4F8BP and curated neighboring subtelomeric duplicon markers. |
| `panel_03_or4f_decay_c8_chr15q` | `.png`, `.pdf` | C8 chr15q OR4F endpoint where the HPRCv2 summary reports 854 OR4F annotations and 852 pseudogene annotations. |
| `panel_04_tar1_c2_chr18p` | `.png`, `.pdf` | C2 chr18p TAR1-rich PHR. TAR1 is drawn only in the repeat lane, separate from nearby TUBB8B and IL9RP4 gene models. |
| `panel_05_acrocentric_c7_p_arm_group` | `.png`, `.pdf` | C7 acrocentric p-arm small multiples using HPRCv2 PanSN relative offsets for chr13p, chr14p, chr15p, chr21p, and chr22p. |

Supporting files:

- `render_gene_browser_panels.R`: reproducible generator.
- `panel_manifest.tsv`: panel list, target IDs, row labels, coordinate systems,
  and track order.
- `render_track_schema.tsv`: fixed visual track schema used across all panels.
- `input_manifest.tsv`: exact input paths observed by the renderer.

## Track Grammar

Every panel uses the same top-to-bottom track order:

1. `community_band`: row-level arm community label with a fixed color.
2. `phr_interval`: dark blue block for the selected PHR interval.
3. `curated_gene_models`: directional gene glyphs with compact labels.
4. `repeat_or_proxy_markers`: TAR1 repeat blocks or DUX4L/D4Z4 proxy ticks.

The panels keep raw gene clutter out of the view. Gene rows are filtered to the
target-locus labels listed in `gene_browser_inventory/target_loci.tsv`; dense
DUX4L copies are collapsed visually as a DUX4L/D4Z4 proxy array.

Community colors are fixed within this directory:

| Community | Color role |
| --- | --- |
| C1 | DUX4/D4Z4 chr4q/chr10q band |
| C2 | TAR1-rich chr18p band |
| C3 | OR4F-rich chr3q band |
| C7 | Acrocentric p-arm group band |
| C8 | OR4F pseudogene endpoint band |

Gene-family colors are likewise fixed across panels: DUX4/DUX4L, OR4F, TAR1,
IL9R/IL9RP, TUBB8, MTCO, FRG/DBET, WASH/DDX/FAM, and
SEPTIN/GTF2I/RPL23A.

## Coordinate Conventions

- CHM13 single-locus panels use CHM13 v2.0 chromosome coordinates from
  `chm13.phrs.bed`. These intervals are BED 0-based half-open.
- Gene models come from `phrs.genes.gff3`. The script converts GFF3
  1-based closed coordinates to 0-based half-open intervals before plotting.
- TAR1 blocks come from the CHM13 RepeatMasker BED and are already 0-based
  half-open.
- The C7 acrocentric panel uses HPRCv2 `CHM13#0` PanSN rows from
  `all-vs-all.1Mb.p95.id95.len.tsv`. PHR intervals are drawn from
  `region_start` to `region_end` as offsets inside each 1 Mb trimmed sequence,
  not by assuming the entire sequence name is the PHR.
- The x-axis is always distance in kb from the panel window start. Row labels
  give the absolute source coordinates for each row.

## Inputs

The generator reads the following inputs directly:

| Input key | Path |
| --- | --- |
| `target_loci` | `slides/v2-review-zoom/_revision_assets/v3/gene_browser_inventory/target_loci.tsv` |
| `inventory_readme` | `slides/v2-review-zoom/_revision_assets/v3/gene_browser_inventory/README.md` |
| `inventory_track_schema` | `slides/v2-review-zoom/_revision_assets/v3/gene_browser_inventory/track_schema.tsv` |
| `chm13_phr_bed` | `chm13.phrs.bed` |
| `phr_gene_gff` | `phrs.genes.gff3` |
| `hprc_all_vs_all` | `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv` |
| `chm13_repeatmasker` | `/moosefs/guarracino/HPRCv2/PHR_III/hprc_repeatmasker/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed.gz` |
| `d4z4_dux4l_by_community` | `/moosefs/guarracino/HPRCv2/PHR_III/plots/d4z4_dux4l_by_community.tsv` |
| `or4f_pseudogene_fraction` | `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv` |
| `community_tar1_by_arm` | `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_by_arm.tsv` |

The same list is machine-readable in `input_manifest.tsv`.

## Reproduction

From the repository root:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v3/gene_browser_panels/render_gene_browser_panels.R
```

Optional environment overrides:

```bash
HPRC_ALL_VS_ALL_TSV=/path/to/all-vs-all.tsv \
CHM13_REPEATMASKER_BED=/path/to/chm13.repeatmasker.bed.gz \
D4Z4_DUX4L_BY_COMMUNITY_TSV=/path/to/d4z4.tsv \
OR4F_PSEUDOGENE_FRACTION_CSV=/path/to/or4f.csv \
COMMUNITY_TAR1_BY_ARM_TSV=/path/to/tar1.tsv \
Rscript slides/v2-review-zoom/_revision_assets/v3/gene_browser_panels/render_gene_browser_panels.R
```

The script writes only into this directory. It does not edit the Typst deck
source.

## Validation

- Ran `Rscript slides/v2-review-zoom/_revision_assets/v3/gene_browser_panels/render_gene_browser_panels.R`.
- Confirmed five PNG/PDF panel pairs were generated at 2933 x 1650 PNG
  resolution plus vector PDF.
- Reviewed the generated PNGs for slide-scale legibility and removed repeated
  lane labels that cluttered the multi-row views.
- Confirmed `panel_manifest.tsv` and `render_track_schema.tsv` record a
  consistent schema across panels.
- Confirmed no deck source files were edited.
