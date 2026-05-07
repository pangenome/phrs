# PGGB graph black layout asset for review zoom v6

Task: `review-zoom-v6-pggb-graph-black`

## Slide asset

- `pggb_graph_2d_black.png`
- 1920 x 1080 PNG, 8-bit colormap
- Shown status: main graph component, not the full graph.
- Component: ODGI layout TSV component `8`.
- Layout-node count in component `8`: 727,156.
- Transform: the existing ODGI X/Y layout coordinates are plotted as Y/X so the
  tall main component fits a 16:9 slide frame. No layout was recomputed.
- Palette change: the v5 blue graph marks were replaced with charcoal graph
  marks (`#111111`, alpha `0.30`) on a white background with a neutral border
  (`#b8c0cc`) for projector contrast.

This is a direct re-render of the prior v5 layout-node asset with a dark neutral
palette. It preserves the prior main-component/full-graph decision and continues
to plot layout nodes rather than extracting and redrawing graph edges from the
9.7 GB `.og` file.

## Source inputs

Primary graph and layout sources are unchanged from
`slides/v2-review-zoom/_revision_assets/v5/pggb_graph_odgi/`:

- ODGI graph:
  `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og`
- ODGI layout binary:
  `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay`
- ODGI layout TSV used for this render:
  `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.tsv`
- Prior slide asset used as the provenance baseline:
  `slides/v2-review-zoom/_revision_assets/v5/pggb_graph_odgi/pggb_graph_2d.png`

Previously inspected PNGs remain unsuitable as direct slide inputs:

- `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.draw.png`
  is 167 x 1000 and too sparse/narrow for direct slide use.
- `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.viz_multiqc.png`
  is 1968 x 157164 and too tall/raw for a clean 16:9 slide crop.

## Reproduction

Run from the repository root:

```bash
./slides/v2-review-zoom/_revision_assets/v6/pggb_graph_black/render_pggb_graph_black.sh
```

The wrapper runs:

```bash
Rscript render_pggb_layout_component8_black.R \
  --layout-tsv /moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.tsv \
  --component 8 \
  --out pggb_graph_2d_black.png \
  --render-log render_log.tsv \
  --point-color '#111111' \
  --point-alpha 0.30 \
  --point-cex 0.12 \
  --border-color '#b8c0cc' \
  --background-color white
```

No `odgi extract`, `odgi layout`, `odgi draw`, or `gfalook` command was
submitted for this v6 recolor. A SLURM job was not needed because the render
reads only the existing 41 MB layout TSV and does not read the full ODGI/GFA
graph. If a future edge-level main-component render is required, that extraction
and render should be submitted to the `tux` SLURM partition.

## Audit

`render_log.tsv` records the source paths, component identity, component bounds,
output dimensions, output byte count, palette parameters, render command, and
the fact that no SLURM job ID was used.
