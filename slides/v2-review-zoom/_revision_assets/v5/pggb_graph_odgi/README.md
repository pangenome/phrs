# PGGB graph ODGI layout asset for review zoom v5

Task: `review-zoom-v5-pggb-gfalook-2d-render`

## Slide asset

- `pggb_graph_2d.png`
- 1920 x 1080 PNG, RGB
- Shown status: main graph component, not the full graph.
- Component: ODGI layout TSV component `8`.
- Layout-node count in component `8`: 727,156.
- Transform: the existing ODGI X/Y layout coordinates are plotted as Y/X so the
  tall main component fits a 16:9 slide frame. No layout was recomputed.

The asset is a slide-oriented rendering of the existing ODGI layout coordinates.
It plots layout nodes rather than extracting and redrawing graph edges from the
9.7 GB `.og` file. This avoided a heavy ODGI graph read on the head node while
still using the PGGB/ODGI layout artifact for the main component.

## Source inputs

Primary graph and layout sources:

- ODGI graph:
  `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og`
- ODGI layout binary:
  `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay`
- ODGI layout TSV used for this render:
  `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.tsv`

Existing PNGs inspected before making the new asset:

- `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.draw.png`
  is 167 x 1000 and too sparse/narrow for direct slide use.
- `/moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.viz_multiqc.png`
  is 1968 x 157164 and too tall/raw for a clean 16:9 slide crop.

## Reproduction

Run from the repository root:

```bash
./slides/v2-review-zoom/_revision_assets/v5/pggb_graph_odgi/render_pggb_odgi_draw.sh
```

The wrapper runs:

```bash
Rscript render_pggb_layout_component8.R \
  --layout-tsv /moosefs/guarracino/HPRCv2/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.6e0e250.11fba48.645f51d.smooth.final.og.lay.tsv \
  --component 8 \
  --out pggb_graph_2d.png \
  --render-log render_log.tsv
```

No `odgi extract`, `odgi layout`, `odgi draw`, or `gfalook` command was
submitted for this final asset. A SLURM job was not needed because the selected
render reads only the 41 MB layout TSV and does not read the full ODGI/GFA graph.
If a future edge-level main-component render is required, that extraction and
render should be submitted to the `tux` SLURM partition.

## Audit

`render_log.tsv` records the source paths, component identity, component bounds,
output dimensions, output byte count, render command, and the fact that no SLURM
job ID was used.

An SVG was not emitted for this asset. A vector rendering of 727,156 plotted
layout nodes would be large and less practical than the 1920 x 1080 raster used
by the deck.
