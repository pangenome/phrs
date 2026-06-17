# Slide 10m.2 clearer 3D visualization audit

Task: `review-zoom-v9-slide10m2-better-3d-viz`

This is a handoff package for the v9 fan-in render. It does not edit the final
Typst deck directly.

## Recommendation

For v9.1, show the same CHM13 PHR contact-space MDS in two steps:

- `chm13_phr_contact_mds_3d_view.png`
- `chm13_phr_contact_mds_3d_view.pdf`
- `best_replacement_chm13_phr_contact_mds.png`
- `best_replacement_chm13_phr_contact_mds.pdf`

The first image is a single oblique 3D view colored with the same sequence-community
palette as the projection slide. The second image shows the same D1/D2/D3 MDS as
two clearer 2D projections. Both are contact-space MDS summaries, not physical
genome reconstructions, but they are a better slide answer than the older
whole-arm/centromere/q-arm 3D scatter:

- It is CHM13, matching the current slide target.
- It uses the PHR/subtelomeric arm-region Hi-C contact-distance matrix from the
  main validation analysis, not the older whole p-arm/centromere/q-arm NOR plot.
- It renders a 3D MDS as two readable projections (`D1-D2` and `D1-D3`) with
  larger arm labels, stronger contrast, and sequence-community colors.
- It keeps the caveat explicit: this is bulk Hi-C contact-space MDS, not a
  physical single-cell 3DG reconstruction.

## Exact caption

Use this caption if the fan-in deck puts the note outside the figure:

> CHM13 Hi-C contact-space MDS over PHR/subtelomeric arm regions from
> `chm13_subtelomeric_regions.bed` at 50 kb resolution, shown as `D1-D2` and
> `D1-D3` projections of a 3D MDS. Colors mark sequence-defined subtelomeric
> communities; nearby points indicate similar bulk Hi-C contact profiles. This
> is contact-space MDS from bulk Hi-C, not a physical single-cell genome
> reconstruction.

The PNG/PDF also embeds this caption at the bottom so the caveat survives if the
image is used as a full-slide raster.

## Candidate inventory

See `candidate_inventory.tsv` for the structured audit. Summary:

| Candidate | Decision |
| --- | --- |
| v8 current CHM13 whole-arm contact-space MDS | Current source; interesting but visually hard to read, mostly whitespace, tiny labels, and not PHR-specific. Keep only as fallback. |
| v8 GM12878 whole-genome 3DG candidate | Clearer physical-coordinate context, but one GM12878 Dip-C cell only. It is not CHM13 Hi-C and not a replacement for slide `10m.2`. |
| HG002 `MDS_3d_coords.png` | Not better: diploid labels are denser and less readable. Still contact-space MDS, not physical reconstruction. |
| `submission_Randiak/images/mds_3d.png` | Not better: older HG002 whole-arm/centromere asset with dense labels. |
| `submission_Randiak/images/mds_3d_q.png` | Cleaner subset but q-arm-only HG002; too narrow for the CHM13 slide. |
| GM12878 and sperm radial-community overlays | Useful validation views already handled by slide 11a; radial statistic, not a 3D structure replacement. |
| CHM13 PHR contact-space MDS at 50 kb | Best replacement. Re-rendered here from local TSVs. |
| HG002 PHR contact-space MDS at 50 kb | Good backup, but less direct than the CHM13 replacement. |
| CHM13 flanking contact-space MDS | Good methods backup for unique flanks, but less direct than PHR/community slide `10m.2`. |
| Cross-organism overlay summaries | Useful robustness summaries, but not a slide-ready 3D visualization. |

## Why not a whole 3D genome structure plot?

A true whole-genome structure plot would need physical 3D coordinates from a
single-cell reconstruction. The current CHM13 source is bulk Hi-C, which measures
contact frequencies across many molecules/cells and supports a contact-space MDS
embedding, not a physical chromosome conformation. Presenting it as a whole 3D
genome structure would overclaim the source data.

There is an existing physical-coordinate candidate,
`gm12878_cell01_whole_genome_3dg_projection.png`, but it is one GM12878 Dip-C
cell. It is useful as optional context for the distinction between physical 3DG
coordinates and contact-space MDS, but it should not replace a CHM13 Hi-C slide.

## Source and reproducibility

Generated with:

```bash
Rscript slides/v2-review-zoom/_revision_assets/v9/slide10m2_better_3d_viz/make_slide10m2_better_3d_viz.R
```

Inputs:

- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_hic.dist_matrix.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_subtelomeric_regions.bed`
- `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_global_test.tsv`
- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_seq_global_test.tsv`

Derived files:

- `candidate_inventory.tsv`
- `best_replacement_chm13_phr_contact_mds.png`
- `best_replacement_chm13_phr_contact_mds.pdf`
- `best_replacement_mds_coords.tsv`
- `best_replacement_metrics.tsv`

Key metrics from `best_replacement_metrics.tsv`:

- `n_arms = 38`
- `D1/D2/D3 positive-eigenvalue variance = 8.0% / 6.4% / 5.9%`
- `Mantel rho = 0.656`, displayed as `p < 1e-4`
- Sequence-community contact enrichment `p = 3.46e-05`

## Validation

- Rendered with `Rscript make_slide10m2_better_3d_viz.R`.
- Confirmed output types with `file`: PNG is `1920 x 1080` RGB; PDF is version
  `1.5`.
- Visually inspected the PNG for contrast, label size, caption presence, and the
  explicit contact-space-vs-physical-reconstruction caveat.
