# Slide patch: v9 slide 10m.2 clearer contact-space MDS

Do not edit the final deck in this task. This file is the handoff for
`review-zoom-v9-fanin-render`.

## Decision

Use `chm13_phr_contact_mds_3d_view.png` as an added slide immediately before
the clearer projection slide, then keep `best_replacement_chm13_phr_contact_mds.png`
as the follow-up slide.

This replaces the v8 whole p-arm/centromere/q-arm CHM13 MDS image with a clearer
CHM13 PHR/subtelomeric arm contact-space MDS rendered from the local 50 kb
validation matrix. It is not a physical single-cell reconstruction; the figure
and caption both state that explicitly.

Do not use the v8 GM12878 whole-genome 3DG projection as the primary replacement.
It is a true physical-coordinate projection, but it is one GM12878 Dip-C cell,
not CHM13 Hi-C and not a PHR/community validation summary.

## Exact caption

> CHM13 Hi-C contact-space MDS over PHR/subtelomeric arm regions from
> `chm13_subtelomeric_regions.bed` at 50 kb resolution, shown as `D1-D2` and
> `D1-D3` projections of a 3D MDS. Colors mark sequence-defined subtelomeric
> communities; nearby points indicate similar bulk Hi-C contact profiles. This
> is contact-space MDS from bulk Hi-C, not a physical single-cell genome
> reconstruction.

The replacement PNG already embeds this caption at the bottom. If the deck uses
an external caption block, preserve the wording above exactly and avoid
duplicating the embedded caption if it makes the page too dense.

## Preferred Typst patch

Replace the current `#captioned-figure-slide("10m.2", ...)` block in
`slides/v2-review-zoom/_typst/zoom_review_deck.typ` with two full-figure slides
so the rendered labels remain large:

```typst
#figure-slide(
  "10m.2",
  "CHM13 Hi-C contact-space MDS in 3D",
  "../_revision_assets/v9/slide10m2_better_3d_viz/chm13_phr_contact_mds_3d_view.png",
  source: "v9/slide10m2_better_3d_viz/chm13_phr_contact_mds_3d_view.png; same D1-D2-D3 contact-space MDS as next slide; colors are sequence communities; not physical 3DG",
)

#pagebreak()

#figure-slide(
  "10m.3",
  "The same 3D MDS is clearer as 2D projections",
  "../_revision_assets/v9/slide10m2_better_3d_viz/best_replacement_chm13_phr_contact_mds.png",
  source: "v9/slide10m2_better_3d_viz; source matrix: chm13_hic.dist_matrix.tsv; regions: chm13_subtelomeric_regions.bed; color: sequence communities; contact-space MDS, not physical 3DG",
)
```

Rationale: the PNG is already a 16:9 slide-ready figure with title, metrics,
legend, and caveat. Using `captioned-figure-slide` would make the plot too small
and duplicate the caption.

## Alternate patch if an external note is required

If fan-in wants the caveat in the deck note style instead of relying on the
embedded caption, use this shorter caption block and consider regenerating a
plot-only image later. Do not use this with the current PNG unless the duplicate
caption is acceptable.

```typst
#captioned-figure-slide(
  "10m.2",
  "Hi-C MDS gives a 3D contact-space view",
  "CHM13 PHR Hi-C contact-space MDS",
  "../_revision_assets/v9/slide10m2_better_3d_viz/best_replacement_chm13_phr_contact_mds.png",
  [
    CHM13 Hi-C contact-space MDS over PHR/subtelomeric arm regions from `chm13_subtelomeric_regions.bed` at 50 kb resolution, shown as `D1-D2` and `D1-D3` projections of a 3D MDS. Colors mark sequence-defined subtelomeric communities; nearby points indicate similar bulk Hi-C contact profiles. This is contact-space MDS from bulk Hi-C, not a physical single-cell genome reconstruction.
  ],
  source: "v9/slide10m2_better_3d_viz; source matrix: chm13_hic.dist_matrix.tsv; regions: chm13_subtelomeric_regions.bed; color: sequence communities",
)
```

## Files for fan-in

- `best_replacement_chm13_phr_contact_mds.png`
- `best_replacement_chm13_phr_contact_mds.pdf`
- `chm13_phr_contact_mds_3d_view.png`
- `chm13_phr_contact_mds_3d_view.pdf`
- `candidate_inventory.tsv`
- `best_replacement_metrics.tsv`
- `best_replacement_mds_coords.tsv`
- `README.md`
