# Slide Patch Recommendation

Task: `review-zoom-v7-hic-3d-plots`

Do not make final deck changes in this task. This note is for `review-zoom-v7-fanin-render`.

## Recommendation

Add **one optional new slide** only if v7 needs a direct visual answer to Erik's question, "can we make 3D plots of the Hi-C stuff?" Otherwise, keep the existing v6 Dip-C/sperm validation slide and do not add more 3D material.

Preferred new slide placement: after slide `10m` ("How 3D contact validates sequence communities") and before slide `10a`.

Suggested title:

> Hi-C MDS gives a 3D contact-space view

Suggested visual:

```typst
#figure-slide(
  "10m.1",
  "Hi-C MDS gives a 3D contact-space view",
  "../_revision_assets/v7/hic_3d_plots/pngs/chm13_hic_mds_3d_coords.png",
  source: "v7/hic_3d_plots README; existing CHM13 Hi-C MDS render copied from PHR_III/HiC/CHM13/plots",
)
```

Required caption/caveat:

> 3D MDS embedding of Hi-C contact frequencies; not a physical single-cell genome reconstruction.

## If Slide Count Cannot Grow

If the deck needs to stay the same length, replace slide `11` rather than `11a`.

Current slide `11`:

```typst
#figure-slide(
  "11",
  "Single-cell 3D: within-community arms are closer than between-community arms",
  "../_revision_assets/v3/11_wb_labels/slide11_explicit_distance_labels_candidate.png",
  source: "v3/11_wb_labels/make_slide11_explicit_distance_labels.R; explicit within-community vs between-community distance labels",
)
```

Replacement option:

```typst
#figure-slide(
  "11",
  "Single-cell 3D: communities share radial position",
  "../_revision_assets/v7/hic_3d_plots/pngs/gm12878_dipc_radial_community.png",
  source: "v7/hic_3d_plots README; existing GM12878 Dip-C radial community PDF/TSV; no 3D rerun",
)
```

Suggested caption language:

> GM12878 Dip-C uses actual single-cell 3D coordinates. Same-community arms have more similar radial nuclear positions than between-community arms.

Keep slide `11a` if possible. It already shows the GM12878 and sperm radial panels side by side from v6 and is still the best concise single-cell 3D validation view.

## Do Not Use Directly

- `HG002/plots/MDS_3d_coords.png`: high-resolution but too label-cluttered for a slide.
- `submission_Randiak/images/mds_3d*.png`: useful provenance, but older and lower-resolution.
- PBMC radial/community PDFs: available, but not a headline result because PBMC community-free analysis is unavailable and the PHR-specific community signal is weak/non-significant.
- RPE-1 MDS-comparison PDFs: contact-derived support plots exist, but no RPE-1 single-cell 3D coordinate asset was found.

## Bottom Line

Use CHM13 MDS only to answer the Hi-C visualization question, and use GM12878/sperm radial plots for actual single-cell 3D structure evidence. Keep those two concepts explicitly separated in the slide text.
