# CHM13/hs1 PHR UCSC Browser Suite

Standalone UCSC Genome Browser slide suite for CHM13/hs1 PHR review.

Rendered outputs:

- `CHM13_PHR_UCSC_browser_suite.pdf`
- `CHM13_PHR_UCSC_browser_suite_by_community.pdf`
- `chm13_phr_ucsc_browser_suite.typ`
- `chm13_phr_ucsc_browser_suite_by_community.typ`
- `manifest.tsv`
- `manifest_by_community.tsv`
- `audit_37_vs_41.tsv`
- `_assets/ucsc/panels/*.png`
- `_assets/ucsc/html/*.html`
- `_assets/ucsc/sprites/*.png`

The deck contains 38 slides total: 1 explicitly labeled title/audit slide plus
37 main UCSC browser-image slides, one for each interval in `chm13.phrs.bed`.

The by-community deck contains 21 pages total: 1 title/audit page plus 20 pages
of full-width browser panels sorted by HPRCv2 Leiden community assignment. It
reuses the cached screenshots, adds visible clickable UCSC URLs for every panel,
and stacks panels only when they fit at the same full-width image scale.

## Rendering

Run from the repository root:

```bash
PATH=/home/erikg/micromamba/bin:/home/erikg/.local/bin:$PATH \
  python3 slides/chm13-phr-ucsc-browser/_scripts/render_ucsc_browser_suite.py --force
```

To regenerate only the community-sorted digest from the cached screenshots:

```bash
python3 slides/chm13-phr-ucsc-browser/_scripts/render_ucsc_browser_suite_by_community.py
cd slides/chm13-phr-ucsc-browser
typst compile chm13_phr_ucsc_browser_suite_by_community.typ CHM13_PHR_UCSC_browser_suite_by_community.pdf
```

The script uses:

- UCSC session `db=hub_3671779_hs1`
- UCSC session `hgsid=3966979908_lGks6rs34CqmdawD8iwY2YCYKVd3`
- `pix=1800`
- a browser-like curl user agent, `Mozilla/5.0`
- hs1 chromosome sizes from
  `https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes`

Each browser page is fetched as real UCSC `hgTracks` HTML. The script verifies
that the `chm13.phrs.bed` custom track is visible in the HTML, extracts the
`../trash/hgt/...png` and `../trash/hgtSide/...png` sprite paths, crops the
visible browser panel using the CSS offsets and dimensions, and appends the side
labels and data image into a single PNG.

## Window Rule

For each BED row, the script computes the terminal-anchored view requested for
review:

- inclusive PHR length: `end0 - start0 + 1`
- view length: `ceil(1.5 * inclusive PHR length)`
- arm inference: choose the smaller terminal gap, `start0` for p or
  `chrom_size - end0` for q
- p-arm browser window: `1..min(chrom_size, view_length)`
- q-arm browser window: `max(1, chrom_size - view_length + 1)..chrom_size`

The chr9q input row renders the expected browser window:

```text
chr9:150279748-150617247
```

## Manifest

`manifest.tsv` has one row per rendered browser panel and records:

- label, chromosome, and inferred arm
- PHR BED start/end and inclusive PHR length
- hs1 chromosome size
- browser start/end and terminal gap
- rendered image path and image dimensions
- cached HTML path
- UCSC custom-track confirmation
- UCSC URL

All main slides include the terminal chromosome end by construction: p-arm
windows start at coordinate 1, and q-arm windows end at the hs1 chromosome size.

## 37 vs 41 Audit

`chm13.phrs.bed` contains 37 intervals. The arm-level audit table
`paper_prep/figures/fig1/architecture_per_arm.tsv` contains 41 arms. This suite
therefore renders 37 browser panels, not 41.

`audit_37_vs_41.tsv` explicitly marks the four architecture-table arms that are
absent from the CHM13 BED and therefore have no UCSC panel here:

```text
chr13_p
chr6_p
chrY_p
chrY_q
```

There are no extra CHM13 BED arms relative to the architecture audit.

## Future Note

The within-population variation plot is not changed by this render suite. A
future review pass may want to use nearest same-superpopulation or
same-population PHR distance rather than centroid distance.
