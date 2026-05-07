# Hi-C and Single-Cell 3D Plot Candidates

Task: `review-zoom-v7-hic-3d-plots`

This directory stages v7 candidates for Erik's question: are there 3D plots we can show for the Hi-C / single-cell Hi-C work?

Short answer: yes, but there are two different meanings of "3D" here.

1. **Hi-C MDS plots** are 3D embeddings of a contact matrix. They are useful contact-space summaries, but they are not physical single-cell genome structures.
2. **Dip-C / sperm single-cell plots** come from reconstructed 3D particle coordinates (`.3dg.gz`). The best slide-ready versions are radial/community summaries, not raw whole-genome 3D particle clouds, because the radial plots explain the validation result without visual clutter.

No heavy contact extraction, `.mcool` querying, hickit modelling, or Dip-C remapping was run for this task. The staged assets are copied from existing renders and existing v6 PNG conversions.

## Staged Candidates

### `pngs/chm13_hic_mds_3d_coords.png`

Source: `/moosefs/guarracino/HPRCv2/PHR_III/HiC/CHM13/plots/MDS_3d_coords.png`

This is the cleanest existing Hi-C 3D MDS view. It places CHM13 chromosome arms/regions in a 3D MDS embedding computed from inter-chromosomal contact frequencies. Red marks the NOR/acrocentric-related entries from the analyzer inputs.

Use it only with explicit caption language:

> 3D MDS embedding of Hi-C contact frequencies; not a physical single-cell reconstruction.

Why it is useful: it directly answers "do we have a 3D plot of the Hi-C stuff?" without requiring new extraction or rendering. CHM13 is more readable than HG002 because it is haploid/reference-like and has fewer overlapping labels.

Why it needs care: the axes are abstract MDS coordinates. It should not be described as actual nuclear positions.

### `pngs/gm12878_dipc_radial_community.png`

Sources:

- PDF: `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_radial_community.pdf`
- PNG conversion copied from v6: `slides/v2-review-zoom/_revision_assets/v6/dipc_validation/pdf_pngs/gm12878_radial_community.png`
- Source tables: `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/gm12878_radial_community.tsv`, `gm12878_summary.tsv`, `gm12878_per_cell.tsv`
- Raw 3D coordinates: `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/3dg/*.3dg.gz`

This is the strongest single-cell 3D candidate. It summarizes actual GM12878 Dip-C 3D coordinates in two ways:

- left: mean normalized radial position by sequence community;
- right: within-community arm pairs have more similar radial positions than between-community arm pairs.

Plain-language caption option:

> In GM12878 single-cell Dip-C structures, same-community subtelomeric arms are closer in radial nuclear position than different-community arms.

The report-backed headline values are W/B = 0.931, Fisher p = 2.4e-05, and Mantel rho = 0.296.

### `pngs/sperm_all20_radial_community.png`

Sources:

- PDF: `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_radial_community.pdf`
- PNG conversion copied from v6: `slides/v2-review-zoom/_revision_assets/v6/dipc_validation/pdf_pngs/sperm_all20_radial_community.png`
- Source tables: `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_radial_community.tsv`, `sperm_all20_summary.tsv`, `sperm_all20_per_cell.tsv`
- Raw 3D coordinates: `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/output_x/3dg/*.3dg.gz` and `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/output_y/3dg/*.3dg.gz`

This is the best sperm single-cell candidate. It is the same radial/community plot type as the GM12878 Dip-C panel, but in 20 haploid sperm cells.

Plain-language caption option:

> The same community-proximity pattern is stronger in sperm single-cell 3D structures, despite the distinct compact sperm nucleus.

The report-backed headline values are W/B = 0.401, Fisher p = 3.9e-51, and Mantel rho = 0.202.

## Survey Notes

### Human Hi-C

Existing MDS PNGs were found for CHM13 and HG002:

- `/moosefs/guarracino/HPRCv2/PHR_III/HiC/CHM13/plots/MDS_3d_coords.png`
- `/moosefs/guarracino/HPRCv2/PHR_III/HiC/HG002/plots/MDS_3d_coords.png`

CHM13 is recommended because it is readable enough for a slide. HG002 is high-resolution but too cluttered to use directly because diploid arm labels overlap heavily.

The Randiak submission also has older HG002 MDS images:

- `/moosefs/guarracino/HPRCv2/submission_Randiak/images/mds_3d.png`
- `/moosefs/guarracino/HPRCv2/submission_Randiak/images/mds_3d_q.png`

These are useful provenance and document the original analysis command in `/moosefs/guarracino/HPRCv2/submission_Randiak/report.md`, but they are not recommended for the v7 deck because they are lower-resolution and superseded by the PHR_III renders.

I did not find a saved MDS coordinate TSV for the CHM13/HG002 renders. The source data appears to be the `.mcool` contact matrix plus the analyzer script; recomputing coordinates would require contact matrix extraction and was not done on the head node.

### Dip-C / GM12878

GM12878 has the strongest set of single-cell 3D assets:

- raw `.3dg.gz` coordinate files under `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/3dg/`;
- summary/radial/Mantel TSVs under `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/`;
- rendered radial and Mantel PDFs in the same directory.

The radial/community plot is recommended over a raw 3D particle cloud because it is already slide-readable and directly tied to the validation question.

### Sperm

Sperm has raw single-cell 3D coordinate files under:

- `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/output_x/3dg/`
- `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/output_y/3dg/`

The corrected enrichment output is under `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/`. The radial/community plot is recommended if the slide needs a second single-cell example after GM12878.

### PBMC

PBMC assets exist under `/moosefs/guarracino/HPRCv2/dipc_t2t/pbmc_hg19/enrichment_corrected/`, including radial and Mantel PDFs plus source TSVs. PBMC is not recommended for the v7 headline because the report says PBMC community-free analysis is unavailable, and the available PHR-specific community-based summary is weak/non-significant.

Inventory row to keep in mind:

- `/moosefs/guarracino/HPRCv2/dipc_t2t/pbmc_hg19/enrichment_corrected/pbmc_radial_community.pdf`

### RPE-1

RPE-1 has contact-derived Pore-C/CiFi validation outputs, including MDS-comparison PDFs and source contact matrices under:

- `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/RPE1/50000bp/`
- `/moosefs/guarracino/HPRCv2/PHR_III/RPE1_subtelo/`

I did not find RPE-1 single-cell `.3dg.gz` coordinate assets. Treat RPE-1 as contact-validation support, not as a single-cell 3D structure candidate for this request.

## Files

- `asset_inventory.tsv`: survey of found assets and source tables, with recommendation status.
- `source_manifest.tsv`: staged-file source paths, SHA256 hashes, and byte counts.
- `pngs/`: 3 slide-ready PNG candidates.
- `source_pdfs/`: source PDFs for the two Dip-C/sperm radial candidates.
- `SLIDE_PATCH.md`: recommended integration approach for the fan-in renderer.
