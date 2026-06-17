# Source Data Citations For Current Review-Zoom Draft

This note records the provenance that should be cited for Hi-C, Pore-C, Dip-C, sperm single-cell Hi-C, and mouse meiotic Hi-C figures in `BoG_2026_review_zoom_v9.pdf`.

## CHM13 bulk Hi-C

- Local analysis source: `/moosefs/guarracino/HPRCv2/code_PHR-and-3D.md`, Human Hi-C/Pore-C section.
- Local download recipe: CHM13 Arima Hi-C from `s3://human-pangenomics/T2T/CHM13/arima/`.
- Paper citation: Nurk, S. et al. The complete sequence of a human genome. Science 376, 44-53 (2022). https://doi.org/10.1126/science.abj6987
- Slide 10m.3 uses a CHM13 Arima Hi-C-derived distance matrix: `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/community_based/50000bp/chm13_hic.dist_matrix.tsv`.

## HPRC sample bulk Hi-C

- Local analysis source: `/moosefs/guarracino/HPRCv2/code_PHR-and-3D.md`, Human Hi-C/Pore-C section.
- Local download recipe: `data_hic_pre_release.index.csv` from the HPRC intermediate assembly repository, then S3 paths listed in that table.
- Consortium citation: Liao, W.-W. et al. A draft human pangenome reference. Nature 617, 312-324 (2023). https://doi.org/10.1038/s41586-023-05896-x

## HG002 Pore-C

- Local analysis source: `/moosefs/guarracino/HPRCv2/code_PHR-and-3D.md`, HG002 Pore-C pipeline.
- Local download recipe: `s3://human-pangenomics/submissions/5b73fa0e-658a-4248-b2b8-cd16155bc157--UCSC_GIAB_R1041_nanopore/HG002_R1041_PoreC/Dorado_v4`.
- Data description in local notes: R10.4.1, Dorado v4 SUP basecalling, 4 flow cells available; analysis used flow cells 1 and 2 before contact extraction.
- Best current citation: HPRC/human-pangenomics data resource plus the GIAB/HG002 source context. No specific paper citation was identified in the local documentation for this exact Pore-C submission.

## HG002 CiFi

- Local analysis source: `/moosefs/guarracino/HPRCv2/code_PHR-and-3D.md`, HG002 CiFi pipeline.
- Local notes describe HG002 CiFi v3.0 as Dennis lab, unpublished, DpnII enzyme, PacBio HiFi.
- Treat as unpublished/collaborator or private resource unless a project-specific citation is supplied.

## RPE-1 Pore-C and CiFi

- Local analysis source: `/moosefs/guarracino/HPRCv2/code_PHR-and-3D.md`, RPE-1 Pore-C/CiFi section.
- Local notes describe collaborator-provided asynchronous and mitotic RPE-1 proximity ligation data.
- Assembly citation in local notes: RPE1v1.1, Nature Communications 2025, DOI `10.1038/s41467-025-62428-z`, BioProject `PRJNA1195024`.
- Treat the proximity-ligation data themselves as collaborator/unpublished unless a project-specific citation is supplied.

## GM12878 Dip-C

- Local analysis source: `/moosefs/guarracino/HPRCv2/code_PHR-and-3D.md`, Dip-C sections.
- Data: GM12878 Dip-C cells from GEO `GSE117109`; raw FASTQ from BioProject `PRJNA473369` for remapping to CHM13v2.
- Paper citation: Tan, L., Xing, D., Chang, C.-H., Li, H. and Xie, X. S. Three-dimensional genome structures of single diploid human cells. Science 361, 924-928 (2018). https://doi.org/10.1126/science.aat5641

## Human sperm single-cell Hi-C

- Local analysis source: `/moosefs/guarracino/HPRCv2/code_PHR-and-3D.md`, Sperm scHi-C section.
- Data: GEO `GSE277040`, SRA/BioProject `PRJNA1160276`; local analysis selected top 10 X-bearing and top 10 Y-bearing QC-passing human sperm cells by clean contact count.
- Paper citation: Xu, H., Chi, Y., Yin, C. et al. Three-dimensional genome structures of single mammalian sperm. Nature Communications 16, 3805 (2025). https://doi.org/10.1038/s41467-025-59055-z

## Mouse meiotic Hi-C

- Local analysis source: `/moosefs/guarracino/HPRCv2/code_PHR-and-3D.md` plus `slides/v2-review-zoom/_revision_assets/hic_methods/README.md`.
- Local data: Zuo 2021 mouse meiotic Hi-C at leptotene, zygotene, pachytene, and diplotene stages.
- Paper citation: Zuo, W. et al. Stage-resolved Hi-C analyses reveal meiotic chromosome organizational features influencing homolog alignment. Nature Communications 12, 5827 (2021). https://doi.org/10.1038/s41467-021-26033-0
