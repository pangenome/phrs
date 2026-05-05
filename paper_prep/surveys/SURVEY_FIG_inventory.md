# SURVEY — Figure inventory across repo

**Scope.** All image files (`.png`, `.pdf`, `.svg`, `.jpg`) under the project root, cataloged for Nature manuscript + 15-min talk planning.

**Search command.**
```
find . -type f \( -iname "*.png" -o -iname "*.pdf" -o -iname "*.svg" -o -iname "*.jpg" -o -iname "*.jpeg" -o -iname "*.gif" -o -iname "*.tiff" \) -not -path "*/.git/*" -not -path "*/.wg*"
```

**Totals.** 156 image files found. 0 SVG / 0 JPG / 0 GIF / 0 TIFF. All are PDF (134) or PNG (22).

**Distribution.**
| Location | Count | Class |
|---|---:|---|
| Repo root (`./`) | 91 | impg/identity coverage plots, Angela GSEA, rendered report PDF |
| `end-to-end-report/pedigree-plots/ceph1463-hifiasm/` | 42 | per-child untangle plots (21 children × 2 haplotypes) |
| `end-to-end-report/pedigree-plots/ceph1463-verkko/` | 20 | per-child untangle plots (10 children × 2 haplotypes) |
| `end-to-end-report/pedigree-plots/washu/` | 3 | per-child untangle plots (PAN027 ×2, PAN028 ×1) |
| **Total** | **156** | |

**Reading the LIKELY USE column.** `main` = Nature main figure candidate; `extended` = Extended Data; `SI` = Supplementary; `talk-only` = useful for the 15-min talk but too dense or off-topic for paper; `not-figure` = a rendered document, not a figure; `not-useful` = exploratory / superseded / draft only.

**Major caveat.** The narrative report at `end-to-end-report/report/*.md` does **not embed any figure paths** — no `.pdf` / `.png` / `.svg` is referenced inline. Most figures cited in the per-section surveys (`paper_prep/surveys/SURVEY_0[1-6]*.md`) live **outside this repo** under `/moosefs/guarracino/HPRCv2/...` and are NOT in this worktree. The gap analysis at the bottom enumerates which report sections currently have no in-repo figure.

---

## 1. Pipeline / inter-chromosomal sequence sharing — impg coverage plots

**Producer:** `plot-impg-coverage.R` (top-level R script).
**Topic mapping:** report §01 (Pipeline), §10 (Limitations — coverage caveats).
**Quality:** Vector PDF + matching high-DPI PNG (4800×2700 for `_wide` PNGs, 1440×1079 for genome-wide). Production-quality.

| PATH | TOPIC | CAPTION HINT | RES/QUALITY | LIKELY USE |
|---|---|---|---|---|
| `p_genome_wide_identity_heatmap.pdf` | Pipeline §01 | Genome-wide stacked identity heatmap (with inset zoom). Per-position max identity across 465 haplotypes per arm, all 24 chromosomes stacked. | Vector PDF, 2.5 MB, with inset draw_plot composition | **main** (top candidate Fig 1) |
| `p_genome_wide_identity_heatmap.png` | Pipeline §01 | Same as above, raster. | 1440×1079 PNG | talk-only / web preview |
| `p_genome_wide_identity_heatmap_no_inset.pdf` | Pipeline §01 | Genome-wide identity heatmap, no inset (clean version). | Vector PDF, 2.5 MB | **main** (alt Fig 1) |
| `p_genome_wide_identity_heatmap_no_inset.png` | Pipeline §01 | Same, raster. | PNG | talk-only |
| `p_genome_wide_numchrom_heatmap.pdf` | Pipeline §01 / §04 | Genome-wide heatmap of NUMBER OF CHROMOSOMES sharing each subtelomeric position (cross-arm count). | Vector PDF, 2.1 MB | **main** (Fig 1 panel B candidate) |
| `p_genome_wide_numchrom_heatmap.png` | Pipeline §01 / §04 | Same, raster. | 604 KB PNG | talk-only |
| `p_combined_alignments.pdf` | Pipeline §01 | Combined alignment metrics (identity / coverage) along subtelomeric position, with BED-region overlays for PHRs. | Vector PDF, 128 KB | extended |
| `p_combined_alignments_wide.pdf` | Pipeline §01 | Wide-format version of `p_combined_alignments`. | Vector PDF, 117 KB | extended (preferred) |
| `p_combined_alignments_wide.png` | Pipeline §01 | Same, high-DPI raster (4800×2700). | High-DPI PNG, 430 KB | talk-only |
| `p_combined_haplo_samples.pdf` | Pipeline §01 | Combined per-haplotype × per-sample alignment metrics with chrX-specific annotation line. | Vector PDF, 192 KB | extended / SI |
| `p_combined_haplo_samples_wide.pdf` | Pipeline §01 | Wide version. | Vector PDF, 193 KB | extended (preferred) |
| `p_combined_haplo_samples_wide.png` | Pipeline §01 | Same, raster. | 4800×2700 PNG, 383 KB | talk-only |
| `p_num_chromosomes.pdf` | Pipeline §01 / §04 | Number-of-chromosomes line plot vs subtelomeric position with PHR BED overlay. | Vector PDF, 104 KB | extended |
| `p_num_chromosomes_wide.pdf` | Pipeline §01 / §04 | Wide version. | Vector PDF, 91 KB | extended (preferred) |
| `p_num_chromosomes_wide.png` | Pipeline §01 / §04 | Same, raster. | 4800×2700 PNG, 345 KB | talk-only |

---

## 2. Pipeline / per-chromosome identity heatmaps

**Producer:** `plot_chromosome_identity_heatmap()` in `plot-impg-coverage.R`.
**Topic mapping:** report §01 (Pipeline) — per-chromosome detail; SI atlas.
**Pattern:** for each of 25 chromosomes (chr1–chr22, chrM, chrX, chrY) up to three views: full-chromosome, zoomed first 1 Mb (`zoom_0-1`), zoomed terminal 1 Mb (`zoom_last1mb`). chrM has no `last1mb`; chrY has no `zoom_0-1`.
**Quality:** Vector PDFs, 90–190 KB each. Production-quality, suitable for an SI atlas. 73 files total.

| PATH (family pattern) | TOPIC | CAPTION HINT | RES/QUALITY | LIKELY USE |
|---|---|---|---|---|
| `identity_heatmap_chr{1..22,M,X,Y}.pdf` (×25) | Pipeline §01 | Per-chromosome identity heatmap, full chromosome, all 465 haplotypes stacked. | Vector PDF | **SI atlas** (1 panel/chromosome) |
| `identity_heatmap_chr{1..22,M,X}.zoom_0-1.pdf` (×24) | Pipeline §01 | Same, zoomed to first 1 Mb (p-arm subtelomere). | Vector PDF | **SI atlas** (zoom set) |
| `identity_heatmap_chr{1..22,X,Y}.zoom_last1mb.pdf` (×24) | Pipeline §01 | Same, zoomed to terminal 1 Mb (q-arm subtelomere). | Vector PDF | **SI atlas** (zoom set) |

**Per-chromosome inventory (file presence matrix).** ✓ = file exists.

| Chr | full | zoom_0-1 | zoom_last1mb |
|---|:-:|:-:|:-:|
| chr1–chr22 | ✓ | ✓ | ✓ |
| chrM | ✓ | ✓ | — (no q-arm) |
| chrX | ✓ | ✓ | ✓ |
| chrY | ✓ | — | ✓ |

---

## 3. Gene enrichment — Angela's GSEA

**Producer:** Angela Gyamfi (Mefford lab, Memphis). Cited in `TODO.md` at lines 60–66.
**Topic mapping:** report §03 (Gene enrichment). Window = 1 Mb around PHR boundaries (caveat: ~10× wider than median PHR of 105 kb), so these capture the *neighborhood*, not strict PHR intervals. A PHR-only re-run is in TODO and would supersede these.

| PATH | TOPIC | CAPTION HINT | RES/QUALITY | LIKELY USE |
|---|---|---|---|---|
| `Figure1_GSEA_BP_vertical.pdf` | Gene enrichment §03 | GO Biological Process GSEA bar/dot plot, vertical layout — 1 Mb PHR-neighborhood window. | Vector PDF | **main** (Fig 3 candidate, pending PHR-only re-run) |
| `Figure_GSEA_MF_vertical.pdf` | Gene enrichment §03 | GO Molecular Function GSEA, vertical layout. | Vector PDF | extended (or main pair with BP) |

---

## 4. Pedigree recombination — untangle plots (WashU, T2T quality)

**Producer:** `scripts/pedigree/plot-pedigree-untangle.R` (per report §14 methods note).
**Topic mapping:** report §14 (Pedigree recombination).
**Quality:** Vector PDFs, 77–270 KB. Production-quality.
**Note.** WashU is the **primary** pedigree per §14 ("WashU is the primary dataset because of its T2T quality"). Plots show child haplotypes colored by parent-of-origin patches with a 50 kb error tolerance.

| PATH | TOPIC | CAPTION HINT | RES/QUALITY | LIKELY USE |
|---|---|---|---|---|
| `end-to-end-report/pedigree-plots/washu/PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf` | Pedigree §14 | PAN027 maternal hap1 untangled against PAN010 (mother). | Vector PDF, 99 KB | **main** (Fig 5 candidate, T2T quality) |
| `end-to-end-report/pedigree-plots/washu/PAN027.paternal_hap2_from_PAN011_father.untangle.pdf` | Pedigree §14 | PAN027 paternal hap2 untangled against PAN011 (father). | Vector PDF, 77 KB | **main** (Fig 5 partner panel) |
| `end-to-end-report/pedigree-plots/washu/PAN028.maternal_hap1_from_PAN027_mother.untangle.pdf` | Pedigree §14 | PAN028 maternal hap1 untangled against PAN027 (3-gen transmission). | Vector PDF, 270 KB | **main** (Fig 5: cross-generation transmission) |

---

## 5. Pedigree recombination — untangle plots (CEPH1463 hifiasm, supplementary)

**Producer:** same as §4.
**Topic mapping:** report §14 (Pedigree recombination, CEPH1463 cross-assembler validation).
**Quality:** Vector PDFs, ~98 KB each. Note the assemblies are **NOT T2T** (~780 fragmented contigs/sample), so per report §14 these are reported only where matched by both hifiasm AND verkko. 42 files.

| PATH (family pattern) | TOPIC | CAPTION HINT | RES/QUALITY | LIKELY USE |
|---|---|---|---|---|
| `end-to-end-report/pedigree-plots/ceph1463-hifiasm/{200081,200082,200084,200085,200086,200087,200101,200102,200103,200104,200106}.{maternal_hap2,paternal_hap1}_from_<parent>_<mother\|father>.untangle.pdf` (×22) | Pedigree §14 (hifiasm) | Per-child untangle plots, F1 generation 200xxx samples (Porubsky et al. 2025). | Vector PDF | SI |
| `end-to-end-report/pedigree-plots/ceph1463-hifiasm/NA{12877,12878,12879,12881,12882,12883,12884,12885,12886,12887}.{maternal_hap2,paternal_hap1}_from_<parent>_<mother\|father>.untangle.pdf` (×20) | Pedigree §14 (hifiasm) | Per-child untangle plots, NA128xx CEPH1463 generation. | Vector PDF | SI |

Total in directory: **42 PDFs** (21 children × 2 haps each). Naming convention: `<sample>.<maternal\|paternal>_<hap1\|hap2>_from_<parent>_<mother\|father>.untangle.pdf`.

---

## 6. Pedigree recombination — untangle plots (CEPH1463 verkko, supplementary cross-validation)

**Producer:** same as §4.
**Topic mapping:** report §14, used solely for cross-assembler validation against §5.
**Quality:** Vector PDFs. Smaller subset (10 children × 2 haps = 20 files; only the NA128xx generation, not 200xxx).
**Note.** Hap-numbering convention is FLIPPED relative to hifiasm: verkko uses `maternal_hap1` / `paternal_hap2`, hifiasm uses `maternal_hap2` / `paternal_hap1`. Cross-checking these requires care.

| PATH (family pattern) | TOPIC | CAPTION HINT | RES/QUALITY | LIKELY USE |
|---|---|---|---|---|
| `end-to-end-report/pedigree-plots/ceph1463-verkko/NA{12877,12878,12879,12881,12882,12883,12884,12885,12886,12887}.{maternal_hap1,paternal_hap2}_from_<parent>_<mother\|father>.untangle.pdf` (×20) | Pedigree §14 (verkko) | Per-child untangle, verkko assembly — used for cross-assembler agreement filter (only patches detected by both hifiasm + verkko are reported in §14). | Vector PDF | SI (cross-validation only) |

---

## 7. Rendered report (not a figure)

| PATH | TOPIC | CAPTION HINT | RES/QUALITY | LIKELY USE |
|---|---|---|---|---|
| `subtelomeric_analysis_report.pdf` | All sections | Rendered PDF of `subtelomeric_analysis_report.md` (the full report). | PDF v1.7, 1.9 MB | **not-figure** (document) |

---

## 8. Gap analysis — sections in the report TOC with NO figure in this repo

The `end-to-end-report/report/` TOC has 14 sections. The per-section surveys at `paper_prep/surveys/SURVEY_0[1-6]*.md` reference many figures, but most of those PDFs live **outside this repo** at `/moosefs/guarracino/HPRCv2/...`. Below is what exists in-repo for each TOC entry.

| Report section | In-repo figures? | Source survey | Status |
|---|---|---|---|
| §01 Pipeline | YES — impg coverage suite (§§1–2 of this inventory) + identity-heatmap atlas | `SURVEY_01_pipeline.md` | **partially covered** — but core community-detection figures (similarity heatmap, MDS, UMAP, UPGMA dendrogram, k-scan, Jaccard matrix) cited in survey live at `/moosefs/guarracino/HPRCv2/.../similarity/` and are NOT in this repo |
| §02 Annotation | **NONE** | `SURVEY_02_annotation.md` | **GAP** — TAR1 prevalence, internal (TTAGGG)n distribution, telomere-tract length plots not in repo |
| §03 Gene enrichment | YES — Angela's GSEA (§3 of this inventory) | (no survey file yet) | **partial** — covered for the 1 Mb window only; PHR-only re-run pending per `TODO.md` |
| §04 Within-arm heterogeneity | **NONE** | `SURVEY_04_heterogeneity.md` | **GAP** — all 12 `within_arm_heterogeneity_*.pdf` figures referenced live elsewhere; survey notes paths but they are NOT in this repo |
| §05 Hi-C / Pore-C validation | **NONE** | `SURVEY_05_hic_validation.md` | **GAP** — all `*_hic_*.pdf` (community heatmap, MDS comparison, bootstrap distributions, per-arm-pair scatter, Mantel scatter) live at `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/...` |
| §06 Dip-C / sperm validation | **NONE** | `SURVEY_06_dipc_validation.md` | **GAP** — all `gm12878_*.pdf`, `sperm_all20_*.pdf`, `pbmc_*.pdf` live at `/moosefs/guarracino/HPRCv2/dipc_t2t/...` |
| §07 Integrated interpretation | **NONE** | (no survey file yet) | **GAP** — likely needs a custom-built schematic figure (D4Z4-CTCF-lamin model, meiotic bouquet diagram); none drafted |
| §08 Mouse T2T | **NONE** | (no survey file yet) | **GAP** — mouse-pipeline community + Hi-C validation across 4 meiotic stages not in repo |
| §09 RPE-1 self-vs-self | **NONE** | (no survey file yet) | **GAP** — RPE-1 self-discovered communities + Hi-C validation not in repo |
| §10 Limitations | n/a | (no figures expected) | OK — text-only |
| §11 Summary of key findings | n/a | (no figures expected) | OK — text-only |
| §12 Literature and novelty | n/a | (no figures expected) | OK — text-only |
| §13 Appendix + References | n/a | (no figures expected) | OK — text-only |
| §14 Pedigree recombination | YES — WashU (3) + CEPH1463 hifiasm (42) + CEPH1463 verkko (20), §§4–6 of this inventory | (no survey file yet) | **covered** for primary panels; a summary "exchanges per generation" / "cross-assembler agreement matrix" figure is NOT yet drafted |

---

## 9. Synthesis recommendations for the figure architect

Based on what is in-repo and the gap analysis above:

**Top main-figure candidates already in-repo (no further generation needed):**
- `p_genome_wide_identity_heatmap.pdf` (or `_no_inset.pdf`) → Fig 1: the inter-chromosomal identity landscape.
- `p_genome_wide_numchrom_heatmap.pdf` → Fig 1 panel B: number of chromosomes sharing each position.
- `Figure1_GSEA_BP_vertical.pdf` (+ MF) → Fig 3 (gene enrichment), with the caveat that PHR-only re-run from `TODO.md` may supersede.
- `washu/PAN027.{maternal,paternal}_hap*.untangle.pdf` + `PAN028.maternal_hap1_from_PAN027_mother.untangle.pdf` → Fig 5: T2T pedigree recombination (3-generation chain).

**Critical gaps that block the manuscript** (must be generated or imported from `/moosefs/guarracino/HPRCv2/...`):
1. **Community-structure figure** — arm-level Leiden communities visualized (heatmap or network or MDS). This is the central "result" of the pipeline section and there is NO in-repo file. Source: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/` (per SURVEY_01).
2. **Hi-C validation panel** — community-blocked B/W ratio heatmap, Mantel scatter. Per SURVEY_05, all source PDFs live at `/moosefs/guarracino/HPRCv2/PHR_III/analysis/human/no_acrocentric/50000bp/`.
3. **Dip-C single-cell 3D figure** — Mantel scatter + radial-by-community for GM12878 and sperm. Per SURVEY_06, source `/moosefs/guarracino/HPRCv2/dipc_t2t/`.
4. **Within-arm heterogeneity panels** — at least the 4 highlighted in SURVEY_04 (separation overview, superpop composition, region length, TAR1 prevalence). Source not specified beyond "exist on disk".
5. **D4Z4-CTCF-lamin model schematic** (§07 integrated) — needs to be drawn, no source.
6. **Mouse meiotic Hi-C validation panel** (§08) — needs to be located or generated.
7. **RPE-1 self-discovered community panel** (§09) — needs to be located or generated.

**Talk-only assets:** all `_wide.png` rasters are well-suited for slides. The 73-PDF identity-heatmap atlas is too dense for a talk but works as an SI scrolling appendix.

**Import action:** before any figure assembly, mirror the relevant `/moosefs/guarracino/HPRCv2/...` PDFs into a `figures/` subdirectory of this repo so the manuscript build can vendor them.
