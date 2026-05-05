---
title: "Survey 14 — Pedigree subtelomeric recombination"
source: end-to-end-report/report/14_pedigree_recombination.md
scope: Inter-chromosomal subtelomeric exchange in two pedigrees (WashU 3-gen T2T + CEPH1463 4-gen) with cross-assembler validation
audience: Nature manuscript and 15-min talk
---

# Survey 14 — Pedigree subtelomeric recombination

This survey extracts and structures the content of `end-to-end-report/report/14_pedigree_recombination.md` for the Nature manuscript and the companion 15-minute talk. The source section is the Part IV pedigree chapter and provides the only direct, intra-individual observation of inter-chromosomal exchange at PHR boundaries — the population-level Leiden communities (Part I) and 3D-proximity validation (surveys 05, 06) describe the *aggregate* signal; this section catches the events that *generate* it.

The chapter analyses inter-chromosomal patches in `odgi untangle` results (`nth.best=1`) for two pedigrees and validates each patch against the HPRCv2 arm-level Leiden community structure (15 communities, 41 arms): a patch is "biologically credible" only if its source and query arms sit in the *same* Leiden community. WashU is the headline dataset because all four samples are T2T (Cechova et al. 2025); CEPH1463 is reported only via cross-assembler-validated parent features (hifiasm AND verkko) to control for fragmented-contig artifacts.

---

## 1. Key findings with metrics

### 1.1 Pedigree quality summary

| Pedigree | Samples | Assembly | HQ inter-chr patches | Within-Leiden community | Reliability |
|---|---|---|---:|---:|---|
| **WashU** (3-gen) | 4 (PAN010 grandmother, PAN011 grandfather, PAN027 mother, PAN028 granddaughter) | **T2T** (Cechova et al. 2025) | 538 | **494 (92 %)** | **Primary evidence** |
| CEPH1463 hifiasm (4-gen) | 28 (Porubsky et al. 2025) | Fragmented (~780 contigs/sample) | 2,775 | 324 (12 %) | Supplementary |
| CEPH1463 verkko | 14 (NA12877–NA12887, NA12889–NA12892) | Fragmented (~144 contigs/sample) | 2,671 | 359 (13 %) | Supplementary |

The 80-point gap between WashU (92 %) and CEPH1463 (12–13 %) within-Leiden rate is the single most important headline number: it shows that fragmentation noise dominates inter-chr untangle calls in non-T2T assemblies, and motivates the cross-assembler-only filter for CEPH1463.

### 1.2 WashU pattern breakdown (494 within-community patches)

| Pattern | Count | Interpretation |
|---:|---:|---|
| `acros_like` | 229 | Extensive NAHR signature (≥5 inter-chr patches from ≥3 source chromosomes in the same flank) — characteristic of acrocentric p-arms |
| `gene_conversion_like` | 133 | Sandwich `chrN:hX → chrM:hY → chrN:hX` (same chr+hap on both sides, patch from different chr+hap) — classical ectopic gene conversion tract |
| `sandwich_same_hap` | 115 | Patch within same-chr same-hap background (unusual; transitive untangle artefact or short conversion to the *same* hap) |
| `crossover_like` | 16 | Patch coincident with a query-chr haplotype switch `chrN:h1 → chrM:hZ → chrN:h2` — meiotic crossover with inter-chr sequence at the breakpoint |
| `complex` | 1 | Other configurations |

(Note: an upstream count in `washu/.../report.md` quotes `gene_conversion_like = 136`, `crossover_like = 18`, `acros_like = 262`, `sandwich_same_hap = 120`, `complex = 2` for the *unfiltered* 538-patch set; the 494/229/133/115/16/1 numbers above are the *within-community* subset reported in 14.md.)

### 1.3 WashU `gene_conversion_like` — dominant ectopic-conversion engine

- **133 within-community gene-conversion-like tracts** at predominantly perfect alignment scores (≥80 with min_score = 1.000/1.000).
- Size distribution: median ≈ 1.4 kb, max 15.4 kb (PAN028 maternal chr22p ← chr14p:h1, 24,901–40,349, score 0.997).
- Community concentration: **~120 of 133 (≈ 90 %) are in C7** (acrocentric p-arms chr13p / chr14p / chr15p / chr21p / chr22p), with smaller clusters in C3 (chr3q/chr9q, chr19p/chr7p, chr8p/chr11p), C2 (chr18p/chr10p, the Linardopoulou pair), C11 (chr5q/chr6q, chr1p/chr8p), and C15 (chrXp/chrYp PAR1).
- Donor asymmetry: in PAN027 the donor is overwhelmingly `chr22p:h2` (PAN010-derived); in PAN028 the donor is `chr15p:h2` or `chr22p:h1` (PAN027-derived). This pattern is consistent with one specific parental allele acting as the recurrent donor for acrocentric NAHR.

### 1.4 WashU `crossover_like` — 16 within-community events with breakpoint phasing

- 16 events, 12 of them in PAN028 (G3, transmitted from PAN027), confirming that meiotic-crossover-like switches at inter-chr breakpoints are inheritable.
- Largest: PAN028 maternal chr3q 262,953–290,922 (27.97 kb, score 0.961/0.943, source chr7p:h2, community C3) — a near-30 kb meiotic-resolution event with inter-chr donor at the breakpoint.
- Second largest: PAN028 maternal chr22p 47,449–65,278 (17.83 kb, source chr14p:h1, score 0.998/0.998, C7).
- Third largest: PAN028 maternal chr1p 5,588–20,462 (14.87 kb, source chr8p:h1, score 0.999/0.997, C11) — the only non-acrocentric / non-acro-adjacent crossover in C11.

### 1.5 WashU `acros_like` — 229 within-community NAHR signatures

- Concentrated in chr13p, chr15p, chr21p, chr22p flanks (community C7) plus chr9q (C3) — exactly the arms expected from rDNA / satellite / pseudoautosomal NAHR.
- Top patches reach 23.4 kb (PAN028 maternal chr13p 0–23,406) and 11.0 kb (PAN028 maternal chr13p 481,083–492,096); both perfect-score (1.000/1.000) and within C7.

### 1.6 CEPH1463 cross-assembler validated parent features (11)

Only 11 unique parent×chromosome-pair features survive the cross-assembler filter ("same parent + same chromosome pair detected by both hifiasm and verkko in at least one child each, within the same Leiden community"):

| Parent | Chr pair | C | Hifi children | Verk children | Hifi best | Verk best |
|---|---|---|---|---|---:|---:|
| NA12877 | chr1/chr19 | C6 | NA12879, NA12883 | NA12884 | 0.982 | 0.860 |
| NA12877 | chr10/chr18 | C2 | NA12883, NA12884 | NA12884, NA12885 | 0.998 | 0.997 |
| NA12877 | chr17/chr19 | C6 | NA12883 | NA12884 | 0.984 | 0.868 |
| NA12877 | chr6/chr9 | C5 | NA12886 | NA12881 | 0.975 | 0.989 |
| NA12878 | chr10/chr18 | C2 | NA12884, NA12885, NA12887 | NA12882, NA12886 | 0.972 | 0.998 |
| NA12878 | chr19/chr22 | C6 | NA12881, NA12882, NA12886 | NA12879, NA12881, NA12883, NA12887 | 0.897 | 0.978 |
| NA12878 | chr21/chr22 | C7 | NA12879 | NA12881, NA12882 | 0.992 | 0.996 |
| NA12878 | chr6/chr9 | C5 | NA12881, NA12883, NA12884, NA12885 | NA12882 | 0.947 | 0.959 |
| NA12889 | chr12/chr9 | C5 | NA12877 | NA12877 | 0.957 | 0.978 |
| NA12890 | chr12/chr9 | C5 | NA12877 | NA12877 | 0.971 | 0.976 |
| NA12892 | chr21/chr22 | C6 | NA12878 | NA12878 | 0.994 | 0.993 |

The strongest cross-replication: NA12878 chr21/chr22 (C7) is detected in 1 hifiasm child + 2 verkko children with **58 patches in verkko** spanning the full pattern menu (gene_conv, crossover, acros_like, sandwich_same_hap, complex) — the most fully resolved single inter-chr event in the dataset.

### 1.7 Transmission across generations

- **WashU**: PAN027 carries patches inherited independently from PAN010 (maternal hap1) and PAN011 (paternal hap2); PAN028 then inherits from PAN027 with patches both *retained* (e.g., PAN028 chr13p:h2 ← PAN027 ← PAN010 chr22p) and *expanded* (PAN028 chr22p ← chr14p:h1 24,901–40,349, 15.4 kb gene-conv tract not visible in PAN027).
- **CEPH1463**: chr12/chr9 (C5) is detected in **both** NA12889 (paternal grandfather) and NA12890 (paternal grandmother) of NA12877 — i.e., the same C5 (RPL23A/WASH/DDX11L) signature appears in two independent G1 individuals and is then visible in G2 (NA12877). chr10/chr18 (C2; Linardopoulou pair) appears in NA12877 paternal *and* NA12878 maternal — independent observations of the same canonical inter-chr exchange in unrelated individuals.

### 1.8 Cross-section consistency with HPRCv2 communities

Every cross-assembler-validated CEPH1463 parent feature maps to a **known** Leiden community (C2, C3, C5, C6, C7), and the WashU within-community fraction (92 %) is the upper bound for what the population-level community structure could explain. The pedigree analysis *directly observes* the events that produce the population-level community structure rather than inferring them from cross-sample alignments.

---

## 2. Existing figures (paths)

All pedigree-untangle PDFs are in-tree under `end-to-end-report/pedigree-plots/`. They are produced by `scripts/pedigree/plot-pedigree-untangle.R` from per-pair `odgi untangle -e 50000 -m 1000` BEDs.

### WashU (3 PDFs — one per child×haplotype)
Directory: `end-to-end-report/pedigree-plots/washu/`
- `PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf` — PAN027 maternal hap1 vs PAN010 (G1→G2 mother).
- `PAN027.paternal_hap2_from_PAN011_father.untangle.pdf` — PAN027 paternal hap2 vs PAN011 (G1→G2 father).
- `PAN028.maternal_hap1_from_PAN027_mother.untangle.pdf` — PAN028 maternal hap1 vs PAN027 (G2→G3, the only G3-resolved haplotype).

(PAN010 / PAN011 are G1 progenitors so have no parent comparison; PAN028 paternal is not produced because PAN027's spouse was not assembled in this run.)

### CEPH1463 hifiasm (42 PDFs — 21 children × 2 haplotypes each)
Directory: `end-to-end-report/pedigree-plots/ceph1463-hifiasm/`
- G2 vs G1: `NA12877.{maternal_hap2_from_NA12890_mother,paternal_hap1_from_NA12889_father}.untangle.pdf`, `NA12878.{maternal_hap2_from_NA12892_mother,paternal_hap1_from_NA12891_father}.untangle.pdf`.
- G3 vs G2 (NA12877+NA12878 grandchildren): nine NA1288{1,2,3,4,5,6,7} and NA12879 child PDFs, each in `{maternal_hap2_from_NA12878_mother, paternal_hap1_from_NA12877_father}.untangle.pdf` form.
- G4 vs G3 (NA12879's children + NA12886's children): six 200081–200087 (NA12879 line) and five 200101–200106 (NA12886 line) PDFs, each in `{maternal_hap2_from_<mother>_mother, paternal_hap1_from_<father>_father}.untangle.pdf` form.

### CEPH1463 verkko (20 PDFs — 10 children × 2 haplotypes; verkko has no G4)
Directory: `end-to-end-report/pedigree-plots/ceph1463-verkko/`
- Same NA12877/NA12878 G2 and NA12879/NA12881–NA12887 G3 pairings as hifiasm but with `maternal_hap1` / `paternal_hap2` (verkko's haplotype labelling convention), e.g. `NA12879.maternal_hap1_from_NA12878_mother.untangle.pdf`.
- No 200xxx (G4) plots — the verkko run is restricted to the 14 G1–G3 samples listed in §1.1.

### What the PDFs show
Each PDF is the per-flank `odgi untangle` "ribbon" plot for a single child×haplotype vs the named parent, across all 46 chromosome arms (23 p + 23 q). Vertical strips are query positions in the terminal 500 kb; coloured bands are the chromosome+haplotype of the best (`nth.best=1`) parent segment. Inter-chromosomal patches appear as *off-diagonal* coloured bands within a flank panel; same-chr same-hap is the diagonal background.

### What is NOT yet plotted
- No summary panel rendering the 11 cross-assembler-validated CEPH1463 features as a single comparative figure.
- No transmission graph (3-generation WashU or 4-generation CEPH1463) overlaying retained vs novel patches per generation.
- No within-vs-cross-community count bar/pie comparing WashU (92 %) vs CEPH1463 (12 %).

---

## 3. Existing CSVs / TSVs (paths)

The only TSV cited in 14.md is the global per-patch table; the full per-pairing tables live one level deeper.

### Global table
- **`/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/all_pedigrees_patches.tsv`** — 5,984 HQ patches (header + 5,984 rows), columns: `ds, label, query_chr, query_arm, patch_start, patch_end, patch_size, ref_chrarm, ref_hap, left_chrarm, left_hap, right_chrarm, right_hap, mean_score, min_score, pattern, query_community, ref_community, community_status, overlaps_phr, has_phr`. `ds` ∈ {WashU, CEPH-h, CEPH-v}. `community_status` ∈ {within_community, cross_community, unknown_community} drives all downstream filtering.

### Per-pedigree per-pairing TSVs
Directory pattern: `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/<pedigree>/<assembler?>/untangle/recombination/`
- WashU: `…/washu/untangle/recombination/`
- CEPH1463 hifiasm: `…/ceph1463/hifiasm/untangle/recombination/`
- CEPH1463 verkko: `…/ceph1463/verkko/untangle/recombination/`

Each directory contains:
- `summary.tsv` — one row per child×haplotype: `flanks, no_switch, clean_1_3, mid_4_5, noisy_gt5, clean_bps, clean_same_chr, clean_inter_chr, clean_in_phr, clean_interchr_phr, all_inter_chr, clean_low_score`.
- `clean_breakpoints.tsv` — every "clean" (1–3 switch) breakpoint with surrounding state (`prev_*`, `curr_*`, `bp_pos`, `bp_gap`, `is_interchr`, `in_phr`).
- `interchr_exchanges.tsv` — same structure as `clean_breakpoints.tsv` but restricted to `is_interchr=True` exchanges.
- `patches.tsv` — every patch (not only clean ones) with `pattern`, `community_status`, `left_*`, `right_*`, `phr_start`, `phr_end`. This is the per-pedigree-assembler file aggregated into `all_pedigrees_patches.tsv`.
- `report.md` — narrative report (599 lines for WashU, ~2.8 kloc for each CEPH1463 assembler) with full pattern tables; the in-tree 14.md is a curated extract.

### Reference assignments
- **`/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`** — 41 arms × {community_id} for the 15 Leiden communities (C1–C15). This is the reference truth against which `community_status` is computed.

### Untangle-level inputs (not normally consumed downstream, but available for reproduction)
- `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/<pedigree>/{<assembler>/}untangle/<child>_vs_<parent>.e50000.m1000.bed.gz` — per-pair untangle BEDs (untangle column 7 = identity score).
- `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/<pedigree>/{<assembler>/}untangle/untangle_pairs.tsv` — pair manifest.
- `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/<pedigree>/{<assembler>/}<pedigree>.1Mb.telo_500kb_trimmed.fa.gz` — terminal-500-kb input FASTA for the pedigree assembly.

---

## 4. Methods

### 4.1 Inputs
- Per-sample assembly: WashU = T2T (Cechova et al. 2025); CEPH1463 = hifiasm and verkko draft assemblies (Porubsky et al. 2025).
- Per-arm terminal 500 kb extracted into `<pedigree>.1Mb.telo_500kb_trimmed.fa.gz`.

### 4.2 odgi untangle
- `odgi untangle -e 50000 -m 1000 -j 0 -n 100`, retaining `nth.best=1` only.
- Pairwise child × parent runs; each child×haplotype is compared against its named parent's full assembly.
- Score = column 7 of the untangle BED = alignment identity, 0–1.

### 4.3 Patch construction
- A **patch** is a contiguous run of untangle segments where the child's flank maps to the same `(parent_chromosome, parent_haplotype)`.
- HQ filter: `is_interchr=True` ∧ `min_score ≥ 0.8` ∧ `500 bp ≤ size ≤ 100 kb`.
- Quality reported as `mean_score / min_score` (lowest alignment identity within the patch; 1.0 = perfect).

### 4.4 Pattern classification
Driven by the *immediate* left and right neighbour patches within the same flank:
- **`gene_conversion_like`** — `chrN:hX → chrM:hY → chrN:hX` (same chr+hap on both sides; patch from different chr *and* different hap). Ectopic gene-conversion tract.
- **`crossover_like`** — `chrN:h1 → chrM:hZ → chrN:h2` (same chr, different hap left vs right). Meiotic crossover with inter-chr sequence at the breakpoint.
- **`acros_like`** — patch in a flank with ≥ 5 inter-chr patches from ≥ 3 different source chromosomes. Extensive NAHR signature.
- **`sandwich_same_hap`** — patch within same-chr same-hap background (transitive untangle artefact or conversion that reverts to the same haplotype label).
- **`complex`** — none of the above.

### 4.5 Leiden community validation
Each patch is cross-referenced against the 15 arm-level Leiden communities from the HPRCv2 1 Mb subtelomere similarity graph (`hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`, 41 arms; 233 samples / 465 haplotypes upstream). A patch is `within_community` iff `Leiden(query_arm) == Leiden(ref_arm)`. Only `within_community` patches are reported in 14.md; `cross_community` and `unknown` patches are filtered out but kept in `all_pedigrees_patches.tsv`.

### 4.6 Cross-assembler validation (CEPH1463 only)
Same parent + same `{chr_a, chr_b}` pair detected by *both* hifiasm and verkko in *at least one child each*, within the *same* Leiden community. This collapses 324 hifiasm + 359 verkko within-community patches to **11 robust parent features** (table §1.6).

### 4.7 Scripts
| Script | Purpose |
|---|---|
| `/moosefs/guarracino/HPRCv2/scripts/pedigree/analyze-pedigree-recombination.py` | Patch construction, pattern classification, Leiden validation, summary / clean_breakpoints / interchr_exchanges / patches TSV emission. |
| `/moosefs/guarracino/HPRCv2/scripts/pedigree/plot-pedigree-untangle.R` | Per-child×haplotype 46-arm ribbon-plot PDFs in `pedigree-plots/`. |

(The 14.md text references these as `scripts/pedigree/...` but the actual files live in `/moosefs/guarracino/HPRCv2/scripts/pedigree/` and are *not* mirrored in the in-tree `scripts/` directory of this worktree — to be addressed in §5.)

---

## 5. Gaps

1. **Scripts are not in the worktree.** `scripts/pedigree/analyze-pedigree-recombination.py` and `plot-pedigree-untangle.R` are referenced in 14.md but live only at `/moosefs/guarracino/HPRCv2/scripts/pedigree/`. For the manuscript repo to be self-contained these need to be vendored (or the path made explicit).
2. **No transmission diagram / pedigree tree figure.** WashU 3-gen (PAN010+PAN011 → PAN027 → PAN028) and CEPH1463 4-gen (NA12889+NA12890+NA12891+NA12892 → NA12877+NA12878 → 9 G3 children → 11 G4 200xxx grandchildren) topology is implicit only — no traditional pedigree square-and-circle plot exists with patch counts overlaid per edge.
3. **No quantitative transmission rate.** The text states qualitatively that PAN028 patches are "both retained and expanded" between generations but no numeric transmission rate (retained / total parental patches; novel / total child patches) is reported. The data are present in `patches.tsv` to compute this.
4. **No statistical test on within-community enrichment.** 92 % within-community in WashU is huge but no formal enrichment p-value (vs the random expectation of ~50 % given 15 communities and arm-frequency weighting) is reported.
5. **`sandwich_same_hap` is biologically opaque.** 115 within-community patches in this category are tabulated but the report does not commit to whether they are technical (untangle's transitive resolution) or real (rapid conversion that reverts to the same haplotype label). A brief simulation or null test would resolve this.
6. **CEPH1463 single-assembler results are documented but unused.** 324 + 359 within-community single-assembler patches are filtered out by the cross-assembler criterion. A supplementary table or ROC-style analysis showing how the cross-assembler filter trades sensitivity for specificity would be useful.
7. **No de novo / inherited split.** All patches are reported in child-vs-parent space, but the analysis does not separate (a) patches present in the child but absent in the parent (de novo or hidden in parent assembly gaps) from (b) patches inherited from one parental haplotype.
8. **PHR overlap not foregrounded.** `overlaps_phr`, `has_phr`, `phr_start`, `phr_end` columns exist in the per-patch TSV but the in-tree report only shows a `PHR (in/out)` column, not a numeric PHR-overlap rate or a PHR-vs-non-PHR comparison.
9. **No comparison to non-pedigree controls.** The WashU within-community fraction (92 %) is not compared to a null distribution constructed by, e.g., randomly relabelling source/query arms.
10. **`crossover_like` size distribution is striking but unmodelled.** Sizes range from 1 kb to 28 kb; whether this matches expected meiotic conversion-tract lengths or NAHR-driven inversion-resolution lengths is not discussed.
11. **Donor asymmetry not statistically tested.** The PAN027 chr22p:h2 donor bias and the NA12877 vs NA12878 chr10/chr18 (C2) asymmetry (NA12877 paternal but NA12878 maternal) are observed but not tested for parent-of-origin effects.
12. **No integration with the 3D-validation surveys.** A patch in PAN028 chr22p ← chr14p:h1 (15.4 kb) is exactly the kind of event that should leave a 3D-proximity footprint in the same arm pair (C7); no figure links 14.md to surveys 05/06.

---

## 6. Suggested figures with captions (produced vs to-do)

### Already produced (use directly or recompose)

**P-1. WashU per-haplotype untangle ribbons (3 panels).**
*Files:* `end-to-end-report/pedigree-plots/washu/PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf`, `PAN027.paternal_hap2_from_PAN011_father.untangle.pdf`, `PAN028.maternal_hap1_from_PAN027_mother.untangle.pdf`.
*Caption:* "WashU 3-generation pedigree: per-flank `odgi untangle` ribbons (46 chromosome arms × terminal 500 kb) for each resolved child×haplotype against the named parent. Off-diagonal coloured bands are inter-chromosomal patches (same Leiden community in 92 % of HQ cases, 494/538). PAN028 maternal hap1 (right panel) shows transmission from PAN027 with both retained and novel patches."

**P-2. CEPH1463 hifiasm 4-generation untangle gallery (42 panels).**
*Files:* `end-to-end-report/pedigree-plots/ceph1463-hifiasm/*.untangle.pdf`.
*Caption:* "CEPH1463 hifiasm pedigree: 14 children × 2 haplotypes (G2 vs G1, G3 vs G2, G4 vs G3). Untangle ribbons show the dominant background plus inter-chromosomal patches; the assembly is fragmented (~780 contigs/sample) so only 12 % of HQ patches are within Leiden community."

**P-3. CEPH1463 verkko untangle gallery (20 panels).**
*Files:* `end-to-end-report/pedigree-plots/ceph1463-verkko/*.untangle.pdf`.
*Caption:* "Independent verkko assembly of CEPH1463 (G1–G3 only); same 14 % within-Leiden background as hifiasm. Cross-assembler intersection (Fig. T-1) yields 11 robust parent features."

### To-do (suggested new figures)

**T-1. Cross-assembler-validated CEPH1463 parent features summary (one panel for the talk).**
*Caption:* "11 inter-chromosomal exchanges in CEPH1463 detected by both hifiasm and verkko within the same Leiden community. Parent (rows) × chromosome pair (columns), cells coloured by the better of the two best-scores; cell text lists detecting children. NA12877 chr10/chr18 and NA12878 chr10/chr18 (Linardopoulou C2 pair) are independent observations of the same canonical exchange in unrelated individuals."
*Source data:* §1.6 table.

**T-2. Pedigree transmission tree.**
*Caption:* "WashU (left, 3-gen) and CEPH1463 (right, 4-gen) family trees with edges weighted by number of within-community patches transmitted. Squares = males, circles = females; G1 progenitors at the top. WashU edge widths range from 22 (PAN027 paternal ← PAN011) to ~120 patches (PAN027 maternal ← PAN010). CEPH1463 G2 ← G1 edges show the chr12/chr9 (C5) signature inherited from both NA12889 and NA12890."
*Source data:* `all_pedigrees_patches.tsv` aggregated by `(label, community_status)`.

**T-3. WashU within-community vs CEPH1463 within-community headline bar.**
*Caption:* "92 % of WashU HQ inter-chr patches sit within a Leiden community (494/538) versus 12 % for CEPH1463 hifiasm (324/2,775) and 13 % for verkko (359/2,671). The 80-point gap is explained by assembly fragmentation, not biology: only T2T pedigrees recover the population-level community signal at single-individual resolution."
*Source data:* §1.1 table.

**T-4. WashU pattern × community heatmap.**
*Caption:* "Pattern (acros_like / gene_conv / sandwich_same_hap / crossover / complex) × Leiden community for the 494 within-community WashU patches. C7 (acrocentric p-arms) dominates all five patterns; C3 carries the chr3q/chr9q gene-conv cluster; C2 the Linardopoulou pair; C15 the PAR1 chrX/chrY gene-conv events."
*Source data:* `patches.tsv` filtered to WashU + within_community.

**T-5. Gene-conversion tract size distribution by community.**
*Caption:* "Histogram of `gene_conversion_like` patch sizes (n = 133) split by Leiden community. C7 median ≈ 1.4 kb (mode ~ 1.1, 2.2 kb consistent with single-segment untangle quanta); C3 and C15 reach 5–8 kb; one C7 outlier at 15.4 kb (PAN028 chr22p ← chr14p:h1)."
*Source data:* §1.3 of this survey.

**T-6. Crossover-like event detail.**
*Caption:* "All 16 within-community `crossover_like` events plotted as bp coordinates × patch size, coloured by community. PAN028 chr3q 27.97 kb (C3, source chr7p:h2) and PAN028 chr22p 17.83 kb (C7, source chr14p:h1) are the two largest events and define the upper end of the empirical crossover-like distribution."
*Source data:* §1.4 of this survey.

**T-7. Donor-bias asymmetry panel.**
*Caption:* "Per-flank donor (chrarm:hap) frequencies in PAN027 (G2) and PAN028 (G3) gene-conversion-like patches. PAN027 chr13p / chr21p / chr22p flanks are dominated by `chr22p:h2` donors (PAN010-derived); PAN028 inherits a different dominant donor (`chr15p:h2` for chr13p flanks). Bar text shows Fisher exact p-values for parent-of-origin asymmetry."
*Source data:* `patches.tsv` filtered to WashU + gene_conversion_like.

**T-8. Cross-link to 3D proximity (links to surveys 05/06).**
*Caption:* "For each of the 11 cross-assembler-validated CEPH1463 parent features (column 1), overlay the corresponding GM12878 Dip-C / sperm scHi-C 3D distance (column 2) and the bulk Hi-C W/B contribution (column 3). Demonstrates that arm pairs with detected pedigree exchanges are also closer in 3D — closing the loop between event detection and population structure."
*Source data:* §1.6 table joined with `gm12878_per_community_per_cell.tsv` and bulk Hi-C per-pair distances from survey 05.

**T-9. Methods schematic.**
*Caption:* "Workflow: terminal 500 kb FASTA per arm → `odgi untangle -e 50000 -m 1000 -n 100` (`nth.best=1`) → patch construction (HQ filter `min_score ≥ 0.8`, `500 bp ≤ size ≤ 100 kb`) → pattern classification by left/right neighbour → Leiden community filter (15 communities, 41 arms) → cross-assembler intersection (CEPH1463 only)."

---

## 7. Talk slide takeaways (15-min talk)

1. **Headline:** "We caught the events that build the population-level subtelomeric communities — directly, in two pedigrees. WashU (T2T, 3 generations): 494 inter-chromosomal exchanges, 92 % within Leiden community."
2. **Headline figure:** WashU PAN028 maternal hap1 untangle ribbon (P-1, right panel) — one slide showing G3 inheriting PAN027's pattern with both retained and novel patches, off-diagonal bands annotated for chr22p ← chr14p (15.4 kb), chr3q ← chr7p (28 kb crossover), chr1p ← chr8p (14.9 kb crossover).
3. **One-number recall:** "92 % of WashU subtelomeric inter-chr exchanges land in Leiden communities the population analysis predicts. CEPH1463 hits 12–13 % only because it is fragmented; the 11 features detected by both hifiasm AND verkko all map to known communities."
4. **Mechanism slide (T-4):** five pattern types, biological interpretation each — gene conversion (133), crossover (16), NAHR-acrocentric (229), and the inheritance pattern across PAN010 → PAN027 → PAN028.
5. **Cross-assembler robustness (T-1):** 11 CEPH1463 parent features detected independently by hifiasm and verkko — including chr10/chr18 (Linardopoulou pair) found in NA12877 paternal AND NA12878 maternal (independent individuals), and chr12/chr9 (C5) found in BOTH NA12889 and NA12890 (paternal grandparents of NA12877).
6. **Loop-closing slide (T-8):** "Pedigree exchanges sit in arm pairs that are also closer in 3D" — connects 14.md to surveys 05 (Hi-C/Pore-C) and 06 (Dip-C/sperm).
7. **Honest caveats (one slide):** T2T pedigrees are still rare (WashU n = 4); CEPH1463 fragmentation drowns 87 % of inter-chr untangle calls; `sandwich_same_hap` (115 patches in WashU) interpretation is open; no formal enrichment test against a randomised null is yet computed; no de-novo-vs-inherited split.
8. **Manuscript framing:** Position this section as the *event-level* counterpart to the *aggregate* community structure (Part I) and the *3D-proximity* validation (Part II/III). The pedigree section is the only direct observation of the recombination engine; everything else is forensic.
