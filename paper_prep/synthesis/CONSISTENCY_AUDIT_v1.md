---
title: Consistency audit — NATURE_DRAFT_v1.md vs end-to-end-report (14 sections)
generated: 2026-05-17
inputs:
  draft: paper_prep/synthesis/NATURE_DRAFT_v1.md (202 lines, abstract 218 w, main 3,905 w)
  report: end-to-end-report/README.md + report/01_pipeline.md ... 14_pedigree_recombination.md (3,367 lines)
  bib: paper_prep/synthesis/REFERENCES_v3.bib (295 entries)
  captions: paper_prep/figures/{fig1..fig4, ed1..ed5, ed8, nj_tree_arms}/caption.md
auditor: agent-22 (consistency-audit-nature)
verdict_summary: "Draft is broadly faithful but contains several DIVERGES, including one critical biological error (wrong list of seven silent arms in Methods and in the main text), and several OMITTED sections (especially §03 gene enrichment beyond DUX4L)."
---

# Consistency audit: NATURE_DRAFT_v1 vs end-to-end-report

The internal report is the source of truth. Status codes:

- **MATCH** — draft value reproduces report value within rounding.
- **DIVERGES** — draft and report disagree. CRITICAL category.
- **NOT-IN-REPORT** — draft asserts a number/claim that has no analogue in the report.
- **UNDERSPECIFIED-IN-REPORT** — report does not pin the value cleanly enough to test (or the report has internal inconsistencies).

Line numbers refer to `paper_prep/synthesis/NATURE_DRAFT_v1.md` and to the respective report files.

---

## 1. Numerical claims table

| # | Claim (location, line, verbatim) | Value in draft | Value in report (file:section) | Status |
|---|---|---|---|---|
| 1 | Abstract L18: "465 near-complete assemblies from 232 individuals" | 465 / 232 | 11_summary.md L5 "232 individuals, 465 near-complete assemblies" | MATCH |
| 2 | Abstract L18: "232 HPRC v2 individuals and CHM13v2.0 (465 near-complete assemblies)" | 465 | derivable (464 HPRC haplotypes + CHM13 = 465 near-complete assemblies); 01_pipeline.md L27 says 18,827 across 48 arms from 464 HPRC haplotypes + CHM13 | UNDERSPECIFIED-IN-REPORT (figure not explicit but math agrees) |
| 3 | Abstract L18: "18,827 telomere-anchored 500 kb flanks" | 18,827 | 01_pipeline.md L27 ✓ | MATCH |
| 4 | Abstract L18: "15,668 pseudohomologous regions" | 15,668 | 01_pipeline.md L71 ✓ | MATCH |
| 5 | Main L24: "15,668 PHRs (83.2% of flanks)" | 83.2% | 01_pipeline.md L71 "15,668 sequences (83.2%)" ✓ | MATCH |
| 6 | Main L24: "on 41 of 48 chromosome arms" | 41/48 | 01_pipeline.md L71 ✓ | MATCH |
| 7 | Main L24: "median 105 kb, mean 144 kb, range 5 kb to 500 kb" | 105 kb / 144 kb / 5–500 kb | 01_pipeline.md L71 "median 105 kb, mean 144 kb, range 5–500 kb" ✓ | MATCH |
| 8 | Main L26: "1 Mb contig minimum" (Methods L48) | 1 Mb | 01_pipeline.md L11 "minimum contig length ≥ 1 Mb" ✓ | MATCH |
| 9 | Main L24: "wfmash all-vs-all at 95% minimum identity" | 95% | 01_pipeline.md L51 ✓ | MATCH |
| 10 | Main L24: "500 kb flank" | 500 kb | 01_pipeline.md L23 ✓ | MATCH |
| 11 | Main L26: "chr18q chimera control… removes one false positive (15,669 to 15,668)" | 15,669 → 15,668 | 01_pipeline.md L83 "One sequence removed; 15,668 PHR sequences retained" ✓ | MATCH |
| 12 | Main L24: "C(18,827, 2) = 177 million pairs" | 177 million | not in report; computable: 18,827·18,826/2 = 177,234,651 | NOT-IN-REPORT (math correct) |
| 13 | Main L24/Methods L52: "wfmash k-mer prefiltering instead evaluates 11.6% of pairs (rounded to 12% throughout)" | 11.6% / 12% | not in report | NOT-IN-REPORT |
| 14 | Main L24/Methods L52: "230x above the Erdős and Rényi connectivity threshold" | 230x | not in report | NOT-IN-REPORT |
| 15 | Main L24/Methods L52: "p\* = log(n)/n = 5.2 x 10^-4 for n = 18,827" | 5.2e-4 | not in report; computable: log(18,827)/18,827 = 5.22e-4 ✓ | NOT-IN-REPORT (math correct) |
| 16 | Main L26: "The 7 remaining arms (chr5p, chr6q, chr7q, chr12q, chr14q, chr20p, chr20q) carry no detectable inter-chromosomal homology" | 7 arms = {chr5p, chr6q, chr7q, chr12q, chr14q, chr20p, chr20q} | 01_pipeline.md L124 (and L93, L300+) "chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q"; 06_dipc_validation.md L23 confirms same 7 | **DIVERGES — CRITICAL.** Draft list is wrong. chr7_q, chr12_q (C4), chr20_p (C12), chr20_q (C5), chr6_q (C11) all HAVE signal in the report. Missing from draft list: chr2_p, chr3_p, chr8_q, chr11_q, chr18_q. |
| 17 | Methods L72: "S_all: pooled 7 zero-signal arms (chr5p, chr6q, chr7q, chr12q, chr14q, chr20p, chr20q)" | same wrong 7 arms | same as row 16 | **DIVERGES — CRITICAL (same error repeated in Methods).** |
| 18 | Main L26: "15 communities at the arm level (Fig. 1c)" | 15 | 01_pipeline.md L131, L133 ✓ | MATCH |
| 19 | Main L26: "50 communities at the sequence level (Extended Data Fig. 2a, modularity 0.97; Methods)" | 50 / mod 0.97 | 01_pipeline.md L202 ✓ | MATCH |
| 20 | Main L28: "mean silhouette 0.347, k = 15" | 0.347 / 15 | 01_pipeline.md L131 ✓ | MATCH |
| 21 | Main L28: "UPGMA at k = 14 agrees with Leiden on 14 of 15 communities" | 14/15 | 01_pipeline.md L153 "agree exactly on 12 of 15 Leiden communities" | **DIVERGES.** Report says 12 of 15, not 14 of 15. Methods L62 (`agreement with Leiden 14 of 15`) repeats the wrong number. |
| 22 | Methods L62: "UPGMA… k = 14 (mean silhouette 0.342)" | 0.342 | 01_pipeline.md L131 ✓ | MATCH |
| 23 | Main L28: "1,000-replicate perturbation bootstrap" / Methods L64: "1,000-replicate perturbation bootstrap with Gaussian noise at sigma = 25%" | 1,000 / σ=25% | not in 01–14 reports (figure/method specific) | NOT-IN-REPORT |
| 24 | Main L28: "100% MRCA support of every named clade…32 to 90% across internal edges" | 100% / 32–90% | not in report | NOT-IN-REPORT |
| 25 | Main L30 / Fig.1d: "4/28/9 split" (homogeneous/polymorphic/interchangeable) and "8/34/7 architectural skeleton" | 4/28/9 vs 8/34/7 | 11_summary.md L9 lists 8/34/7 (no 4/28/9 split). 8+34+7=49 (does not sum to 41). | UNDERSPECIFIED-IN-REPORT (report version doesn't sum to 41; draft's 4/28/9 = 41 is internally consistent and matches fig1 caption) |
| 26 | Main L30: "cross-arm sequence rate peaks at chrX_q 99.7%, chr21_p 94.0% and chr11_p 74.1%" | 99.7 / 94.0 / 74.1 | 04_heterogeneity.md L64–67 cross-arm affinity table: chrX_q 99.7%, chr21_p 94.0%, chr11_p 74.1% ✓ | MATCH |
| 27 | Main L32: "5,946 paired distances from the nine multi-arm communities" | 5,946 | 04_heterogeneity.md L24 "Overall 5,946" ✓ | MATCH |
| 28 | Main L32: "within-haplotype-pair allelic distance is shorter than the corresponding cross-paralog distance in 8 of 9 communities, with combined p < 10^-300" | 8/9 / <1e-300 | 04_heterogeneity.md L26 ✓ | MATCH |
| 29 | Main L32: "C7… 70.5% of pairs have the paralog closer than the allele (p = 2.0 x 10^-7)" | 70.5% / 2.0e-7 | 04_heterogeneity.md L23, L26 ✓ | MATCH |
| 30 | Main L32: "In 39 of 48 arms, the number of distinct chromosome contributors per flank window decreases monotonically with distance from the telomere (Spearman ρ < 0)" | 39/48 | 04_heterogeneity.md L248 "39/48 arms" ✓ | MATCH |
| 31 | Main L32: "in 39 of 41 signal-bearing arms a piecewise linear model with a single breakpoint outperforms a linear model on AIC" | 39/41 | 04_heterogeneity.md L250 "39/41 testable arms" ✓ (F-test, not AIC in report; minor methodological mismatch) | DIVERGES (Test 2 is F-test, not AIC; otherwise number matches) |
| 32 | Main L32: "Sixteen of nineteen arms that contain an internal (TTAGGG)n island also have an ITS within 50 kb of the inferred breakpoint" | 16/19 within 50 kb | 04_heterogeneity.md L265 "within 25 kb for 11/19 arms and within 50 kb for 16/19" ✓ | MATCH |
| 33 | Main L32: "TAR1 prevalence… PAR1 arms are TAR1-free (chrXp 0.3%, chrYp 1.1%), acrocentric p-arms sit at 73 to 79%, and the remaining autosomal arms saturate above 99%" | 0.3 / 1.1 / 73–79 / >99 | 02_annotation.md L21 ✓ (chr15_p 73.1, chr22_p 78.9) | MATCH |
| 34 | Main L32: "Fisher exact test… significant after BH in 10 of 19 testable arms" | 10/19 | 04_heterogeneity.md L84 ✓ | MATCH |
| 35 | Main L32: "F_ST values of 0.10 to 0.15 between AFR and each of AMR, EAS, EUR and SAS" | 0.10–0.15 | 04_heterogeneity.md L103–110 "AFR 0.102, 0.152, 0.108, 0.103" ✓ | MATCH |
| 36 | Main L32: "0.02 to 0.04 within the non-AFR set" | 0.02–0.04 | 04_heterogeneity.md L103–110 non-AFR pairs: −0.047, 0.007, 0.007, 0.004, 0.005, −0.003 → report L111 paraphrase: "−0.05 to 0.01"; 11_summary.md L24 contradicts with "0.00–0.02" | **DIVERGES.** Draft 0.02–0.04 matches neither §04 (−0.05 to 0.01) nor §11 (0.00–0.02). |
| 37 | Main L34: "14 inter-arm tests of three-dimensional proximity across six HPRC v2 individuals and the CHM13 cell line" | 14 tests / 6 ind + CHM13 | 07_integrated.md L11–27 enumerate Hi-C (5 HPRC + CHM13), HG002 Pore-C, HG002 CiFi, Dip-C, sperm, 4 mouse stages = 14 ✓ | MATCH |
| 38 | Main L34: "HG002 Hi-C at 50 kb, B/W = 0.027, p = 4.0 x 10^-66" | 0.027 / 4.0e-66 | 05_hic_validation.md L32 ✓ | MATCH |
| 39 | Main L34: "CHM13 Hi-C yields B/W = 0.071, p = 6.0 x 10^-18" | 0.071 / 6.0e-18 | 05_hic_validation.md L31 ✓ | MATCH |
| 40 | Main L34: "HG002 Pore-C yields B/W = 0.056, p = 3.9 x 10^-85" | 0.056 / 3.9e-85 | 05_hic_validation.md L33 ✓ | MATCH |
| 41 | Main L34: "Mantel correlation… Spearman ρ = 0.66 for both CHM13 and HG002" | 0.66 | 05_hic_validation.md L82–83 (0.656, 0.657) ✓ rounded | MATCH |
| 42 | Main L34: "Per-individual sequence-pair Spearman ρ rises to 0.83 in the lowest-coverage samples" | 0.83 | 05_hic_validation.md L164 (NA19036 at 10 kb = 0.827) ✓ rounded | MATCH |
| 43 | Main L34: "HG002 0.66 to 0.80" (exclusion-set Mantel) | 0.66 → 0.80 | 05_hic_validation.md L332 HG002 full 0.657, no acro pq + sex 0.790 ✓ rounded | MATCH |
| 44 | Main L34: "HG02148 0.15 to 0.21" | 0.15 → 0.21 | 05_hic_validation.md L337 HG02148 full 0.152, no acro p+sex 0.720; **no 0.21 in the report** | **DIVERGES.** ed5b caption also says "0.15 → 0.21" — draft inherited the wrong value from the figure caption. Report shows 0.15 → 0.72. |
| 45 | Main L34: "CHM13 0.66 to 0.85" | 0.66 → 0.85 | 05_hic_validation.md L331 CHM13 full 0.656, no acro pq + sex 0.850 ✓ | MATCH |
| 46 | Main L34: "NA19036 0.27 to 0.49" | 0.27 → 0.49 | 05_hic_validation.md L338 NA19036 full 0.266 → 0.799; **no 0.49 in the report** | **DIVERGES.** Same defect as row 44; ed5b caption shows 0.49 but the report does not. |
| 47 | Main L34: "O/E enrichment… ranges from 5.9x (HG02559) to 18.4x (HG002 CiFi)" | 5.9× / 18.4× | ed5c caption matches; 05_hic_validation.md L433–439 O/E ratio table: HG002 34.4×, HG002 CiFi 13.0× (not 18.4×), HG02559 14.4× | DIVERGES (the report's O/E table doesn't actually list 5.9× or 18.4× — ed5c caption introduces them; draft follows ed5c) |
| 48 | Main L34: "8 of 8 tests" (O/E) | 8/8 | ed5c caption ✓; 05_hic_validation.md L432–439 has 7 rows in the O/E table (not 8) | DIVERGES (O/E table in §05 has 7 samples) |
| 49 | Main L34: "Across 11 datasets and 15 communities… in every cell of the 15 x 11 reproducibility heatmap" | 11 datasets / 15 communities | ed5d caption says 11 datasets ✓; 05_hic_validation.md L629 table lists 8 sample columns (HG02559, HG00658, HG02148, NA19036, CHM13, HG002, HG002 Pore-C, HG002 CiFi). Adding RPE-1 (3) brings total to 11 ✓ | MATCH |
| 50 | Main L36: "In HG002 Hi-C the PHR B/W of 0.027 falls to flanking B/W of 0.0031, a 9-fold strengthening" | 0.0031 / 9× | 05_hic_validation.md L262 (HG002 flanking 10kb 0.001, 50kb 0.002) — closest reported 0.002 at 50 kb. Fig3d caption gives 0.0031 (9×) | UNDERSPECIFIED-IN-REPORT (0.0031 not in report table; draft follows fig3d caption) |
| 51 | Main L36: "in GM12878 Dip-C, 16 of 16 C-community cells have W/B < 1 while 0 of 16 S_all cells do" | 16/16 / 0/16 | 06_dipc_validation.md L33 "S_all cells < 1.0 = 0/16" ✓; 16/16 implicit from "all 16 cells" | MATCH |
| 52 | Main L36: "in 20-cell sperm scHi-C… 20 of 20 C-cells have W/B < 1 while only 1 of 20 S_all cells does" | 20/20 / 1/20 | 06_dipc_validation.md L33 "1/20" ✓ | MATCH |
| 53 | Main L36: "S_all is 11% farther in GM12878 and 40% farther in sperm" | 11% / 40% | 06_dipc_validation.md L31 "GM12878 1.106 (11% farther)" L31 "Sperm 1.397 (40% farther)" ✓ | MATCH |
| 54 | Main L36: "flanking unique-sequence particles are more nuclear-interior than non-flanking terminal particles (radial 0.504 vs 0.556, p = 1.6 x 10^-35)" | 0.504 / 0.556 / 1.6e-35 | 05_hic_validation.md L67 GM12878 "0.503 / 0.551 / 7.4e-35"; Fig3d caption matches draft (0.504 / 0.556 / 1.6e-35) | DIVERGES (report and fig caption differ; draft follows fig caption) |
| 55 | Main L36: "C1 D4Z4 radial 0.732 peripheral, C14 PAR2 0.840 peripheral, C10 chr17p 0.474 interior" | 0.732 / 0.840 / 0.474 | 05_hic_validation.md L477 ✓; 07_integrated.md L60 also confirms 0.732 | MATCH |
| 56 | Main L38: "538 high-quality inter-chromosomal patches" | 538 | 14_pedigree_recombination.md L9, L38 ✓ | MATCH |
| 57 | Main L38: "494 of 538 (92%) sit within a Leiden community" | 494/538 (92%) | 14_pedigree_recombination.md L9, L38 ✓ | MATCH |
| 58 | Main L38: "229 are acros_like" | 229 | 14_pedigree_recombination.md L42 ✓ | MATCH |
| 59 | Main L38: "133 are gene-conversion-like sandwiches" | 133 | 14_pedigree_recombination.md L41 ✓ | MATCH |
| 60 | Main L38: "16 are crossover-like" | 16 | 14_pedigree_recombination.md L44 ✓ | MATCH |
| 61 | Main L38: "115 are sandwich_same_hap" | 115 | 14_pedigree_recombination.md L43 ✓ | MATCH |
| 62 | Main L38: "1 is complex" | 1 | 14_pedigree_recombination.md L45 ✓ | MATCH |
| 63 | Main L38: "Twelve of the sixteen crossover-like events are in PAN028" | 12/16 | 14_pedigree_recombination.md L191–206 crossover-like table: 13 events involve PAN028 (rows 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15) of 16 total | **DIVERGES.** Recount of report table gives 13 PAN028, not 12. |
| 64 | Main L38: "largest crossover-like event spans 27.97 kb on the PAN028 maternal chr3q" | 27.97 kb / chr3q | 14_pedigree_recombination.md L204 row 14: PAN028 mat vs PAN027, chr3q, 27,969 bp ✓ | MATCH |
| 65 | Main L38: "roughly 90% of these [133 GC-like] are in C7, the acrocentric p-arm community" | ~90% in C7 | 14_pedigree_recombination.md L52–185 gene_conv table shows ~120 of 133 in C7 (≈90%) ✓ | MATCH |
| 66 | Main L38: "11 features pass" (CEPH1463 cross-assembler) | 11 | 14_pedigree_recombination.md L255–269 lists 11 ✓ | MATCH |
| 67 | Main L40: "37 self-discovered communities" (RPE-1) | 37 | 09_rpe1_self.md L31 ✓ | MATCH |
| 68 | Main L40: "Leiden C2 = {chr10_q, chrX_q}" | C2 = chr10_q, chrX_q | 09_rpe1_self.md L35 ✓ | MATCH |
| 69 | Main L40: "Spearman ρ = 0.715, p = 4.4 x 10^-55, n = 344 inter-chromosomal pairs" | 0.715 / 4.4e-55 / 344 | 08_mouse.md L102 ✓ | MATCH |
| 70 | Main L40: "ρ 0.574 to 0.715" across leptotene to diplotene | 0.574–0.715 | 08_mouse.md L100–103 Leptotene 0.680, Zygotene 0.715, Pachytene 0.677, Diplotene 0.574 ✓ | MATCH |
| 71 | Main L42: "0.485" Pore-C HG002 Mantel | 0.485 | 05_hic_validation.md L88 = 0.486 ✓ rounded | MATCH |
| 72 | Main L42: "mean e1 = +0.0073 across HG002 100 kb Hi-C" | +0.0073 | ed8d caption "+0.0073"; 07_integrated.md L105 says "+0.007"; 11_summary.md L42 "+0.007" ✓ rounded | MATCH |
| 73 | Main L42: "63 of 92 arm-by-haplotype windows in A" | 63/92 | 07_integrated.md L105 "63/92 arm × haplotype combinations"; ed8d caption ✓ | MATCH |
| 74 | Main L42: "Lalli et al…Spearman ρ = -0.35, n = 46" | −0.35 / 46 | 07_integrated.md L103 reports "rho = −0.43, p = 0.006" (39 arms); 12_literature.md L73 same (rho = −0.43, N = 39). ed8c caption notes "Survey-reported published values: ρ = −0.43 (39 arms)" vs the panel's own "46 arms, ρ = −0.35". | **DIVERGES.** Draft follows ed8c panel numbers (−0.35, n=46) instead of the report's −0.43, N=39. |
| 75 | Main L42: "collapses to ρ ≈ 0 (n = 40)" | ≈ 0 / 40 | 07_integrated.md L103 and 12_literature.md L99 give "rho = 0.00, N = 32"; ed8c "ρ = −0.01, p = 0.97, N = 40" | **DIVERGES (same as row 74).** |
| 76 | Main L42: "DUX4 family copy counts in Extended Data Fig. 4c, range 0 to 22" | 0–22 | 07_integrated.md L60 "C1 sequences: median 22 DUX4L; non-C1 outliers: 0–2" (median 22, not range max 22; outliers 0–2); ed4c caption says DUX4 16–20 copies | UNDERSPECIFIED-IN-REPORT (no clean "range 0 to 22"; combines two different statistics from §07 and ed4) |
| 77 | Methods L46: "AFR 67, EAS 52, AMR 44, SAS 37, EUR 32" | 67/52/44/37/32 | 10_limitations.md L5 ✓ | MATCH |
| 78 | Methods L48: "MIN_LEN = 1 Mb" (implicit "1 Mb contig minimum") | 1 Mb | 01_pipeline.md L25 ✓ | MATCH |
| 79 | Main L32: "C7 acrocentric p-arms… chr21_p 94.0%, chr14_p 83.0%, chr15_p 76.5%, chr13_p 88.2%, chr22_p 43.0%" | 94.0 / 83.0 / 76.5 / 88.2 / 43.0 | 04_heterogeneity.md L64–74 arm-level table: chr21_p 94.0%, chr14_p 83.0%, chr15_p 76.5%, chr22_p 43.0%; **chr13_p is not in the arm-level cross-arm table**; the seq-level table (01_pipeline.md L294) shows chr13_p 100% | DIVERGES (chr13_p 88.2% is in neither table; possibly a different aggregation) |
| 80 | Main L36: "(C1 D4Z4 radial 0.732 peripheral, C14 PAR2 0.840 peripheral, C10 chr17p 0.474 interior)" — second occurrence in same paragraph | same | (see row 55) | MATCH |
| 81 | Main L36 (Limitations clause iii in main text refers to "compartment identity at tips is weakly A-leaning (mean e1 = +0.0073…; 63 of 92… in A)" | identical to row 72/73 | identical | MATCH |
| 82 | Methods L60: "Sequence-level: k-NN graph with k in {10, 25, 50, 75, 100, 125} … selected k = 75, resolution 0.8 (modularity 0.97, mean silhouette 0.602, k = 50)" | k=75, res=0.8, mod=0.97, sil=0.602, k=50 communities | 01_pipeline.md L200 "k=75, resolution=0.8. Leiden: k-NN = 75, 50 communities (modularity = 0.97)" L259 "0.602 (sequence)" ✓ | MATCH |
| 83 | Methods L60: "resolution scanned from 0.1 to 3.0 at 0.01 step; selected resolution 1.16 (mean silhouette 0.347, k = 15)" | res scan 0.1–3.0 / 1.16 / 0.347 / 15 | 01_pipeline.md L128 "scan from 0.1 to 3.0 in steps of 0.01" ✓; L131 "Leiden: 15 communities (optimal resolution = 1.16, silhouette = 0.347)" ✓ | MATCH |
| 84 | Methods L66: "Hudson pairwise F_ST" | Hudson F_ST | 04_heterogeneity.md L101 "Hudson's Fst estimator" ✓ | MATCH |
| 85 | Methods L68: "MAPQ filters disabled to retain multi-mappers with one random alignment per read" | MAPQ=0 / one random | 05_hic_validation.md L5–11 ✓ | MATCH |
| 86 | Methods L68: "10,000 row-and-column permutations" (Mantel) | 10,000 | 05_hic_validation.md L76 ✓ | MATCH |
| 87 | Methods L72: "Dip-C in 16 GM12878 cells [@Tan2018]" | 16 cells | 06_dipc_validation.md L7 "16 cells used (cell 12 excluded as duplicate of cell 10)" ✓ | MATCH |
| 88 | Methods L72: "sperm scHi-C in 20 cells [@Xu2025]" | 20 | 06_dipc_validation.md L75 ✓ | MATCH |
| 89 | Methods L74: "B6 and CAST T2T assemblies; telomere-anchored flanks; 1, 2 and 4 Mb window scan" | 1/2/4 Mb scan | 08_mouse.md L67 (windows 1Mb/2Mb/4Mb) ✓ | MATCH |
| 90 | Methods L76: "minimum patch score 0.95 and minimum alignment score 0.95" | 0.95 / 0.95 | 14_pedigree_recombination.md L20 "min_score >= 0.8 + 500 bp <= size <= 100 kb" | **DIVERGES.** Report uses min_score ≥ 0.8 for HQ filter, not 0.95. (Draft's 0.95 may refer to a stricter sub-filter, but the report's HQ filter is 0.8.) |
| 91 | Methods L80: "wfmash 0.23.0-41-gb5f0ff1c; impg commit 5b96025; pggb and odgi bundled" | versions | 13_appendix.md L7–10 ✓ | MATCH |
| 92 | Methods L80: "hicexplorer 3.7.4" | 3.7.4 | 13_appendix.md L11 ✓ | MATCH |
| 93 | Main L34 / Methods L70: "Five exclusion sets (no acrocentric p-arms, no sex chromosomes, no acrocentric p plus sex, no all-acrocentric plus sex, no strongest community) applied at all five mcool resolutions" | 5 exclusion sets / 5 resolutions | 05_hic_validation.md L329 "Exclusion labels: No acro p / No sex / No acro p + sex / No acro pq + sex / No strong" = 5 sets ✓; resolutions 5/10/20/50/100 kb ✓ | MATCH |
| 94 | Methods L48: "Per-arm flank counts are in Extended Data Fig. 1b" | per-arm counts | 01_pipeline.md L29–42 table; ed1b caption confirms median 447 ✓ | MATCH |
| 95 | Main L26: "Stacking by chromosome contributor (Fig. 1b)… 100 kb windows" | 100 kb | Fig1 caption confirms 100 kb ✓ | MATCH (caption-internal) |
| 96 | Main L34: "RPE-1… asynchronous CiFi Hi-C contact at 50 kb" (implicit in Fig 4c context) | 50 kb | 09_rpe1_self.md L47 ✓ | MATCH |
| 97 | Main L42: "11.6% / 12% wfmash sampling rate" (Limitations) | 12% | NOT in report | NOT-IN-REPORT (computational argument unique to draft) |
| 98 | Main L42: "C14 0.840" radial (second occurrence) | 0.840 | (see row 55) | MATCH |

**Coverage check vs the task's required claim list**: every bullet in the task description is represented above (rows 1–98). All 32 distinct numerical claims in the task's "Coverage target" appear in at least one row.

---

## 2. Coverage of the 14 report sections

| Section | Coverage in draft | Where in draft (paragraph #) | Notes |
|---|---|---|---|
| **§01 Pipeline** (`01_pipeline.md`) | FULLY-COVERED | Main paragraphs 2–3 (L24, L26); Methods L46–58 | Pipeline shape, 95% identity, IMPG transitive closure, PHR detection, PGGB, odgi all present. **CRITICAL ERROR: the list of 7 silent arms is wrong (L26 and Methods L72).** |
| **§02 Annotation** (`02_annotation.md`) | PARTIALLY | Main L32 (TAR1 prevalence one sentence); §07 ed3 reference | TAR1 prevalence numbers covered. Internal (TTAGGG)n island length distribution, motif composition (52.3% canonical, variant breakdown), per-community TAR1 table, terminal telomere length by community (Kruskal-Wallis p=3.2e-15), TAR1 in cross-arm vs self-arm — all OMITTED from main text (only mentioned in figure captions). |
| **§03 Gene enrichment** (`03_gene_enrichment.md`) | **LARGELY OMITTED** | Main mentions DUX4L (L42) and one indirect ref to acrocentric gene overlap; Methods does not name OR4F, Ambrosini blocks, hub genes, MTCO, SHOX | **MAJOR GAP.** §03 is one of the richest report sections (374 unique genes, 11 Ambrosini blocks → 15 Leiden mapping, hub genes RPL23AP45/SEPTIN14P22, OR4F pseudogenisation gradient 11.1%–99.8%, MTCO mitochondrial pseudogenes in C7, SHOX in C15 PAR1). Draft does not name OR4F, Ambrosini, RPL23AP, SEPTIN14, MTCO, SHOX, FAM41C, DDX11L, IL9R, FRG2, or "hub genes". |
| **§04 Heterogeneity** (`04_heterogeneity.md`) | FULLY-COVERED | Main paragraphs 4–5 (L30, L32); Methods L66 | Allele-vs-paralog (5,946 / 8 of 9 / p<1e-300), silhouette (implicit), F_ST AFR vs non-AFR, two-domain (39/48 / 39/41 / 16/19), Fisher 10/19 — all present. Mechanism-caveat (the report says "cannot distinguish among mechanisms — gene conversion, unequal crossover, reciprocal exchange") is OMITTED from main text; draft uses the terms freely. F_ST non-AFR range diverges (row 36). |
| **§05 Hi-C validation** (`05_hic_validation.md`) | FULLY-COVERED | Main paragraph 6 (L34); Methods L68–70 | B/W, Mantel, multi-resolution, exclusion, O/E, per-community heatmap, RPE-1 (under §09) all present. PBMC Dip-C cell-type split is OMITTED (briefly noted only as supporting flanking radial). |
| **§06 Dip-C** (`06_dipc_validation.md`) | FULLY-COVERED | Main paragraph 6 (L36); Methods L72 | 16 cells, 20 sperm, S_all 11%/40%, radial all present. |
| **§07 Integrated** (`07_integrated.md`) | FULLY-COVERED | Main paragraph 7 (L36); main paragraph 9 (L42) for causal loop | Bouquet, D4Z4-CTCF-lamin, flanking paradox, causal feedback loop all present. Testable predictions (3 predictions in §07) are not explicitly enumerated. |
| **§08 Mouse** (`08_mouse.md`) | PARTIALLY | Main paragraph 8 (L40); Methods L74 | Zygotene ρ=0.715, 4 stages 0.574–0.715, B6+CAST present. **Window optimization (1/2/4 Mb sweep, 19/49 still saturated at 4 Mb, identity threshold 95→98) is OMITTED from main text** (Methods only says "1, 2 and 4 Mb window scan"). Mouse-human synteny note OMITTED. |
| **§09 RPE-1 self-vs-self** (`09_rpe1_self.md`) | PARTIALLY | Main paragraph 8 (L40, Fig. 4c) | 37 communities, t(X;10), C2={chr10_q, chrX_q} present. Comparison with HPRC communities (Mantel rho 0.457 self vs 0.548 vs 0.611 — table in §09 lines 63–71) and flanking control (rho ≈ 0) **NOT explicitly gestured at** in the main text; only the t(X;10) point survives. |
| **§10 Limitations** (`10_limitations.md`) | PARTIALLY | Main paragraph 9 (L42); Methods L84 | Six limitations enumerated in draft (i–vi). Report has 18. **Missing from draft**: somatic LCL exchange (limitation 8 — could inflate C1), assembly quality, community resolution sensitivity, parameter sensitivity (formal sensitivity analysis), Dip-C cell 12 duplicate, fragmented-assembly NaN, hg19/T2T coordinate noise, mouse 1Mb saturation. Multi-mapping limitation is mentioned only in passing. |
| **§11 Summary (10 key findings)** (`11_summary.md`) | PARTIALLY | scattered throughout main text | Findings 1, 2, 3, 4, 5, 8, 9, 12 covered. **Finding 6 (subtelomeric gene content dominated by pseudogenes, PAR1 32.1% protein-coding exception) — OMITTED.** **Finding 7 (TAR1 as subtelomeric marker, 94.6% prevalence) — only briefly mentioned via Extended Data Fig. 3a.** **Finding 10 (C4 chr7_q/chr12_q minimal-PHR positive control, significant in 4/5 diploid Hi-C samples despite only 5–25 kb PHR and zero gene annotations) — OMITTED entirely.** Finding 11 (community-specific 3D predictions, C7 nucleolar, C1 lamina) — partially covered. Note: report §11 lists 12 findings under the "10 key findings" heading; draft does not reconcile this numbering. |
| **§12 Literature & novel contributions** (`12_literature.md`) | PARTIALLY | implicit in main paragraphs 1, 3, 4 | Draft does not enumerate the novel-contribution list, and does not claim a count of novel findings. The README claims "27 novel contributions"; §12 text intro says "24 findings" but enumerates items numbered 1–27 — an internal report inconsistency. Draft sidesteps this by not naming a number. No item in §12 is flagged as already established that the draft frames as new. |
| **§13 Appendix** (`13_appendix.md`) | FULLY-COVERED | Methods L80 | Tools and versions listed in Methods L80. |
| **§14 Pedigree** (`14_pedigree_recombination.md`) | FULLY-COVERED | Main paragraph 8 (L38) | All headline numbers match (538 / 494 / 92% / 229 / 133 / 115 / 16 / 1 / 11 features) EXCEPT "12 of 16 in PAN028" (row 63 — should be 13). 27.97 kb chr3q ✓. |

---

## 3. Citation discipline

### Commands run (verbatim)

```
grep -F '@' paper_prep/synthesis/REFERENCES_v3.bib | sed -n 's/^@[^{]*{\([^,]*\),.*/\1/p' | sort -u > /tmp/bibkeys.txt
grep -oE '\[@[^]]+\]' paper_prep/synthesis/NATURE_DRAFT_v1.md | tr -d '[]@' | tr ';' '\n' | sed 's/^ *//;s/ *$//' | sort -u > /tmp/draftkeys.txt
comm -23 /tmp/draftkeys.txt /tmp/bibkeys.txt
```

### Output of `comm -23` (draft keys NOT in REFERENCES_v3.bib)

```
(empty — zero leakage)
```

### Totals

| Quantity | Count |
|---|---:|
| REFERENCES_v3.bib universe | 295 unique keys |
| Draft-used keys | 101 unique keys |
| Overlap (draft ∩ bib) | 101 unique keys |
| Leaked keys (draft \\ bib) | 0 |
| Unused-in-draft bib keys | 194 |

**Citation discipline: CLEAN.** Every `[@key]` in the draft resolves to an entry in `REFERENCES_v3.bib`.

The "References" section at the bottom of the draft (lines 88–188) lists 101 keys (one per line); this matches the 101 unique keys cited inline. The References list and the inline citations are consistent.

---

## 4. Figure discipline

### Callouts found in draft

Extracted via `grep -oE '(Extended Data Fig\.|Fig\.)[^,.;()]+'` and de-duplicated:

| Figure | Panel callouts in draft | caption.md panels | Status |
|---|---|---|---|
| Fig. 1 | 1a, 1b, 1c, 1d | a, b, c, d | MATCH |
| Fig. 2 | 2a, 2b (top, bottom), 2c (left, right), 2d | a, b, c, d | MATCH |
| Fig. 3 | 3a, 3b, 3c, 3d | a, b, c, d | MATCH |
| Fig. 4 | 4a, 4b, 4c, 4d | a, b, c, d | MATCH |
| Extended Data Fig. 1 | 1a, 1b, 1c, 1d | a, b, c, d | MATCH |
| Extended Data Fig. 2 | 2a, 2b, 2c | a, b, c, d | MATCH (2d exists but uncited) |
| Extended Data Fig. 3 | 3a | a, b, c, d | MATCH (3b, 3c, 3d exist but uncited) |
| Extended Data Fig. 4 | 4c | a, b, c, d | MATCH (4a, 4b, 4d exist but uncited) |
| Extended Data Fig. 5 | 5a, 5b, 5c, 5d | a, b, c, d | MATCH |
| Extended Data Fig. 8 | 8a, 8b, 8c, 8d | a, b, c, d | MATCH |
| Extended Data nj_tree_arms | "paper_prep/figures/nj_tree_arms/" (path callout in Main L28) | (no caption.md exists) | AMBIGUOUS — no caption.md in `nj_tree_arms/` to verify panel structure; draft refers to it via path rather than a numbered Extended Data Fig. |

### ed6 / ed7 callouts

Search: `grep -E 'ed6|ed7|Extended Data Fig\. 6|Extended Data Fig\. 7' paper_prep/synthesis/NATURE_DRAFT_v1.md` → **NONE FOUND**. ✓

The draft does not call any non-existent Extended Data Fig. 6 or 7.

### Panel letter plausibility

- Fig. 1: 1a–1d match Fig1 caption ✓
- Fig. 2: 2a (allele/paralog), 2b (two-domain top/bottom), 2c (Fisher + Fst left/right), 2d (out-of-Africa) — all panels exist; "top/bottom" and "left/right" sub-references match the caption structure ✓
- Fig. 3: 3a (HG002 Pore-C contact matrix), 3b (forest plot), 3c (S_all negative control), 3d (flanking paradox + Dip-C radial inset) — all match ✓
- Fig. 4: 4a (WashU 538/494/92%), 4b (CEPH1463 11 features), 4c (RPE-1 t(X;10)), 4d (mouse zygotene rho=0.715) — all match ✓
- ed1a–d, ed2a–c, ed3a, ed4c, ed5a–d, ed8a–d all match caption.md panel letters ✓

**No AMBIGUOUS-PANEL flags raised** for numbered figures.

---

## 5. Methods completeness

Comparing draft Methods (L46–84) against the algorithmic steps in `01_pipeline.md` and `05_hic_validation.md`:

| Algorithmic step in report | Mentioned in draft Methods? | Notes |
|---|---|---|
| Contig classification (`classify_contigs.py`, p/q/pq) | NO | Draft L48 just says "extracted… from every contig classified as a p- or q-arm telomere-bearing contig". The p/q/pq three-way classification is not described. |
| Telomere trimming (telomere repeat tract removed from 500 kb window) | NO | Not in draft Methods. |
| wfmash flags `-p 95 -t 48 --quiet` | YES (L50) | ✓ |
| wfmash version `0.23.0-41-gb5f0ff1c` | YES (L50, L80) | ✓ |
| IMPG transitive closure `query -x` | YES (L52) | ✓ — but the parameter expansion / sliding window step within IMPG is glossed |
| PHR detection: `find-multichr-regions-incremental.py` sliding window 5 kb step, ≥5 alignments per chromosome, ≥2 chromosomes, ≥3 kb output | YES (L54, condensed) | ✓ — draft compresses but covers the main parameters |
| Chr18_q chimera exclusion | YES (Main L26, not Methods) | mentioned in main text only |
| PGGB `-p 95 -D /scratch` | YES (L56) | ✓ |
| `odgi similarity --all -P` | YES (L56) | ✓ |
| Leiden edge weight `w_ij = exp(-d_ij / median(d))` | YES (L60) | ✓ |
| Leiden modularity tuning details | PARTIAL (L60 mentions silhouette scan; doesn't say modularity vs silhouette balance) | the report (§01 L131) explicitly distinguishes modularity-optimized internally vs silhouette-optimized for k selection — draft conflates this |
| UPGMA `hclust(..., method = "average")` | YES (L62) | ✓ |
| NJ tree via `ape::nj()` | YES (L64) | ✓ |
| f7501 direct alignment validation (`minimap2 -x asm20`, ≥30 kb threshold) | NO | Report has a full sub-section (01_pipeline.md L157–192) — draft doesn't mention f7501 at all in Methods or Main. |
| Hi-C: HiC-Pro + Bowtie2 with MIN_MAPQ=0, RM_MULTI=0 | PARTIAL (L68 "MAPQ filters disabled") | ✓ in spirit; vendor-specific names omitted |
| Pore-C minimap2 + pairtools; CiFi minimap2 + pairtools | NO | Not named in draft Methods (draft just says "Hi-C, Pore-C and CiFi pipeline") |
| Dip-C BWA-MEM2, hickit, dip-c impute3 (4 rounds), NA12878 SNP panel | NO | Draft L72 just says "Dip-C in 16 GM12878 cells [@Tan2018] remapped to T2T-CHM13v2.0" |
| Bootstrap permutation: 10,000 permutations for B/W; W/B definition | YES (L68) | ✓ |
| Mantel: 10,000 row+column permutations | YES (L68) | ✓ |
| ARI (Adjusted Rand Index) for Hi-C-vs-sequence community comparison | NO | Not mentioned in draft Methods; report devotes a whole sub-section (§05 L107–127) |
| Per-PHR-pair vs Mantel vs sequence-pair correlation (three complementary tests) | PARTIAL | Methods L68 mentions Mantel and "W/B ratio computed by bootstrap"; sequence-pair correlation (§05 L149–178, Spearman rho 0.64–0.83 across samples) is not separately described |
| O/E inter-chromosomal normalization (E_ij = row_sum_i × col_sum_j / total_inter) | YES (L68 "observed-over-expected inter-chromosomal normalisation") | ✓ named but not formula-defined |
| Exclusion sets ×5 / resolutions ×5 | YES (L70) | ✓ |
| Single-cell 3D: hickit `--min-mapq=0`, sperm haploid mode `run_dipc_cell.sh` | NO | Not in draft |
| Pedigree `odgi untangle -e 50000 -m 1000 -j 0 -n 100 nth.best=1` | PARTIAL (L76 just says `nth-best=1`) | flag expansion absent |
| HQ filter `is_interchr=True + min_score >= 0.8 + 500 bp <= size <= 100 kb` | INCONSISTENT (L76 says 0.95/0.95) | **DIVERGES** — report uses 0.8 for HQ, not 0.95 (see numerical row 90) |
| Pattern classification script `scripts/pedigree/analyze-pedigree-recombination.py` | YES (L76) | ✓ |
| Cross-assembler intersection (hifiasm AND verkko) | YES (L78) | ✓ |

**Methods gaps (flagged):**

1. f7501 (L78442 cosmid) direct alignment validation is omitted from both Main and Methods, although it's a major confirmation of Mefford & Trask (2002) at population scale.
2. Pedigree HQ filter threshold mismatch (0.95 in draft vs 0.8 in report).
3. ARI test (independent Hi-C community detection) is not described.
4. Multi-mapping handling specifics (one randomly-chosen alignment per read) — the policy is restated, but the per-technology aligner-flag detail is lost.
5. Dip-C reconstruction pipeline (BWA-MEM2 + hickit + impute3) not named.
6. The "per-PHR-pair vs Mantel vs sequence-pair" three-test framework is not laid out; readers cannot tell which test gave which Spearman number.

---

## 6. Synthesis

### VERDICT

The draft is **broadly faithful** to the report and uses the correct framing (community structure, 3D mirroring, pedigree-resolved exchanges). Citation discipline is **clean** (0 leaked keys among 101 inline cites against a 295-key bib), figure callouts are clean (no ed6/ed7, all panel letters match captions), and the headline pedigree numbers all reproduce §14 exactly. However, the draft contains one **critical biological error** (wrong list of seven silent arms, repeated in both the Main text and the Methods), several **numerical divergences** that reflect uncritical inheritance of figure-caption numbers that no longer match the report (Lalli ρ, ed5b exclusion deltas, Dip-C flanking radial p-value), and one **structural omission** (§03 gene enrichment — DUX4L is the only specific subtelomeric gene named in the main text, with no mention of OR4F, Ambrosini blocks, hub genes, MTCO, SHOX, or RPL23AP). A reader of the Nature draft alone would have no way to discover the report's substantial gene-content story or the C4 minimal-PHR positive control finding.

**Top 5 highest-severity issues (ranked):**

1. **DIVERGES (CRITICAL): wrong list of 7 silent arms.** Draft L26 says `{chr5p, chr6q, chr7q, chr12q, chr14q, chr20p, chr20q}`. Report consistently says `{chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q}` (01_pipeline.md L124, L93, L300+ and 06_dipc_validation.md L23). Five of the seven names in the draft (chr6q, chr7q, chr12q, chr20p, chr20q) are arms with PHR signal that belong to existing communities (C11, C4, C12, C5). The same wrong list is repeated in Methods L72 as the S_all negative control composition.
2. **DIVERGES: 14 of 15 UPGMA agreement.** Draft L28 and Methods L62 say "14 of 15". Report 01_pipeline.md L153 says "12 of 15". 12/15 is the value supported by the UPGMA comparison table in the report. Either the draft is wrong or the report needs updating.
3. **OMITTED: §03 gene enrichment.** Draft does not mention olfactory receptors, Ambrosini's 11 subtelomere-specific duplicon blocks, the 6 subterminal families, hub genes (RPL23AP45 across 10 communities, SEPTIN14P22 across 9), MTCO pseudogenes in C7, SHOX in C15 PAR1, the 374 unique genes, or the biotype composition (32.1% protein-coding in C15 vs ≤9% elsewhere). DUX4L is the only specific gene name in the main text.
4. **DIVERGES: Lalli cM/Mb correlation (ρ=-0.35, n=46 → ρ≈0, n=40).** Report (07_integrated.md L103, 12_literature.md L73, L99) says ρ=−0.43 with N=39 collapsing to ρ=0.00 with N=32. The ed8c caption introduces the (46/40) panel-internal numbers, and the draft mirrors the caption instead of the report's stated values. The numbers should reconcile.
5. **DIVERGES: F_ST non-AFR range 0.02–0.04.** Report 04_heterogeneity.md L103–111 gives non-AFR pairwise Fst from −0.047 to 0.007 (paraphrased as −0.05 to 0.01). 11_summary.md L24 paraphrases as 0.00–0.02. Draft's 0.02–0.04 matches neither.

**Lower-severity but worth fixing**: HG02148 exclusion delta `0.15→0.21` (report has 0.15→0.72); NA19036 exclusion delta `0.27→0.49` (report has 0.27→0.79); 12 of 16 crossover-like in PAN028 (report table contains 13); O/E "8 of 8" tests (report O/E table has 7 rows); pedigree HQ filter 0.95 in draft vs 0.8 in report; chr13_p 88.2% cross-arm rate (not in the report's tables); piecewise model selection via AIC (report uses F-test).

### REVISION TODO (concrete edits)

- **Main L26:** Replace `(chr5p, chr6q, chr7q, chr12q, chr14q, chr20p, chr20q)` with the correct set from `01_pipeline.md L124`: `{chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q}`. (Citing 01_pipeline.md L124 and 06_dipc_validation.md L23.)
- **Methods L72:** Make the same replacement in the S_all definition. (Citing 06_dipc_validation.md L23.)
- **Main L28:** Change "UPGMA at k = 14 agrees with Leiden on 14 of 15 communities" to "agrees on 12 of 15", or reconcile with 01_pipeline.md L153 if the 14/15 figure is now correct under a revised computation. (Citing 01_pipeline.md L153.)
- **Methods L62:** Same fix as above (`agreement with Leiden 14 of 15` → `12 of 15`). (Citing 01_pipeline.md L153.)
- **Main L40 (after the RPE-1 / mouse paragraph) or new paragraph 6½:** Add a 3–5 sentence summary of §03 gene enrichment naming OR4F (the canonical Mefford prediction), hub genes (RPL23AP45 / SEPTIN14P22), MTCO-pseudogene enrichment in C7, SHOX as PAR1 protein-coding gene, and the dominant pseudogene+ncRNA composition. (Citing 03_gene_enrichment.md L8, L31, L33–39, L94–100, L104–120.)
- **Main L42 (Limitations clause iii):** Change `Spearman ρ = -0.35, n = 46` to `ρ = −0.43, n = 39` and `ρ ≈ 0 (n = 40)` to `ρ ≈ 0 (n = 32)` to align with the report. Update ed8c caption accordingly. (Citing 07_integrated.md L103 and 12_literature.md L99.)
- **Main L32 (third clause):** Change `0.02 to 0.04 within the non-AFR set` to the report's range. (Citing 04_heterogeneity.md L103–111.)
- **Main L34:** Fix the exclusion-Mantel quartet for HG02148 (0.15→0.72) and NA19036 (0.27→0.79). Update ed5b caption accordingly. (Citing 05_hic_validation.md L329–339.)
- **Main L38:** Change "Twelve of the sixteen crossover-like events are in PAN028" to "Thirteen" (or recount; the table at 14_pedigree_recombination.md L191–206 shows 13 PAN028 rows). (Citing 14_pedigree_recombination.md L191–206.)
- **Methods L76:** Change `minimum patch score 0.95 and minimum alignment score 0.95` to `min_score ≥ 0.8 with 500 bp ≤ size ≤ 100 kb` to match the report's HQ filter, OR provide explicit justification for a stricter draft-specific filter. (Citing 14_pedigree_recombination.md L20.)
- **Main L40 (or Methods L74):** Add one sentence on mouse window optimization (1Mb saturates 30/49 flanks → 2Mb 25/49 → 4Mb 19/49) to support the choice of 1Mb in the published rho=0.715 number. (Citing 08_mouse.md L124–128.)
- **Main L42 (Limitations):** Add three of the report's missing limitations: (a) somatic LCL exchange could inflate C1 (cell-line origin); (b) Dip-C cell 12 duplicate explicitly excluded; (c) hg19/T2T coordinate noise for Tan 2018 Dip-C. (Citing 10_limitations.md limitation 8, 17, 12.)
- **Main L34 (or new sentence in the integrated section):** Add the C4 (chr7_q/chr12_q) finding — "minimal-PHR positive control": only 5–25 kb shared, zero gene annotations, yet significant 3D enrichment in 4/5 diploid Hi-C samples. (Citing 11_summary.md finding 10 and 05_hic_validation.md per-community enrichment table.)
- **Methods (new clause):** Add the f7501 direct-alignment validation step (L78442 cosmid, minimap2 -x asm20, ≥30 kb threshold) given the prominence of f7501 in Mefford & Trask (2002). (Citing 01_pipeline.md L157–192.)
- **Methods L68 (Hi-C):** Add ARI sub-test (Leiden on O/E contact matrix; Adjusted Rand Index vs sequence partition; 8 samples with ARI 0.06–0.54). (Citing 05_hic_validation.md L107–127.)
- **Internal consistency note (not a draft edit, but worth flagging upstream):** Report §11 lists 12 numbered findings under the header "10 key findings" — clean up the header before the draft cites a count. The README's "27 novel contributions" vs §12 text's "24 findings" (with 27 numbered items) should also be reconciled.

### MISSING-FROM-DRAFT (ranked by importance for a Nature Article)

1. **§11 finding 6 (pseudogene-dominated subtelomeric gene content with PAR1 protein-coding exception, 32.1%).** This is one of the cleanest population-scale observations and connects directly to the telomere-position-effect literature (Baur 2001, cited in Mefford & Trask 2002). The PAR1 exception is the kind of observation a Nature reviewer would want to see.
2. **§11 finding 10 (C4 chr7_q/chr12_q minimal-PHR positive control).** A 5–25 kb tip-only shared region produces detectable 3D enrichment in 4/5 diploid samples — a strong falsifier of the "3D signal is gene-content-driven" alternative.
3. **§03 olfactory receptor and Ambrosini duplicon block story.** The 11 Ambrosini blocks → 15 Leiden community mapping is one of the most direct population-scale extensions of pre-T2T duplicon biology.
4. **§03 hub genes (RPL23AP45 across 10 communities, 21 arms; SEPTIN14P22 across 9; DDX11L16 across 9).** These are the molecular markers of the inter-chromosomal duplicon backbone.
5. **§08 mouse window-size sweep.** The fact that 19/49 mouse PHRs still saturate the window at 4 Mb (vs all human PHRs ≤500 kb) is biologically interesting and contextualizes the choice of 1 Mb for the published correlations.
6. **§09 RPE-1 self vs HPRC-community comparison.** The HPRC-vs-self ARI/Mantel comparison (Mantel 0.457 self async CiFi vs 0.548 HPRC async CiFi) shows the population partition transfers to a single individual — currently the draft only mentions t(X;10).
7. **§05 PBMC Dip-C negative result (W/B=0.983, p=0.305).** A real negative control that strengthens the GM12878 T2T finding by contrast; absent from the draft.
8. **§07 testable predictions 1 (LINC complex SUN1 W151R), 2 (haplotype-resolved 3D contacts).** Predictions 1 and 2 are testable with existing data and would help reviewers see the next steps; only prediction 3 (the Lalli analysis) made it into the draft.
9. **§02 internal (TTAGGG)n island motif composition (52.3% canonical, 47.7% degenerate variants).** Quantifies at population scale Ambrosini's variant-motif observation; draft mentions ITS islands only as boundary markers.
10. **§10 limitations 8 (somatic LCL exchange), 11 (GM12878 EBV-transformed cell line caveat), 17 (Dip-C cell 12 duplicate exclusion).** Standard hygiene items a Nature reviewer would expect.

---

*End of audit. Generated by agent-22 (consistency-audit-nature) on 2026-05-17. All status assignments backed by line citations to both files.*
