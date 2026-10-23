# Interleaved review response — every point, with status

Sources, both intact and unmodified:
- **`paper_prep/synthesis/PEER_REVIEW_v1.md`** — the actual peer review, dated **2026-05-17** (M1–M12 majors, m1–m18 minors, a ~28-row statistical table, figure adequacy, citation discipline). Verbatim gists are quoted below; the file is the authority.
- **`paper_prep/synthesis/NARRATIVE_EXTRACT.md`** — the talk Q&A (Q1, Q2, Q2b, Q3 by audience members).
- **`REVIEWS.txt`** — our *evidence log* for the pedigree objections (R1/R2/R3), created 2026-10-15. It never contained the review text; `git log -- REVIEWS.txt` shows its first commit already summarising the pedigree objections. Nothing was destroyed.

Status key: **DONE** · **DONE-EARLIER** (resolved before this round) · **PARTIAL** · **RESOLVED-BY-REMOVAL** (the claim was cut) · **NOT ADDRESSED** · **DEFERRED** (stated as future work).

Counts: majors **7 done / 1 earlier / 2 partial / 2 resolved-by-removal** · minors **2 done / 1 resolved / 4 removed / 11 not addressed-or-unverified** · statistical rows **6 addressed / 22 not addressed** · figures **2 closed** · citations **closed**.

---

## A. Major concerns (M1–M12)

| # | concern (verbatim gist) | response / status | evidence |
|---|---|---|---|
| M1 | "Reference list is unformatted … bare bibkeys" | **DONE-EARLIER.** Manuscript is a real bibtex build with three bibliographies (main/Meth/Supp); compiles to 0 undefined refs. | `submission/paper.tex`, `compile.sh` |
| M2 | "Word count is over the Nature Article cap (4158 vs ~3000)" | **PARTIAL.** Main text ≈ 3,340 words including legends; **abstract is 223 words vs the 200 cap** — still needs 23 words cut. | `submission/paper.tex` |
| M3 | "flanking-paradox claim contradicts itself … B/W strengthens at flanking; Mantel ρ weakens … paper reports only the favourable statistic" | **RESOLVED-BY-REMOVAL, and this is a gap we should close.** The flanking claim is no longer in `paper.tex` at all (`grep -c flanking` = 0), yet **our Q1 answer relies on it** (0.027 vs 0.002). The manuscript must state the flanking control *and* both metrics, or the letter cites evidence the paper doesn't contain. | action needed |
| M4 | "92% within Leiden has no null baseline" | **DONE.** 10,000-permutation null preserving per-arm marginals: published 0.872 vs null 0.677 [0.664,0.693]; corrected 0.945 vs 0.840 [0.832,0.852]; p < 1e-4; Wilson CIs. | `ANALYSIS_D_M4.json`, `scripts/pedigree/monte_carlo_null_d_m4.py`, Methods |
| M5 | "mouse ρ = 0.715 computed on 344 non-independent pairs" | **DONE.** Arm-level Mantel, 10k row-and-column permutations: ρ = 0.41/0.58/0.63/0.33 (20 kb), all p < 1e-4, null ±0.11; peak pachytene; 50 kb flat. | `scripts/hic/mantel_mouse_d_m5.py`, `mantel_mouse_d_m5.json`, Results + Methods |
| M6 | "F_ST 0.10–0.15 indistinguishable from background; compute matched genome-wide control" | **DONE.** Matched continental baselines: ratios 1.44 / 1.06 / 0.72 / 0.94, all *equivalent*; claim removed, comparison installed; block heterogeneity reported. | `scripts/popgen/matched_fst_from_blocks.py`, `matched_fst_d_m6.tsv`, Results + Methods |
| M7 | "four-link causal loop is not closed in human … abstract writes as if it is" | **DONE-EARLIER.** Claim softened in the V7 pass. Remaining "bouquet" mentions (2) are narrative, not loop claims. | `paper.tex` |
| M8 | "sample selection is opaque … which one, on what criterion" | **DONE-EARLIER.** Methods carries an "Exclusions and flank extraction" subsection. | `paper.tex` |
| M9 | "the 100% bootstrap support claim is not a real bootstrap" | **DONE, with a consequence.** Character bootstrap, B = 1000: PAR1 100%, PAR2 100%, 10p/18p 100%, acrocentric p-arms 54.4% NJ / 87.8% UPGMA; **TIGHT_q and 4q/10q are not monophyletic in the reference tree**, so they are now described as Leiden community assignments, not clades. | `char_bootstrap_d_m9.R`, `d_m9/` |
| M10 | "ED Figures 6 and 7 are missing without comment" | **DONE (and improved).** ED set now runs 1, 2, 3 with no gap (new ED3 module maps). | `paper.tex`, `ED_Fig3_module_maps.pdf` |
| M11 | "'concerted evolution' loosely defined, Dover 1982 missing" | **DONE-EARLIER.** Dover 1982 cited (Discussion); the term is used in the NAHR sense with that framing. | `paper.tex` |
| M12 | "no confidence intervals on key correlations" | **PARTIAL.** Pedigree (Wilson), F_ST (block jackknife), Mantel ρ (bootstrap), mouse arm-level Mantel (permutation) all carry intervals. **Still missing:** CIs on B/W ratios, O/E enrichment factors, and per-sample tables. | `scripts/ci/results_d_m12.json` |

---

## B. Minor concerns (m1–m18)

| # | concern | status |
|---|---|---|
| m1 | abstract 214 words vs 200 cap | **NOT ADDRESSED** — measured 223 words now |
| m2 | PAR2 length framing needs a homogenisation caveat | **DONE** — caveat present |
| m3 | cite-stack overflow (10 and 7 citations at paragraph ends) | **NOT ADDRESSED** |
| m4 | DUX4 "0–22" needs median and IQR | **NOT ADDRESSED** |
| m5 | PBMC Dip-C null (N=18) oversells a negative result | **RESOLVED-BY-REMOVAL** — PBMC/DiP-C no longer in the manuscript (0 mentions) |
| m6 | "no community-specific gene signature" floats without a takeaway | **NOT VERIFIED** |
| m7 | sandwich_same_hap and complex categories unexplained in main text | **RESOLVED-BY-REMOVAL** — pedigree pattern categories are gone with the child-side panels |
| m8 | per-arm "cross-arm sequence rate" lacks its denominator | **NOT VERIFIED** |
| m9 | Methods software list omits the IMPG version | **DONE** — IMPG v0.4.1 appears in Methods |
| m10 | Methods cites a private script path | **RESOLVED** — no `scripts/pedigree/…` path in the manuscript |
| m11 | RPE-1 t(X;10) "rediscovered" overstates novelty | **RESOLVED-BY-REMOVAL** — RPE-1 is no longer in the manuscript (0 mentions) |
| m12 | OR4F 62.1% in ED but not main text | **PRESENT** — OR4F numbers appear in text |
| m13 | "33% of mouse PHRs saturate the 1 Mb window" needs Methods + sensitivity panel | **NOT ADDRESSED** |
| m14 | "Seven limitations" but only six listed | **RESOLVED-BY-REMOVAL** — the enumeration is gone; one "limitation" mention remains, unnumbered |
| m15 | citation density (131 keys) too high | **NOT ADDRESSED** |
| m16 | "bouquet" used 5× but never defined | **NOT ADDRESSED** — 2 mentions remain, neither defined |
| m17 | 8/34/7 vs 4/28/9 architectural split needs attribution | **RESOLVED-BY-REMOVAL** — 8/34/7 no longer appears |
| m18 | Fig 1b sentence unhelpful | **NOT VERIFIED** |

---

## C. Statistical rigor table — what was addressed

**Addressed (6):**
- pedigree 92% → Monte Carlo null + Wilson CI (**DONE**)
- mouse zygotene ρ → arm-level Mantel with permutations (**DONE**)
- F_ST → matched baseline + jackknife CI (**DONE**)
- Mantel ρ = 0.66 → bootstrap/permutation interval (**DONE**)
- mouse stage series (0.574–0.715) → arm-level Mantel per stage (**DONE**)
- flanking B/W vs Mantel contradiction → claim removed (**RESOLVED-BY-REMOVAL**, to be reinstated properly per M3)

**Not addressed (22 rows), grouped:**
- **CIs on ratios**: HG002/CHM13/Pore-C B/W p-values, O/E enrichment 5.9–18.4×, per-sample sequence-pair ρ = 0.83 ("sounds cherry-picked … report all samples"), CEPH-derived comparisons.
- **Multiplicity**: per-arm Spearman ρ (39/48) without Bonferroni; piecewise-vs-linear F-test (39/41) without correction; multi-sample Mantel; C7's 156-pair subset (CI requested).
- **Nulls**: ITS-within-50 kb co-localisation (16/19) needs a permutation null over random breakpoints; C1 D4Z4 Mann-Whitney groups undefined; TAR1 PAR-vs-autosomal needs a test.
- **Independence**: Dip-C radial (2,927 vs 7,267 particles) must model intra-cell correlation; sign tests 16/16 and 20/20 should report exact p (1.5e-5, ~2e-6).
- **Denominator**: the CEPH1463 "11 features" is uninterpretable without the candidate count.
- **Reproducibility**: the 15×11 heatmap should state exact cell counts at q < 0.05 and q < 0.001; the 66.5k-pair Wilcoxon should confirm pair independence.

---

## D–F. Figures, citations, audience questions

- **Figures (§6):** ED gap **DONE**; new ED3 adds the module architecture. Remaining asks in that section were mostly about panels that no longer exist after the V7 cut.
- **Citation discipline (§7):** **DONE** — `Vollger2023` and `concerted_evolution_nahr_Vollger2023` are different papers (Science 376:eabj6965 vs Nature 617:325); Smith 1976 was already present; Cech 2004 added to the synthesis bib.
- **Q1 (MAPQ0 placement):** **DONE** — reframed around the dilution theory (uniform placement attenuates by ~1/k per end, symmetric across classes, so PHR enrichment is a *lower bound*); MAPQ-strict re-run explicitly rejected as a degenerate control (it deletes the region); placement-seed stability + orthogonal assays specified instead.
- **Q2 (RNA-soup):** **DONE** — kept as a Discussion hypothesis with direct precedent (Liang et al. 2023, complementary Alu sequences mediate enhancer–promoter selectivity) plus the sparse-participation argument.
- **Q2b (segmental-duplication generality):** **DEFERRED** — roadmap statement with the analysis required.
- **Q3 (p/q asymmetry):** **DONE and extended** — P–P 19 / Q–Q 22 / P–Q 17 against arm-type permutations (Q–Q p = 0.019, P–Q depleted p = 0.005), C6 all-Q-Q, C7 all-P-P, plus the module architecture (98.6% same telomere-relative orientation; repertoire 0.83 vs 0.57; order not conserved).

---

## G. What is genuinely still open (the honest list)

1. **M3 / flanking control**: the manuscript no longer contains the flanking comparison that our Q1 answer depends on. **Reinstate it with both metrics.** Highest priority — it's an internal inconsistency between letter and paper.
2. **Abstract length** (223 → ≤200) and **word count** (~3,340 incl. legends).
3. **CIs on the remaining ratios** (B/W, O/E, per-sample ρ) — the reviewer's blanket M12 ask is only half met.
4. **Multiplicity and nulls** for the per-arm tests, ITS co-localisation, TAR1, and the Dip-C intra-cell independence model. These are real statistical asks we have not touched.
5. **CEPH1463 denominator** for "11 features".
6. **Minors**: m1, m3, m4, m6, m8, m13, m15, m16, m18 — prose-level, mostly cheap but genuinely unaddressed.
7. **Build**: `b362563` is committed but **build-unverified** (the apt texlive install is the current blocker; `/gnu/store` is hung).
