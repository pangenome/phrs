# Point-by-point response tracker (real review)

Working document. One section per referee point, with the **verbatim ask**, our **current best
answer**, a **resolution level**, and inline `TODO:` notes describing the work needed. This is the
file to edit as work lands; the letter is assembled from it at the end.

Review: `PEER_REVIEW.md` (Nature editor + five referees). Refs #2 and #4 co-reviewed and submitted
no independent text.

## Resolution levels

| level | meaning |
|---|---|
| **L0** | unanswered |
| **L1** | answered by argument using existing evidence (no new computation) |
| **L2** | answered by re-analysing data already on disk |
| **L3** | answered by new analysis / new computation |
| **L4** | L3 plus manuscript text and/or figure changed, build verified |

Status flag: `OPEN` · `IN PROGRESS` · `DONE` · `DISPUTED` (we disagree, with reasoning).

---

# Master ranked list

Ordered by (impact on acceptance) × (feasibility with data in hand). Every item that is logically
achievable is listed.

| rank | id | ask (short) | level now | target | effort |
|---|---|---|---|---|---|
| 1 | **R1-M1** | concerted evolution is never tested | **L2 → partly done today** | L4 | 2–3 d |
| 2 | **R1-M2** | pedigree can't tell exchange from inherited haplotype | **L3 (core done)** | L4 | 1 d to finish |
| 3 | **R1-M7** | mouse needs dependency-preserving null; "zygotene bouquet" claim | **L4** | L4 | done |
| 4 | **R3-M7** | Jaccard captures similarity not order | **L4** | L4 | done (ED3) |
| 5 | **R1-M4 / R3-M3 / R3-M4** | random multi-mapping placement; flank-control logic; B/W paradox | **L1 in letter, L0 in manuscript** | L4 | 1 d |
| 6 | **R1-M5 / Q2b** | non-subtelomeric duplication control (matched identity/length/copy) | **L0** | L3 | 1–2 w |
| 7 | **R1-M8** | 500 kb right-censoring; 2–5 Mb re-extraction | **L0** | L3 | 2–3 d |
| 8 | **R1-M9** | gene enrichment rests on K=187 and one SEPTIN14 duplicon | **L0** | L3 | 3–5 d |
| 9 | **R1-M6** | mouse p-arm flanks are satellite-dominated | **L1 partial** | L3 | 3–4 d |
| 10 | **R3-M6** | 41/48 arms depend on thresholds | **L0** | L2 | 1–2 d |
| 11 | **R5-L193** | why only 2,830 pairs of 15,668 PHRs | **L0** | L1/L2 | 2 h |
| 12 | **R1-M3** | implied exchange rate vs Lemmers 2010 | **L1 partial** | L4 | 1 d |
| 13 | **R1-m11** | no population-stratified description | **L1 (F_ST done)** | L4 | 1 d |
| 14 | **R5-L86** | sensitivity of PHR count to window size | **L0** | L3 | 3–5 d |
| 15 | **R1-M10** | place the compartment in a comparative frame | **L0** | L1 | 1 d (lit) |
| 16 | **R1-m2/4/5/6/7/8/12/13, R3-m1/2, R5-L152/L229/L348** | figure/text fixes | mixed | L4 | 2–3 d total |
| 17 | **Editor compliance** | Data/Code availability, ethics, ED format, Source Data, word count, tracked changes, checklists | **L0** | L4 | 2 d |

---

# Referee #1

### R1-M1 — "Concerted evolution is never tested" · L2 (partial) · OPEN

**Verbatim ask.** "It requires that paralogous copies at different ends be more similar than their
divergence times predict, i.e. actively homogenized; high identity alone does not show this, since
recent duplication gives the same signature. The data allow a direct test: duplicon phylogenies
(D4Z4/DUX4L, TUBB8/TUBB8B, OR4F, SEPTIN14) asking whether copies cluster by end or in…"

**What we already have (stronger than we have claimed).** The allele-versus-paralog comparison in
`end-to-end-report/report/04_heterogeneity.md` §"Allele vs paralog distance" is exactly the
homogenisation signature: for each individual with two haplotypes in a multi-arm community, the
allelic distance (maternal vs paternal copy of *the same* arm) is compared with the paralog
distance (to the closest copy on a *different* arm of the same community). If paralogs are closer
than alleles, cross-arm copies are being homogenised faster than alleles diverge — which
high-identity alone cannot produce. Result: alleles closer in 8 of 9 communities (overall
p < 1e-300), and **C7 is the exception with paralogs closer than alleles in 70.5% of individuals
(p = 2e-7)** — i.e. in the rDNA-bearing acrocentric arms, homogenisation has already overtaken
allelic variation. That inversions the reviewer's premise for the strongest community.

`TODO:` promote this from report §04 into the main text as the explicit concerted-evolution test.
`TODO:` per-family version for the four families the referee names (D4Z4/DUX4L, TUBB8/B, OR4F,
SEPTIN14): cross-arm paralog identity vs within-arm allelic identity across the 232 individuals.
Data: `gene_family_copy_counting_analysis.md`, PHR gene lists, `data/chm13.phrs*.bed`.
`TODO:` duplicon phylogeny per family (copies × individuals, NJ with character bootstrap) and
report whether copies cluster **by arm/end** (homogenisation) or **by individual** (drift/descent).
`TODO:` state the dating logic explicitly (see "Perspective" below).

**Perspective on "they look like they exchanged yesterday".** That appearance *is* the prediction of
concerted evolution, so it cannot be evidence against it. What we can deliver is bounds, not dates:
(i) if the homogenised state is shared across all 232 individuals, the last homogenisation predates
the sample's common ancestor, i.e. it is not a recent event in any one meiosis; (ii) the residual
site-level differences between paralogs give a lower bound on the number of homogenisation events
since duplication (each conversion erases differences); (iii) the clock that survives is in the
flanking single-copy sequence and at breakpoints — which is exactly why the pedigree work (R1-M2)
asks about junctions rather than identity. Our module-architecture result (repertoire shared,
orientation conserved at 98.6%, order scrambled) is the expected signature of repeated conversion
acting on a shared module set rather than of one recent duplication.

### R1-M2 — pedigree: exchange vs inherited haplotype · L3 · DONE (finish to L4)

**Verbatim ask.** "The pedigree evidence does not distinguish a new exchange from an inherited
haplotype… An exchange claim needs the tract shown absent from both paternal copies of the recipient
[and: (1) chr5q/chr1p called both father-to-daughter and daughter-to-granddaughter; (2) all three
Fig 5C calls on the same paternal haplotype; (3) identity margin ~1% with comparable dips outside
the called tract]. We ask for the comparison against both parental haplotypes, a genome-wide
distribution of identity margins locating the candidates, a false-positive estimate from permuted
parent-child pairs, and confirmation that each tract lies in an unbroken contig."

**What we did.** Both-parental-haplotype comparison (father's chr3q h2 at 99.79–99.84% with 100%
coverage; his chr9q h1 excluded at 7.3–7.5% coverage); per-chromosome maternal-copy resolution
(the earlier child-side calls used the paternal copy — 73% of them); threshold sweep f10–f8192;
HiFi flagger cleanliness (1/538). Their three sub-points are all confirmed rather than refuted:
chr5q/chr1p is ancestral homology, the "two events" is transmission, and margins are ≤1%.

`TODO:` genome-wide distribution of identity margins with the candidates located on it.
`TODO:` permuted parent–child false-positive rate (shuffle which child haplotype is compared).
`TODO:` contig-breakpoint check for the retained tract (unbroken contig across the junction).

### R1-M7 — mouse null and the "zygotene bouquet" claim · L4 · DONE

Arm-level Mantel, 10,000 row-and-column permutations: ρ = 0.41/0.58/0.63/0.33 (20 kb),
all p < 1e-4, null 0 ± 0.11. Peak moved to **pachytene**; the bouquet claim was removed from the
abstract/introduction and replaced by the synapsis framing. `TODO:` none, except R1-M6 below.

### R1-M4 / R3-M3 / R3-M4 — random placement, flank-control logic, B/W paradox · L1 · OPEN

`TODO:` reinstate the flanking unique-sequence control in the manuscript **with both metrics**
(B/W 0.027 vs 0.002 at 50 kb; Mantel ρ 0.657→0.520) and state the dilution logic explicitly:
uniform placement attenuates a specific contact by ~1/k per end (~1/k² for PHR–PHR), symmetrically
across community classes, so PHR-window enrichment is a **lower bound** and the artefact hypothesis
predicts the reverse of what we observe.

### R1-M5 / R3-M3 — matched non-subtelomeric duplication control · L0 · OPEN

`TODO:` pericentromeric + interstitial duplications matched on identity, length and copy number,
through the identical pipeline; report whether subtelomeric ρ = 0.66 exceeds the all-SD value.

### R1-M8 — 500 kb right-censoring · L0 · OPEN

`TODO:` re-extract saturated ends with 2–5 Mb flanks; report the censored fraction; restate 3.51 Mb
as a lower bound.

### R1-M9 — gene enrichment on K=187 · L0 · OPEN

`TODO:` rework around the 215 unannotated copies; quantify the SEPTIN14 dependence of the 25-fold
spermatid term; report enrichment conditional on the largest family.

### R1-M6 — mouse p-arm satellite confound · L1 partial · OPEN

`TODO:` satellite-masked (minor/major) mouse flank re-analysis; sensitivity of ρ to masking.

### R1-M10 — comparative frame · L0 · OPEN

`TODO:` short comparative paragraph (apes, yeast, Phaseolus, Plasmodium VAR) — literature only.

### R1-M3 — rate vs Lemmers 2010 · L1 partial · OPEN

`TODO:` explicit calculation showing our conclusion (one ancestral arrangement, transmitted) is
consistent with four 4q/10q transfers over human evolution.

### R1 minors

- **m2** scope statement: which analyses are pangenome-wide vs CHM13-anchored. `TODO:`
- **m4** Fig 1A vs 1B disagreement at chr3p/chr8q. `TODO:` reconcile or add a PHR track.
- **m5** distinguish rDNA-driven acrocentric PHRs; restore chrY/14/21 centromere tracks. `TODO:`
- **m6** quantify "concentrated at ends" against a segdup-density null. `TODO:`
- **m7** Fig 4A vertical-axis points; per-assay B/W labelling; pair-inclusion rule. `TODO:`
- **m8** Fig 5 number inconsistencies (PAR1 138 vs 92 kb; 28 kb vs other length). `TODO:` — panels
  were rebuilt; reconcile numbers rather than delete.
- **m11** population-stratified descriptives (length, community membership by superpopulation).
  `TODO:` cheap; our F_ST work already provides the analysis for the differentiation part.
- **m12** engage Xiao 2026 (D4Z4) and Karimian 2024 (telomere length). `TODO:`
- **m13** title: "unorthodox" → ectopic/non-allelic. `TODO:`

---

# Referee #3

- **M1** novelty framing vs Linardopoulou 2005 / Guarracino 2023. `TODO:` framing paragraph.
- **M2** Fig 1A positional-profile test. `TODO:` quantify profile similarity across ends.
- **M3** random assignment could create contacts — see R1-M4. `TODO:` same fix.
- **M4** B/W direction paradox — see R1-M4. `TODO:` same fix.
- **M5** somatic vs meiotic; human meiosis gap. `TODO:` consistent limitation statement.
- **M6** 41/48 arms depend on thresholds. `TODO:` threshold sweep for the PHR call (identity,
  window, alignment counts) — data exists.
- **M7** similarity vs order/colinearity. **DONE** — order concordance, orientation, ED3.
- **M8** call pedigree findings candidate exchanges consistently. **DONE** — one structurally
  supported event; recurrence withdrawn.
- **minors** Fig 3 crop at high resolution; define NOR. `TODO:`

---

# Referee #5

- **L86** window-size sensitivity of the 15,668 count. `TODO:` re-run calling at 50 kb.
- **L152** stability of copy-number-aware GO across assemblies. `TODO:`
- **L189** which Hi-C/Pore-C datasets and cell types. `TODO:` explicit Methods table.
- **L193** why 2,830 pairs from 15,668 PHRs. `TODO:` state the filter; 2 h.
- **L229** Fig 5 legend dark-grey exchanges undefined. `TODO:` moot after the panel rebuild —
  re-check the current legend.
- **L348** why 11.6% sampling rather than 3.5%. `TODO:` one paragraph.

---

# Editor's compliance list (not scientific, all mandatory)

`TODO:` Data availability statement; Code availability statement (DOI-minting repository, e.g.
Zenodo); ethics statement; Extended Data as `.jpg/.tif/.eps` (currently PDF) and ≤10 items;
Source Data per figure; final word count in the cover letter; tracked-changes version; updated
Reporting Summary and Code/software checklist.

---

# Flow — the order this actually happens in

**Step 0 (now).** Restore the as-submitted manuscript as the base (`cc7921d`), then re-apply our
seven change-groups as isolated hunks. Reason: referees quote text we deleted (flanking control,
Fig 5 numbers, pattern categories), so a splice must come before any further text editing.

**Step 1 (cheap, high impact — ~3 days total).** R5-L193 filter explanation; R1-m11 superpopulation
descriptives; R1-M2's three completions; R3-M6 threshold sweep; editor compliance statements;
reinstatement of the flank rationale (R1-M4).

**Step 2 (the two heavyweight intellectual asks — ~1 week).** R1-M1 concerted-evolution test
(allele-vs-paralog promoted to main text, per-family trees, dating bounds) and R1-M8 re-extraction
at 2–5 Mb.

**Step 3 (controls that referees treat as decisive — 2–3 weeks).** R1-M5/R3-M3 matched SD control;
R1-M6 satellite-masked mouse; R1-M9 gene-enrichment rework.

**Step 4 (text/figures and the letter).** All figure and minor fixes, then assemble the point-by-point
letter from this file, with tracked changes.

**Step 5 (verification).** `bash compile.sh` → exit 0, `grep -c undefined paper.log` = 0, and
`python3 scripts/ci/audit_response_numbers.py` for every quoted number.

---

# Lane results (2026-10-23) — five parallel lanes, all completed

Notes: `response/notes/{M1_concerted_evolution,M2_pedigree_completions,R5_pairs_filter_and_R3_thresholds,D1_superpopulation_descriptives,D2_editor_compliance,M6_mouse_satellite}.md`

## R1-M1 — concerted evolution · now **L3**, and the answer is negative for the named families

Substrate definitions pinned before any number: **P** = the telomere-anchored PHR window itself (carries the duplicons); **F** = the adjacent 100 kb of *unique* centromere-proximal sequence. Because F's "paralog" pair is two non-paralogous unique loci, F cannot test homogenisation at all — its 0.18 distance is a floor. **Only P is informative.**

Substrate P, alleles vs paralogs: D4Z4/DUX4L paralog-closer 36/85 (0.42, p=0.19); TUBB8/B **0/197** (p=9.96e-60); OR4F 272/1144 (0.24, p=1.08e-73); SEPTIN14 843/3127 (0.27, p=8.25e-152). On the referee's own formulation, within-individual cross-arm copies are **2.22× / 11.45× / 2.95× / 4.57×** farther apart than same-arm copies from *other* individuals. Sketch distances validated against FastGA alignments (Spearman +0.77, p=8e-25, n=120; ordering preserved, sign p=0.0026).

**Conclusion to write up:** no homogenisation signature for these four families at whole-window resolution. This does **not** contradict the paper: the aggregate's single homogenisation signal is community **C7** (paralog closer in 70.5%, p=2e-7), and **none of the four families has a copy on a C7 arm**. The M1 answer is therefore "homogenisation is demonstrated in the rDNA-bearing acrocentrics, and is not detectable in these gene families" — stated with the dating caveat (no exchange is dated; shared state bounds it as pre-dating the individuals' common ancestor).

## R1-M6 — mouse satellite confound · L3, but the stratified result is RETRACTED (void axis)

> The satellite-masked axis is broken (QLEN/QSTART swap at `mouse_satellite_sensitivity.py:706`; 4,228/4,228 records reversed, all derived Jaccard = 0). The "pachytene peak does not survive in satellite-poor windows" claim is withdrawn. What stands: the per-arm satellite FRACTIONS (p-end 0.26-0.99, mean 0.885; q-end 0.004) from annotation, and Corr(satellite, published Jaccard) = +0.790 on the real axis. Original (superseded) text follows.

- p-end flanks are **0.26–0.99 satellite** (38/39 ≥ 0.5, mean 0.885); q-end flanks **0.002–0.008**. The referee is right.
- Major satellite `GSAT_MM` is a real Dfam annotation; **minor satellite is a proxy** (no Mus family in the Dfam library) and the proxy recovers ~0, so "minor satellite does not contribute" is **unmeasured**.
- Stratified: the association **survives in satellite-poor windows** (rho +0.49 to +0.64 at 20 kb, p 0.005–0.025) but the **pachytene peak does not** — the low-satellite subset is flat across stages (0.536 / 0.539 / 0.486 / 0.486).

**MANUSCRIPT CONSEQUENCE (action):** the pachytene-peak framing we added to the Introduction, Results, Fig4c caption and Discussion is a satellite-rich-p-arm phenomenon. It must be qualified: report the whole-flank profile, state that it is dominated by satellite-rich p-arm windows, and report that the satellite-poor subset is flat across stages. This supersedes the "sequence similarity translates into ectopic proximity precisely during the window in which crossovers form" sentence.

## R1-M2 — pedigree completions · **L4-eligible**

Margin distribution: 177,998 windows; only 1.02% have margin > 0; the retained tract sits at the 99.0–99.7th percentile but **558 windows exceed its max margin**, and the withdrawn chr5q/chr1p call sits at 99.1–99.6th — margin does not discriminate (confirms the referee, supports our structure-based argument). Permutation FP: null A p<1e-4; null B genome-wide p<1e-4 but candidate-level p=0.09 (chr9) / 0.32 (chr5). Contig: 23/23 lengths identical, zero N-runs, tract 21,149 bp inside `PAN027.chr9.paternal`.

## R5-L193 + R3-M6 — **answered**

2,830 = C(77,2) − 96 same-chromosome pairs (object identity, not similarity filtering), verified two ways and matching `paper.tex`'s n and ρ. Threshold gate reproduced the published call exactly (15,668 PHRs, 41 arms, all 18,826 region sets identical): consecutive-fail 1–4 → 41 arms; flank 25–400 kb → 41 arms; chromosomes-per-window 2/3–5/6–12/15 → 41/40/38/37; no inter-chr rule → 48 arms (the seven the paper names).

## R1-m11 + editor compliance — **D1 done, D2 drafted**

D1: per-superpopulation PHR count, bp, median length + IQR, community shares, and the 500 kb right-censoring share; superpopulation mapping is repo-local (`slides/.../nearest_same_superpop_mds_distances.tsv`), reproduction gate md5-identical. D2: Data/Code availability, ethics split into attestable vs **AUTHOR TO SUPPLY**, format checklist (ED are 3 items but all `.pdf` → need `.tif/.eps`), and it flags a repo URL inconsistency (`ekg/phrs` vs `pangenome/phrs`).
