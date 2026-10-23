# L5 — Data consistency reconciliation + LaTeX patch snippets

> Provenance for this note: see `response/PROVENANCE.md` (tool versions: `response/notes/tool_versions_and_pins.md`).

STATUS: COMPLETE. Four deliverables written; appending log kept in step order.

Scope: I own the data-consistency reconciliation and three LaTeX patch files.
I do NOT edit submission/paper.tex. No git operations. No invented numbers.

Deliverables:
1. response/notes/L5_data_reconciliation.md  (this file)
2. response/patches/flanking_control.tex
3. response/patches/pedigree_pattern_categories.tex
4. response/patches/mouse_satellite_caveat.tex

---

## Step 0 — log of actions (append-only)

1. Skeleton created; explored the repo layout (`ls`, `PEER_REVIEW.md`, `response/POINT_BY_POINT.md`).
2. Re-derived the BED-vs-published-call discrepancy from the tables themselves (step 1).
3. Located every consumer of `data/chm13.phrs.bed` and measured the downstream consequence (steps 1.4, 1.5, 4).
4. Re-verified every flanking-control number against the raw `global_test.tsv` files (step 2).
5. Wrote the three patch files and mechanically balance-checked them (step 3).
6. Wrote the figure-by-figure consequence table, the recommendation and the residual risks (steps 4-6).

## Step 1 — re-derived the vendored-vs-published discrepancy from scratch

All numbers below are from commands I ran in `/home/erikg/phrs` on 2026-10-23; no number is copied
from another lane's note unless I say so.

### 1.1 The two call sets

| artefact | rows | source |
|---|---|---|
| `data/chm13.phrs.bed` | 37 intervals | vendored; matches the **earlier** `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv` call. I recounted that file myself: `rows=18225 signal=15088 cols=6` (6 columns, no `arms_involved`) |
| published call | 38 intervals | reconstructed from `data/all-vs-all.1Mb.p95.id95.len.tsv`: `NR>1 && $1 ~ /^CHM13/ && $4!="."` → 38 rows |

`awk -F'\t' 'NR>1 && $1 ~ /^CHM13/ {n++} END{print n}' data/all-vs-all.1Mb.p95.id95.len.tsv`
→ **46** CHM13 flank rows; of those **8** have `region_start/end == "."` →
**38** CHM13 PHR intervals in the published call.

### 1.2 Exact per-arm comparison (published interval = `self_chr_start + region_start` .. `+region_end`)

Two systematic classes, and they must not be conflated:

* **Class A — coordinate convention only.** 22 of the 37 shared intervals differ by exactly 1 bp at
  the end (`e1`). The vendored BED mixes conventions (`start = st0+rs`, `end = st0+re-1` for some
  rows and `st0+re` for others). This is *not* a different call and has no biological content.
* **Class B — genuinely different interval.** 8 published intervals are not faithfully represented,
  **7** of them among the 37 shared arms plus **1** published arm with no vendored row at all:

| arm | vendored BED | published call | Δ length |
|---|---|---|---|
| chr13_p | **absent** | 2544–502544 (500,000 bp) | **+500,000** |
| chr14_p | 2075–372075 (370,000) | 2075–502075 (500,000) | +130,000 |
| chr15_p | 3258–408258 (405,000) | 3258–503258 (500,000) | +95,000 |
| chr7_q | 160560602–160565601 (5,000) | 160520602–160565602 (45,000) | +40,000 (start moves 40 kb proximal) |
| chr17_p | 1679–16679 (15,000) | 1679–46679 (45,000) | +30,000 |
| chr8_p | 1881–11881 (10,000) | 1881–31881 (30,000) | +20,000 |
| chr21_p | 2505–487505 (485,000) | 2505–502505 (500,000) | +15,000 |
| chr6_q | 171978883–172123882 (145,000) | 171983883–172123883 (140,000) | −5,000 (start moves 5 kb distal) |

**Correction to the inherited count.** The prompt / `response/notes/R5_pairs_filter_and_R3_thresholds.md`
say "6 of 37 CHM13 intervals differ". My re-derivation gives **8 of 38** (7 shared arms + chr13_p
which is missing entirely). The two named examples in that note (chr14 end 372,075 vs 502,075;
chr8 11,881 vs 31,881) reproduce exactly, so the comparison is the same one — the count of 6 is an
undercount that misses chr13_p and the chr6_q / chr7_q start shifts.

### 1.3 Aggregate size

`awk` sums of interval lengths in the *same* convention:

* vendored: n=37, **6,014,981 bp** (6.015 Mb) — reproduces the hard-coded
  `EXPECTED_TARGET_BP = 6_014_981` in `paper_prep/_brainstorming/chm13_copy_enrichment/build_inputs.py:29`
* published: n=38, **6,839,962 bp** (6.840 Mb in true 1-based inclusive length)
* **Δ = +824,981 bp, +13.72 %**

Length distribution (sorted, index-based quantiles):

| | n | min | q1 | median | q3 | max | rows at the 500 kb flank ceiling |
|---|---|---|---|---|---|---|---|
| vendored | 37 | 4,999 | 25,000 | 124,999 | 235,000 | 500,000 | **2** (chr22p, chrXp) |
| published | 38 | 9,999 | 44,999 | 124,999 | 254,999 | 499,999 | **6** (chr13p, chr14p, chr15p, chr21p, chr22p, chrXp) |

The saturated fraction goes **2/37 = 5.4 % → 6/38 = 15.8 %**: the vendored BED understates
right-censoring at the 500 kb flank by about 3×. That is directly Referee #1 major point 8.

### 1.4 Which consumers actually read `data/chm13.phrs.bed`

I grepped every non-worktree `.py/.R/.sh/.md` for `chm13.phrs`, then read each hit.

| consumer | uses the BED for | does the plotted/headline number change? |
|---|---|---|
| `submission/scripts/figures/make_fig2bc_jaccard_heatmaps.R:28` | only the side-table `arm_inclusion_audit.tsv` (`chm13_phr_bed_interval_count`, status strings) | **No.** The Fig 2b/2c heatmap values come from `data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv` (md5 `3db1f608…`, byte-identical to the moosefs `similarity/` product, 41 arms) + the Leiden table. Neither is built from the BED. Only the audit annotation changes; the arm assignment of every interval is unchanged. |
| `submission/scripts/figures/make_fig5_pedigree_supporting_material.py:34` | the `chm13_query_phr_bed_interval` / `chm13_donor_phr_bed_interval` text columns of `submission/supporting_material/Fig5_pedigree_supporting_material.tsv` | **Yes, 3 of 5 distinct arms.** Current strings: chr15_p `chr15:3258-408258`, chr21_p `chr21:2505-487505`, and **chr13_p `no called interval in data/chm13.phrs.bed`** → published gives `chr13:2544-502544`. The chr13_p row is the acrocentric rDNA exchange — the one Referee #5 (L229) asked to have explained — and the stale BED currently makes it read as an arm with no PHR. chr5_q, chr9_q, chrX_p are unaffected. |
| `paper_prep/_brainstorming/chm13_copy_enrichment/build_inputs.py:632` | the entire copy-number-aware gene-enrichment input universe | **Yes.** See 1.5. |
| `map_all_protein_genes.py:177` | `chm13.phrs.no_acro.bed` (29 rows, drops whole acrocentric chromosomes) | Yes, via the same 5 affected chromosomes/arms; the acro-excluded set also loses the newly added chr13_p entirely. |
| `scripts/pedigree/crossover_rate_f34.py:23` | "CHM13 PHR length (chm13.phrs.bed)" = the per-meiosis rate denominator | Yes in principle: 6.015 → 6.840 Mb, i.e. **the rate falls by 12.1 %**. **But the script as committed cannot run**: `PHR_BED = REPO_ROOT / "chm13.phrs.bed"` and `ls chm13.phrs.bed` → *No such file or directory*. The vendored file is at `data/chm13.phrs.bed`. So this provenance defect is currently masked by a path bug. |
| `scripts/popgen/matched_fst_d_m6.py:229` | a comment inside an unexecuted pseudocode stub | No. |
| `permutation_based_enrichment_investigation.md` (repo root) | a *proposed* `bedtools shuffle` protocol; never executed | No, but see 1.5: the real permutation nulls are the burden nulls inside `COPY_engine.py`, which consume the tables built from the BED. |

`make_fig1a_genomewide.R` reads `hprc25272-wf.CHM13.100kb-xm5-id098-l50000.tsv.gz`;
`make_fig1b_lengths.R` reads `data/all-vs-all.1Mb.p95.id95.len.tsv` (default `data_path`).
**Fig 1a and Fig 1b are therefore already on the published call** — they are not affected.

### 1.5 Measured consequence for the gene-enrichment universe (Fig 3 / ED2 / paper.tex:696)

I ran the real builder against both BEDs (`build_inputs.py` imported as a module; `load_phrs` and
`load_gene_loci` with the hard-coded expectations relaxed; GFF3 `data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz`,
cytobands `data/chm13v2.0_cytobands_allchrs.bed`):

```
VENDORED  n_phrs=37 total_bp=6014981 n_arms=37
          gene_loci_total=61312 midpoint_assigned=402 overlap_assigned=412
PUBLISHED n_phrs=38 total_bp=6840000 n_arms=38
          gene_loci_total=61312 midpoint_assigned=404 overlap_assigned=414
          shared midpoint loci 402, only-vendored 0, only-published 2
          shared overlap  loci 412, only-vendored 0, only-published 2
```

The vendored run reproduces the manuscript's numbers **exactly** (402 / 412 / 61,312), which pins the
provenance: `paper.tex:696`'s "402 gene copies with midpoints in CHM13 PHRs" *is* the vendored-BED
count. Published-call consequences:

* 402 → **404** midpoint-assigned (+0.50 %), 412 → **414** any-overlap (+0.49 %)
* the two added loci: `SCGB1C2` (chr17:43009-44511, chr17_p, **protein_coding**, from chr17_p
  15 kb → 45 kb) and `RPL23AP53` (chr8:5120-29123, chr8_p, **transcribed_pseudogene**, from chr8_p
  10 kb → 30 kb). No locus is lost: the only shrink (chr6_q, −5 kb) removes nothing.
* Result held under both BED end-conventions (`end = st0+re` and `end = st0+re-1`): 404 / 414 both
  ways, so this is not a ±1 artefact.
* `K = 187` (the annotated subset of the 402) therefore changes by 0, 1 or 2 depending on whether the
  pipeline's pseudogene filter admits a *transcribed* pseudogene. **I did not re-derive K**: the
  `analysis_ready/` and `results/v1/` products are not vendored in this worktree
  (`ls paper_prep/_brainstorming/chm13_copy_enrichment/analysis_ready` → does not exist), so I record
  the bound and do not invent the split.
* The 99,999-permutation burden nulls were **not** re-run. Reason, quoted from
  `paper_prep/_brainstorming/chm13_copy_enrichment/REPRODUCIBILITY_README.md`: "about 18.08 GiB for
  the transient batch products", "moosefs worker environment with `guix` … and Slurm `sbatch`", and
  "Do not run the full permutation workload on a login node". Those nulls are spatial samplers keyed
  on the PHR intervals, so they are *sensitive* to the +13.7 % target-bp change and to the
  2→6 saturated arms, not only to the +2 loci. That is a real open item, not a resolved one.

### 1.6 Additive check

The 8 Class-B deltas sum to +500,000 +130,000 +95,000 +40,000 +30,000 +20,000 +15,000 −5,000 =
**+825,000**, against the aggregate **+824,981** measured in 1.3 (the −19 residual is 19 of the 37
Class-A ±1 rows). The two independent routes agree, so the per-arm table in 1.2 is complete.

## Step 2 — flanking control numbers re-verified (before writing the patch)

Sources: `end-to-end-report/report/05_hic_validation.md` (tables at :29–41, :83–92, :95–105,
:257–279) plus the raw `*_global_test.tsv` files on moosefs and the vendored
`data/hg002_hic_{,flanking_}global_test_20000bp.tsv`.

Two facts I had to establish before I could write an honest sentence, because the brief's
phrasing and the report's tables disagree in two places:

1. **The manuscript's current pair (0.0013 vs 0.028, "21-fold") is the 20 kb row, not 50 kb.**
   `data/hg002_hic_global_test_20000bp.tsv` gives within_mean 0.10514836, between_mean
   0.00294190 → B/W = 0.027979; `data/hg002_hic_flanking_global_test_20000bp.tsv` gives
   0.11078114 / 0.00014418 → 0.0013014. Those are exactly `paper.tex:294`'s 0.028 and 0.0013.
2. **The 50 kb fold is 12.4, not 13.5.** 0.026672/0.0021438 = 12.44. The 13.5 in the letter draft
   (`response/02_letter_draft_v1.md:40`) and in the brief is 0.027/0.002 from *rounded* inputs.
   Likewise "10 kb: 0.029 vs 0.001, 29-fold": 0.029 is the **5 kb** PHR value
   (`05_hic_validation.md:269–275`), the 10 kb PHR value is 0.028, and the exact 10 kb fold is
   26.4. I quote exact numerators/denominators at a stated resolution instead of rounded folds.

Full sweep, B/W = between_mean/within_mean, all from the 12 `*_global_test.tsv` files
(n_within/n_between identical between PHR and flank windows: 256/2,670 HG002, 53/650 CHM13):

| res | sample | B/W PHR | B/W flank | fold | Mantel PHR | Mantel flank | Δρ |
|---|---|---|---|---|---|---|---|
| 10 kb | HG002 | 0.028429 | 0.0010758 | 26.4 | 0.7052 | 0.3789 | 0.326 |
| 20 kb | HG002 | 0.027979 | 0.0013014 | 21.5 | 0.6484 | 0.4484 | 0.200 |
| 50 kb | HG002 | 0.026672 | 0.0021438 | 12.4 | 0.6573 | 0.5203 | 0.137 |
| 10 kb | CHM13 | 0.072496 | 0.071384 | 1.02 | 0.7166 | 0.3420 | 0.375 |
| 20 kb | CHM13 | 0.074439 | 0.057047 | 1.30 | 0.7209 | 0.4748 | 0.246 |
| 50 kb | CHM13 | 0.070569 | 0.056652 | 1.25 | 0.6563 | 0.5217 | 0.135 |

Mantel p = 0.0 (i.e. below the 1e-4 resolution of 10,000 permutations) in all twelve.

**The instability I must report as the finding.** The B/W flank control *does not reproduce
across haplotypes*: 12.4× in HG002, 1.25× in CHM13 at 50 kb. The rank metric does reproduce:
flank/PHR ρ = 0.792 (HG002) and 0.795 (CHM13) at 50 kb. So the claim "the flanks are the more
within-community-enriched of the two" is HG002-specific, and a conclusion that holds under one
metric and not the other is exactly the kind of thing the brief told me to surface rather than
smooth over.

**Independent corroboration of the report's absolute levels.** I also read the absolute means,
which the report does not print: at 50 kb HG002 the within-class mean contact is *higher* in the
flanks (0.1076) than at PHR coordinates (0.0583), while the between-class mean is 6.7× *lower*
(2.3e-4 vs 1.6e-3). So the flank B/W of 0.002 is carried by the between-class term, not by the
within-community contact being larger in unique sequence. That is the mechanism the patch states,
and it is why B/W cannot be read as "PHR windows are more enriched".

## Step 3 — patch files written

Three files, all new, nothing outside `response/` touched, no git commands run.

| file | replaces | status |
|---|---|---|
| `response/patches/flanking_control.tex` | `submission/paper.tex:291–296` (verbatim excerpt in the comment block) + optional Methods replacement for `:764–767` | written; brace/dollar balance checked |
| `response/patches/pedigree_pattern_categories.tex` | nothing verbatim — the terms `sandwich_same_hap` / `complex` do **not** occur anywhere in `paper.tex` (grep), so this is a re-insertion at a named anchor (`:346–349`) | written |
| `response/patches/mouse_satellite_caveat.tex` | `paper.tex:306–312` (Results), `:81–83` (summary paragraph), `:545–547` (Fig 4c caption), each quoted verbatim | written |

Verification done on all three: `$` parity even, brace depth 0, and the only macros used are
`\emph`, `\rho`, `\times` (i.e. nothing that needs a package the manuscript does not already load).
No LaTeX binary exists on this box (`ls /usr/bin | grep latex` → nothing; `submission/compile.sh`
expects `pdflatex`), so I could not compile-test the snippets; I checked balance mechanically instead.
The `mouse_satellite_caveat.tex` R3 snippet carries one intentionally unbalanced `}` because its
excerpt ends by closing the `\caption{...}` group — flagged in the file.

### 3.1 Definitions authority for the pattern categories

I did not paraphrase from the report's short glosses; I read the classifier.
`/moosefs/guarracino/HPRCv2/scripts/pedigree/analyze-pedigree-recombination.py:236–249`:

```
if not left_is_query or not right_is_query:  -> 'complex'
elif left_hap == right_hap:
    if ref_hap == left_hap                   -> 'sandwich_same_hap'
    else                                     -> 'gene_conversion_like'
else                                         -> 'crossover_like'
```

plus the second pass at `:265–269` that relabels `complex`→`acros_like` inside acros-like flanks, so
the surviving `complex` count is post-relabelling.

Counts: `sandwich_same_hap` 115, `complex` 1, both of **494** within-community WashU HQ inter-chr
patches (`end-to-end-report/report/14_pedigree_recombination.md:35,43–45`; the five categories sum
to 494, checked; and `scripts/pedigree/monte_carlo_null_d_m4.py:503–505` reaches the same census
independently: 30 + 199 = 229 acros_like).

### 3.2 Mouse caveat — what I refused to write

I recomputed the satellite fractions from `scripts/hic/mouse_satellite_sensitivity.json`
(`satellite_annotation.per_flank`, 78 flanks) rather than trusting the note, and three obvious
sentences turned out to be false, so the patch says something weaker and true instead:

1. **"the p-arm satellite is major satellite"** — false alone. GSAT_MM covers a mean 0.2886 of the
   p-end flank; the union of eight Dfam families covers 0.8848. By the union bound the other seven
   families cover **at least** 0.596, which is what the text says (it does not claim exact
   additivity, because the union is not a sum).
2. **"masking removes / preserves the association"** — no masking was done. `response/notes/M6_mouse_satellite.md` §5a and the JSON `not_established` both record that the contact axis
   cannot be masked and that the re-derived sequence axis is void (a `qs > qe` bug made every
   intersection empty). What exists is a *stratification* by flank satellite fraction. The patch
   says "setting aside the satellite-rich windows" and never "masked".
3. **swap pachytene for zygotene** — refused. The satellite-poor 20 kb profile is
   leptotene 0.5364 / zygotene 0.5390 / pachytene 0.4857 / diplotene 0.4857, so its nominal top
   stage is zygotene by **0.0026**, smaller than the spread introduced by a 200-permutation null
   (§5c), and the *same* flip happens in the satellite-rich subset (zygotene 0.6515 vs pachytene
   0.6394). Swapping one unstable rank for another would repeat the original error in the opposite
   direction, so the patch states a flat profile and no maximum.

Also recorded in the file: the "satellite-poor" subset **is** the q–q subset, 21 pairs among 7 arms,
because arm satellite fraction is bimodal with nothing between 0.10 and 0.25 — so it is not an
independent sample, and the text says so.

## Step 4 — downstream consequence, measured, consumer by consumer

### 4.1 `submission/supporting_material/Fig5_pedigree_supporting_material.tsv` (Fig 5 supporting table)

Six data rows, five distinct arms. Current annotation columns, and what the published call gives:

| query arm | current `chm13_query_phr_bed_interval` | published call | change |
|---|---|---|---|
| chr13_q | `chr13:113553247-113563246` | same (±1 bp) | none |
| chr13_p | `no called interval in data/chm13.phrs.bed` | `chr13:2544-502544` | **absent → 500 kb interval** |
| chr15_p | `chr15:3258-408258` | `chr15:3258-503258` | +95 kb |
| chr21_p | `chr21:2505-487505` | `chr21:2505-502505` | +15 kb |
| chr5_q | `chr5:181848925-182043924` | same (±1 bp) | none |
| chr9_q | `chr9:150390025-150615024` | same (±1 bp) | none |
| chrX_p | `chrX:1579-501579` | same (±1 bp) | none |

Three of the five distinct arms change, and the chr13_p row is the acrocentric rDNA exchange — the
dark-grey exchange Referee #5 (L229) asks us to name. The stale BED currently makes that row read as
an arm with no PHR at all.

### 4.2 `submission/fig/MainFigures/arm_inclusion_audit.tsv` (committed artefact)

The mismatch is **already printed in a committed file**, so no computation of mine is required to
see it: summing the two columns gives `chm13_phr_bed_interval_count` = **37** against
`projected_chm13_nonempty_rows` = **38**, and exactly one arm disagrees —
`chr13_p bed_count=0 projected_nonempty=1 "non-empty CHM13#0 projected PHR interval present"`.
Only the per-arm counts and status strings change; the Fig 2b/2c cell values are the 41×41 matrix
and are untouched (see 1.4).

### 4.3 ED2 term table — exact sensitivity, not an estimate

`submission/supporting_material/ED_Fig2_supported_ontology_terms.tsv` has **209** data rows and every
one carries `K_phr_total = 187`, `nonphr_total = 31,779`, `N_eligible_total = 31,966` — i.e. the
whole table is keyed on the vendored call. I reproduced all 209 `p_exact_upper` values from
`a_T, M_T_genome_burden, K_phr_total, N_eligible_total` with the one-sided hypergeometric upper
tail to a maximum absolute error of **1.8e-17**, which validates the recomputation; then re-evaluated
each term under the published call (N' = N+e, K' = K+e, e ∈ {1,2} = the two new eligible copies):

| scenario | max abs Δp over 209 terms |
|---|---|
| e = 1, new copy carries no term | 1.74e-07 |
| e = 2, new copy carries no term | 3.69e-07 |
| e = 1, new copy carries every term (upper bound) | 1.44e-06 |
| e = 2, new copy carries every term (upper bound) | 1.47e-06 |

The largest p in the table is **1.47e-06** and the smallest is 2.87e-136, so under every scenario all
209 terms stay below 1e-5 and none changes status. **The enrichment conclusions are invariant; the
numbers are not.**

### 4.4 `scripts/pedigree/crossover_rate_f34.py` — not affected, and currently fail-silent

I ran it: `python3 scripts/pedigree/crossover_rate_f34.py` exits 0 and prints
`CHM13 PHR length (chm13.phrs.bed) = 0.000 Mb per haploid`. Reading `:101-112` and `:243-266`, that
value is a **diagnostic printout only** — the three denominator definitions use
`PHR_LENGTH_PER_HAPLOID_MB = 3.5` (pangenome, `:51`) and `N_ARMS_SIGNAL = 41`, so the reported rates
(0.762 / 1.524 / 0.263 per Mb per meiosis) do not depend on the BED at all. Two defects to fix
separately: `PHR_BED = REPO_ROOT / "chm13.phrs.bed"` (`:23`) points at a file that does not exist
(`ls chm13.phrs.bed` → No such file), and `phr_length_from_bed` returns **0.0** instead of raising
(`:103-104`). This violates our own fail-loud rule; the fix is `REPO_ROOT / "data/chm13.phrs.bed"`
plus a `raise FileNotFoundError`, and then the printed diagnostic becomes 6.015 Mb (vendored) or
6.840 Mb (published).

### 4.5 Fig 1a / 1b / 2b / 2c / 5A–D are clean

* Fig 1a reads `data/hprc25272-wf.CHM13.100kb-xm5-id098-l50000.tsv.gz`.
* Fig 1b's default input is `data/all-vs-all.1Mb.p95.id95.len.tsv` — the **published** table.
* Fig 2b/2c cell values are `data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv` + the Leiden table.
* Fig 5A is a hard-coded pedigree topology; 5B–D come from the `fig5_donor_recipient_ribbon_draft`
  and `fig5_extended_maternal_whole_genome_ribbons` tables.

---

## 5. RECOMMENDATION

**Regenerate the vendored BED from the published call and re-run the three deterministic consumers;
and state the call set in Methods either way.** Not the cheaper option of adding a caveat alone.

Grounds, in order of weight:

1. **The defect is already visible to a referee.** `arm_inclusion_audit.tsv` is a committed
   submission artefact that prints 37 next to 38 and names chr13_p (4.2). A referee tracing a number
   back to that file finds the inconsistency without any help from us.
2. **Cost is minutes, not days.** The conversion is one `awk` over a table already in the repo
   (giving `chr, start0, end0, chrs_involved` filtered to bare chromosome names); the two deterministic
   consumers then rerun in **≈14 s** (`build_inputs.py`, measured) and ≈1 s
   (`make_fig5_pedigree_supporting_material.py`). No moosefs, no Slurm, no guix.
3. **One consequence is qualitative, not cosmetic.** The number of CHM13 PHRs sitting at the 500 kb
   flank ceiling goes **2/37 (5.4 %) → 6/38 (15.8 %)**, a three-fold understatement of the
   right-censoring that Referee #1 major point 8 makes "the most important correction". chr14_p alone
   moves from 370 kb to 500 kb and chr15_p from 405 kb to 500 kb — i.e. the vendored BED hides two
   saturated arms.
4. **The scientific conclusions survive**, so this is a hygiene fix, not a reopening: the enrichment
   term statuses are invariant (4.3), and every locus in the published set is a superset-or-shift of
   the vendored set (0 loci lost, 2 gained, 1.2 §1.5).
5. **Do NOT regenerate alone, though** — the 99,999-permutation burden nulls inside `COPY_engine.py`
   sample spatially against the PHR intervals, so a new BED invalidates the frozen run identity
   (`results/v1/...`, 18.08 GiB of transient products). Regenerating the *inputs* without re-running
   the *nulls* would be worse than either extreme, because the observed counts and the null would then
   come from different BEDs. So the sequence must be: swap BED → rebuild `analysis_ready` → re-run
   the burden nulls on a moosefs worker → re-emit ED2. If that full chain cannot be scheduled, then
   the fallback below is mandatory, not optional.

### 5.1 Minimum acceptable Methods statement (if regeneration cannot be scheduled)

> PHR calls were made once, on the 18,826 retained 500 kb flanks, giving 15,668 PHRs
> (`data/all-vs-all.1Mb.p95.id95.len.tsv`). Gene-copy assignment, the ontology enrichments
> (Extended Data Fig. 2), the pedigree support table and the arm-level similarity matrix
> behind Fig. 2 use the CHM13 interval set derived from that call. Figure panels were
> generated from the deposited tables named in the Data availability statement, and each
> such table names the call it derives from in its header.

That last clause is the part that actually earns the referee's trust: make every deposited table
carry its provenance in-band.

### 5.2 Two follow-ups I recommend regardless of the choice above

* Fix `scripts/pedigree/crossover_rate_f34.py:23` and make `phr_length_from_bed` raise on a missing
  file (4.4). It currently prints `0.000 Mb` and exits 0.
* Add a CI gate that recomputes `data/chm13.phrs.bed` from
  `data/all-vs-all.1Mb.p95.id95.len.tsv` and byte-compares, the way
  `build_inputs.py --check` byte-compares `analysis_ready`. That is the structural fix: this class of
  drift cannot recur silently if the vendored derivative is checked against its own source.

---

## 6. Residual risks

1. **I did not re-run the permutation nulls**, so I cannot report the effect of the BED swap on the
   null distributions — only on the observed counts and the exact hypergeometric p's (4.3). The
   sensitivity bound there assumes the null is unchanged, which is exactly what would be false after a
   `COPY_engine` rerun.
2. **`K`'s new value is a bound, not a measurement** (187, 188 or 189), because the term maps that
   decide eligibility are not vendored in this worktree.
3. **The vendored BED mixes coordinate conventions** (`start = st0+rs`, end sometimes `st0+re` and
   sometimes `st0+re-1`): 22 of 37 rows are ±1 against the published intervals. Any consumer that
   intersects the BED with 1-based feature coordinates inherits a possible 1 bp boundary error. I did
   not chase this; it is pre-existing and independent of the call-set drift.
4. **I read `/moosefs` ~20 times with no hang and no I/O error**, so nothing here is degraded by a
   partial read; but the `*_global_test.tsv` values I quote for CHM13 and for the 10/20/50 kb sweep
   exist only on moosefs, not in the repo. Only the 20 kb HG002 pair is vendored (`data/hg002_hic_*`),
   which I verified matches moosefs exactly.
5. **No LaTeX binary on this box**, so the three patch files are balance-checked but not compiled.

---

## 7. Where each deliverable lives (final)

| path | contents | verified how |
|---|---|---|
| `response/notes/L5_data_reconciliation.md` | this file: the re-derivation, the per-consumer consequence table, the recommendation, the residual risks | every number from a command in this file's steps |
| `response/patches/flanking_control.tex` | Results replacement for `paper.tex:291–296` (both metrics, 50 kb + 20 kb, the metric-disagreement mechanism, the multi-mapping dilution argument) + optional Methods replacement for `:764–767` | `$` parity even; brace depth 0; macros limited to `\emph`, `\rho`, `\times` |
| `response/patches/pedigree_pattern_categories.tex` | one sentence each for `sandwich_same_hap` and `complex`, counts 115 and 1 of 494 | definitions read from the classifier source, not paraphrased; census cross-checked in two in-tree files |
| `response/patches/mouse_satellite_caveat.tex` | Results replacement `:306–312`, summary-paragraph replacement `:81–83`, Fig 4c caption replacement `:545–547`, plus a shorter fallback | satellite fractions recomputed from the JSON; three plausible-but-false framings explicitly refused |

Nothing outside `response/` was written. `git status --porcelain | grep -v '^??'` → 0 lines, i.e. **no
tracked file modified**, and no git commands were run.
