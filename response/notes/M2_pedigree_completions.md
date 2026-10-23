# R1-M2 (3 completions) — genome-wide identity margins, permuted parent-child FP rate, contig check

Status: three completions delivered as artifacts + measurements. Written 2026-10-16 by
`deepseek-flash-writer`. All numbers below were produced by the commands shown, on the
repo-local inputs named. Nothing is estimated or inferred.

Authority followed: `PEER_REVIEW.md` (real review) and `response/POINT_BY_POINT.md` §R1-M2
(`TODO:` "genome-wide distribution of identity margins with the candidates located on it" /
"permuted parent–child false-positive rate" / "contig-breakpoint check for the retained tract").
`paper_prep/synthesis/SIMULATED_PEER_REVIEW_v1.md` was **not** consulted.

Referee #1 verbatim (PEER_REVIEW.md): *"We ask for the comparison against both parental
haplotypes, a genome-wide distribution of identity margins locating the candidates, a
false-positive estimate from permuted parent-child pairs, and confirmation that each tract lies
in an unbroken contig."* The both-parental-haplotype comparison is already done (POINT_BY_POINT
R1-M2 "What we did"); this note covers the other three.

---

## 0. Inputs (all repo-local, no moosefs read for items 1 and 2)

| file | what it is |
|---|---|
| `data/fig5_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz` | IMPG class-winner table (layer provenance per `response/01_evidence_log_pedigree.txt`: sweepga 0.1.1 FastGA f32 10:10 → IMPG 0.4.1), 352,940 rows, 177,998 distinct 2 kb windows |
| `data/fig5_PAN027pat_vs_PAN011_joint.query.fa.fai` | recipient (query) contigs — PAN027 paternal |
| `data/fig5_PAN027pat_vs_PAN011_joint.target.fa.fai` | donor/competitor (target) contigs — PAN011 both haplotypes |
| `/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies_v11/PAN027_pat_HiFi_element_final_pat.v1.1.raw.fa.gz.fai` | assembly truth for item 3 (read fine, no timeout) |

Query = `PAN027#2#chrN`. Independent confirmation that `PAN027#2` is the **paternal** (SELF.p)
copy, not a haplotype index coincidence:

* `fig5_PAN027pat_...` names it `pat`;
* every query contig length is **identical** to the `PAN027.chrN.paternal` contig in the
  assembly `.fai` (all 23, zero mismatches) while `PAN027.chr9.maternal` is 129,071,589 bp vs
  the query's 136,209,149 bp — i.e. it is not the maternal chr9.

### Identity-column semantics (important for every number below)

`estimated.identity` in this table **is** the IMPG window-level dice/coverage identity, not a
base-level % identity recomputed from raw alignment:

```
python3 scripts/pedigree/identity_margin_distribution.py   # audits this itself
# TSV row: summary  max_abs_delta_estimated_identity_minus_dice  1e-07
```

(max |estimated.identity − dice.similarity| over all 352,940 rows = 1e-07, floating-point only.)
Consequence: a "1% margin" in this metric is a 1% difference in coverage-weighted estimated
identity, not 1% base divergence. Every statement below inherits that caveat.

---

## 1. Genome-wide identity-margin distribution

**Definition used** (POINT_BY_POINT R1-M2's own logic, and the referee's):

```
margin(window) = identity(best interchromosomal parental match)
               - identity(best same-chromosome/homologous parental match)
```

Both terms are the per-window class winners already used by the Fig. 5 pipeline
(`winner_class` ∈ {`interchrom`, `same_chrom`}; exactly one row per window per class, which the
script asserts).

**Script:** `scripts/pedigree/identity_margin_distribution.py`
**Artifact:** `scripts/pedigree/identity_margin_distribution.tsv` (96 data rows; tidy 6-column
TSV `record_type/label/value/value2/value3/note`, sections `summary`, `quantile`, `hist`,
`candidate_run`, `candidate_window`, `candidate_context` — no genome-wide per-window dump is
written).

### Command and verbatim stdout

```
$ cd /home/erikg/phrs && python3 scripts/pedigree/identity_margin_distribution.py
windows_total	177998
margin_defined	174942
interchrom_only	3052
margin>0	1778 (1.0163%)
median	-0.339722	p95	-0.004251	p99	0.000250
qualifying_runs	156
candidate	PAN027#2#chr9	136168000	136188000	20000 bp	donor=chr3	mean_id=0.998775	rank=61/156
candidate	PAN027#2#chr5	182052000	182080000	28000 bp	donor=chr1	mean_id=0.998267	rank=41/156
wrote	scripts/pedigree/identity_margin_distribution.tsv
```

Wall: `real 0m2.096s / user 0m1.857s / sys 0m0.238s` (single core, repo-local gz input,
nothing else of mine running on the box; treat as indicative, box load not otherwise controlled).

### Window accounting

| class | windows |
|---|---|
| total windows present in the table | 177,998 |
| both classes present → **margin defined** | 174,942 |
| `interchrom` only (no homologous parental match in the window) | 3,052 |
| `same_chrom` only (no non-homologous parental match) | 4 |
| margin exactly 0 (paralogue and homologue tie) | 839 |

### Distribution (n = 174,942 margin-defined windows)

Quantiles (nearest-rank, from the TSV `quantile` rows):

| q | margin | q | margin |
|---|---|---|---|
| p0 | −0.999001 | p90 | −0.017985 |
| p1 | −0.974384 | p95 | −0.004251 |
| p5 | −0.894439 | p97.5 | −0.001250 |
| p10 | −0.822969 | **p99** | **+0.000250** |
| p25 | −0.637980 | p99.5 | +0.005500 |
| **p50** | **−0.339722** | p99.9 | +0.256410 |
| p75 | −0.081459 | p100 | +0.993018 |

mean −0.376486, sd 0.304060, min −0.999001, max 0.993018.

Threshold counts (TSV `summary` rows):

| margin > | windows | % of margin-defined |
|---|---|---|
| 0 | 1,778 | 1.0163 |
| 0.001 | 1,466 | 0.8380 |
| 0.002 | 1,200 | 0.6859 |
| 0.005 | 911 | 0.5207 |
| 0.01 | 716 | 0.4093 |
| 0.02 | 572 | 0.3270 |
| 0.05 | 391 | 0.2235 |
| 0.10 | 310 | 0.1772 |

Histogram: 26 bins in the TSV (`hist` rows), from `[-1.0000,-0.9000)` (9,151 windows, 5.2309%)
through `=0` (839, 0.4796%) to `(0.5000,1.0000)` (98, 0.0560%). **The distribution is centred
strongly negative**: the homologous parental copy is the better match in the overwhelming
majority of the genome, and a positive margin is by itself a top-1% event.

### Where the candidates fall

Candidate windows are **derived, not hard-coded**: on each marker query chromosome, the longest
run inside the terminal 100 kb satisfying the scan's own run filter (see §2). This independently
reproduces the coordinates already reported by hand in
`response/01_evidence_log_pedigree.txt` (chr5q 182.052–182.080 Mb, 28 kb; chr9q telomere-first
~13–57 kb block), which is a useful cross-check on the two artifacts.

Per-window margins and their genome-wide percentiles (`candidate_window` rows in the TSV;
percentile = % of the 174,942 margin-defined windows strictly below that margin):

**(a) Retained chr9q←chr3q tract** — `PAN027#2#chr9:136,168,000–136,188,000`, 20,000 bp,
10 windows, single donor `PAN011#joint#h2_chr3`, mean identity 0.998775, no homolog-absent
window. Margins 0.000250–0.021347 (mean 0.006840) → percentiles **99.010 / 99.536 / 99.681**
(min / mean / max). Length rank **61 of 156** qualifying runs genome-wide.

Broader chr3-donor census in the terminal 60 kb (one-off command, recorded because the strict
run filter breaks the block):

```
18 windows in chr9 terminal 60 kb have best interchrom donor chr3 (all PAN011#joint#h2_chr3):
first at 136,152,000, last at 136,204,000–136,206,000, so the covered interval is
136,152,000–136,206,000 (54,000 bp). 0 homolog-absent windows.
```

i.e. the class-winner layer sees the same ~54 kb chr3-derived terminal block the FastGA
analysis reports (evidence log rev 12: "100% of the terminal 60 kb covered at 99.79–99.84%");
the 20 kb figure is only what survives *uninterrupted* same-donor contiguity, because at
136,154 / 136,158 / 136,160 / 136,188 / 136,192 / 136,194 / 136,196 / 136,200 / 136,202 a
**different paralogue** (chr7 / chr11 / chr19) wins by ≤0.6% inside the same near-identical
family. That competition is the substance of referee sub-point (3), not an artifact of the
filter.

**(b) Withdrawn chr5q/chr1p candidate** — `PAN027#2#chr5:182,052,000–182,080,000`, 28,000 bp,
14 windows, single donor `PAN011#joint#h2_chr1`, mean identity 0.998267.
Margins 0.000500–0.007944 (mean 0.003196) → percentiles **99.080 / 99.411 / 99.556**.
Length rank **41 of 156**.

### What the margins actually say (this is a confirmation of the referee, not a rebuttal)

The two candidates are **not separable by margin**:

| | min margin | mean margin | max margin |
|---|---|---|---|
| retained chr9q←chr3q | 0.000250 | 0.006840 | 0.021347 |
| withdrawn chr5q/chr1p | 0.000500 | 0.003196 | 0.007944 |

and the retained tract is not even in the top tail of the genome:

```
558 of the 1,778 positive-margin windows genome-wide have margin > 0.021347
777 have margin > 0.007944   (the chr5 candidate's maximum)
1609 have margin > 0.000500   (all 1778 are > 0 by construction)
```

**Therefore: the identity margin does not discriminate the retained from the withdrawn call.**
The point-by-point's existing concession ("margins are ≤1%") is correct and is now quantified
genome-wide. Discrimination for the chr9q case has to come from the structural evidence
(module present on the father's chr9 h2, absent from his chr9 h1; coverage exclusion of the
homologue), exactly as `response/01_evidence_log_pedigree.txt` rev 11/12 argues — not from the
margin. This note should be cited for the distribution, not as support for a margin-based call.

---

## 2. Permuted parent–child false-positive rate

**Script:** `scripts/pedigree/permuted_parent_child_fp.py`
**Artifact:** `scripts/pedigree/permuted_parent_child_fp.json` (contains the full null
distribution summaries, the observed run list, per-chromosome marginals and a `limitations`
block).

**Scan rule applied to both observed and permuted data** (identical code path, imported from
`identity_margin_distribution.py`, single source of truth):

* eligible window = best interchromosomal parental identity ≥ 0.95 **and** (no same-chromosome
  parental match **or** interchromosomal identity > same-chromosome identity);
* tract = adjacent touching eligible windows, single donor **chromosome**, ≥ 10 kb (≥ 5
  windows), mean identity ≥ 0.95.

### Command and verbatim stdout

```
$ cd /home/erikg/phrs && time python3 scripts/pedigree/permuted_parent_child_fp.py --reps 10000
observed qualifying runs	156	(3886000 bp)
candidate PAN027#2#chr9	136168000-136188000	20000 bp	donor=chr3	mean_id=0.998775	runs_overlapping=1
   null A runs_overlapping>=obs 0/10000 p=<1/10000 | eligible_inside obs=10 null_mean=0.185 p=<1/10000
   null B runs_overlapping>=obs 908/10000 p=0.0908 (null_mean=0.091 null_max=2)
candidate PAN027#2#chr5	182052000-182080000	28000 bp	donor=chr1	mean_id=0.998267	runs_overlapping=1
   null A runs_overlapping>=obs 0/10000 p=<1/10000 | eligible_inside obs=14 null_mean=0.081 p=<1/10000
   null B runs_overlapping>=obs 3207/10000 p=0.3207 (null_mean=0.340 null_max=3)
null A total runs: mean=0.0221 sd=0.1484 max=2  reps>=obs=0/10000 p=<1/10000
null B total runs: mean=48.4953 sd=5.1243 max=69  reps>=obs=0/10000 p=<1/10000
wrote	scripts/pedigree/permuted_parent_child_fp.json
```
Wall: `real 0m50.655s / user 0m49.917s / sys 0m0.737s` (single core; my only running job).

Determinism checked: two independent `--reps 200` runs produce byte-identical JSON
(`cmp` clean), seeds are `--seed 20261016` (null A) and `--seed 20261017` (null B).

### The transmitted-copy permutation the referee suggested is not available from this input

`permutation_audit` in the JSON:

```
same_chrom_winner_haplotype_counts = {"PAN011#joint#h1": 68416, "PAN011#joint#h2": 106530}
transmitted_copy_permutation_feasible = false
```

The class-winner table retains only the winner *within* the `same_chrom` class, i.e. it is
already the per-window maximum over the father's two haplotypes (both appear, 68,416 / 106,530
windows). The per-haplotype identities needed to swap which haplotype is "the transmitted copy"
are not in this file. We therefore use the alternative the referee also allows — shuffle window
labels preserving per-chromosome (per-arm) window counts.

### null A — permute the whole per-window record within each query chromosome

Per chromosome this preserves: number of scanned windows, number of eligible windows, and the
joint (donor, identity) distribution of eligible windows. It destroys spatial order, hence both
the clustering of candidate windows and single-donor contiguity.

Implemented exactly rather than approximately: a uniform permutation of the N_c records over N_c
window positions puts the K_c eligible records on a uniformly random *ordered* K_c-subset of the
positions, so positions are sampled without replacement and the shuffled eligible records are
assigned to them.

| statistic | observed | null A (10,000 reps) | p (≥ obs) |
|---|---|---|---|
| total qualifying runs genome-wide | **156** (3,886,000 bp) | mean 0.0221, sd 0.1484, median 0, p95 0, p99 1, **max 2** | **< 1e-4** (0/10,000) |
| eligible windows inside chr9 candidate interval (10 windows) | **10** | mean 0.1854, sd 0.4281, max 3 | **< 1e-4** (0/10,000) |
| eligible windows inside chr5 candidate interval (14 windows) | **14** | mean 0.0812, sd 0.2853, max 2 | **< 1e-4** (0/10,000) |
| qualifying runs overlapping chr9 candidate | **1** | mean 0.0000, max 0 | **< 1e-4** (0/10,000) |
| qualifying runs overlapping chr5 candidate | **1** | mean 0.0000, max 0 | **< 1e-4** (0/10,000) |

Reading: given this many windows per chromosome in which the paralogue beats the homologue, and
given their identity/donor marginals, runs of the required kind essentially never form. The
observed 156 runs and the local clustering at both candidate intervals are far outside null A.

### null B — keep eligible/not status and identity, shuffle only the donor label per chromosome

This conditions on *where* the paralogous winners are and tests only the remaining claim, that a
run shares one donor more often than chance.

| statistic | observed | null B | p (≥ obs) |
|---|---|---|---|
| total qualifying runs genome-wide | **156** | mean 48.4953, sd 5.1243, min 32, **max 69** | **< 1e-4** (0/10,000) |
| qualifying runs overlapping chr9 candidate | **1** | mean 0.091, max 2 | **0.0908** (908/10,000) |
| qualifying runs overlapping chr5 candidate | **1** | mean 0.340, max 3 | **0.3207** (3,207/10,000) |

Reading — the honest split:

* **Genome-wide, donor contiguity is highly non-random** (156 vs 48.5 ± 5.1, p < 1e-4).
* **But within the two candidate intervals it is not.** Once you condition on the eligible
  windows being there, the fact that they share one donor chromosome is unremarkable: p = 0.09
  for the retained chr9q tract and p = 0.32 for the withdrawn chr5q/chr1p candidate. The
  candidate-level improbability lives in the *clustering of paralogous winners*, not in the
  donor label. This is consistent with both candidates sitting inside a single near-identical
  paralogous family where any one member could have won each window.

### Context for the observed run population (one-off commands, recorded here)

```
qualifying runs per query chromosome:
  chr1 14, chr13 67, chr14 1, chr15 10, chr18 7, chr21 19, chr22 30, chr5 3, chr9 3, chrX 2
runs wholly inside the first 10 Mb of an acrocentric p arm (chr13/14/15/21/22): 127
runs wholly inside the first 1 Mb of chrX (PAR):                                 2
runs containing >=1 homolog-absent window: 111   (all windows homolog-absent: 76)
runs >= 28 kb: 41        runs >= 20 kb: 61
longest runs: 142, 136, 104, 102, 102, 96, 90, 78, 76, 74 kb
```

The run population is dominated by the acrocentric short arms (127/156) and PAR (2), i.e. the
expected positive controls; the two pedigree candidates are in the residual 27.
`homolog_absent_windows = 71` for the longest (chr15p↔chr13) run: the father's assembly has **no
chr15-vs-chr15 alignment at all** in those windows, which for an acrocentric short arm means the
paternal p-arm copy is not on his `chr15` contig — that is an assembly-representation statement,
not evidence of exchange, and it is why the "homolog absent" arm of the eligibility rule must not
be read as a stronger claim than it is.

### Sensitivity: the margin requirement is doing work

One-off variant with the `margin > 0` requirement dropped (eligible = interchrom identity ≥ 0.95
regardless), **not part of the emitted artifacts**:

```
SENSITIVITY (margin requirement dropped): qualifying runs = 1244, total bp = 34,276,000
```

156 → 1,244 runs. The scan's discriminating power comes from the margin requirement, exactly as
the referee suspected; without it the filter is essentially "any 10 kb of ≥95% inter-chromosomal
match" and returns 34 Mb of the genome.

---

## 3. Contig-breakpoint check for the retained chr9q/chr3q tract

Requirement: confirm the tract lies inside an unbroken contig, using
`data/fig5_PAN027pat_vs_PAN011_joint.query.fa.fai` plus the assembly `.fai`. The moosefs
directory was readable (no I/O error, no timeout) — commands below were run with `timeout 60`.

### Exact contig and coordinates

| item | value |
|---|---|
| recipient contig (query/scan names) | `PAN027#2#chr9` |
| recipient contig (assembly name) | **`PAN027.chr9.paternal`** |
| recipient contig length | **136,209,149 bp** |
| retained tract (strict same-donor run) | **`PAN027#2#chr9:136,168,000–136,188,000`** (20,000 bp) |
| distance from tract end to contig end | **21,149 bp** |
| chr3-donor block span in terminal 60 kb | `136,152,000–136,206,000` (54,000 bp covered interval), 3,149 bp from contig end |
| donor contig | `PAN011#joint#h2_chr3` = **`PAN011.chr3.haplotype2`**, 205,578,940 bp (assembly `.fai`: 205,578,940 — identical) |

### Evidence 1 — the query fasta is the assembly chromosome contig, one record per chromosome

```
$ M=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131/assemblies_v11
$ timeout 60 cat $M/PAN027_pat_HiFi_element_final_pat.v1.1.raw.fa.gz.fai > /tmp/pan027_pat.fai
# 23 contigs, one per chromosome, named PAN027.chrN.paternal
$ awk '$1 ~ /chr9/' /tmp/pan027_pat.fai
PAN027.chr9.paternal	136209149	1572482401	60	61
```

Comparing all 23 contigs (query vs assembly): **23/23 lengths identical, 0 mismatches**, and the
byte-offset difference is a constant −7 per successive contig — exactly the length difference
between the query header `PAN027#2#chrN` (13 chars) and the assembly header `PAN027.chrN.paternal`
(20 chars). Same sequences, one record per chromosome, renamed only. **There is exactly one chr9
contig and therefore no contig boundary anywhere on chr9**, let alone inside the tract.

```
PAN027#2#chr1  PAN027.chr1.paternal  247769562 247769562 len_ok off_delta -7
...
PAN027#2#chr9  PAN027.chr9.paternal  136209149 136209149 len_ok off_delta -63
...
PAN027#2#chrX  PAN027.chrX.paternal  154091348 154091348 len_ok off_delta -161
mismatches: []
```

(The same check on the target side: all 23 `PAN011#joint#h2_chr*` contigs are length-identical to
`PAN011.chrN.haplotype2`, offset deltas ~3,079,396,803 — the size of the haplotype-1 block that
precedes them in the joint fasta. The donor contig is likewise unbroken.)

### Evidence 2 — no assembly gap (N-run) anywhere in the contig

```
$ timeout 120 samtools faidx $M/PAN027_pat_HiFi_element_final_pat.v1.1.raw.fa.gz \
      PAN027.chr9.paternal:136148000-136209149 > /tmp/tract9.fa      # exit 0
# 61,150 bp extracted (= 136209149-136148000+1), 0 N-runs
$ samtools faidx ... PAN027.chr9.paternal > /tmp/chr9_full.fa        # exit 0
chr9 contig length: 136209149
N runs total: 0
total N bp: 0
N runs inside 136148000-136209149: []
```

The whole `PAN027.chr9.paternal` contig contains **zero** N runs, and the 61,150 bp window that
contains the entire tract contains zero N runs. So both the scaffold-join question and the
gap question come out clean for the recipient.

### What said nothing

The assembly fasta headers carry only `>PAN027.chr9.paternal`; there is no join/stitch annotation
to read, so this check establishes gaplessness and one-contig-per-chromosome, not the absence of
a historical mis-join that the assembler did not record. That is the residual limit of an
`.fai`-based check and should be stated as such.

---

## What this does NOT establish

1. **The permutation null is not an orthogonal assay.** No permutation of labels can confirm or
   refute an inter-chromosomal exchange. It bounds the false-positive rate *under the null it
   enforces*: null A bounds it against "the window-level marginals could have produced this",
   null B against "donor identity is chance given where the winners are". Both nulls are
   compatible with exchanged sequence, with inherited paralogy, and with assembly/mapping
   artefacts; they cannot separate those three. An orthogonal assay (read-level support, an
   independent aligner, a population-level footprint) is still required and is not in this note.
2. **Null A is anti-conservative in a specific direction.** It destroys spatial autocorrelation,
   which is precisely the structure being tested, so its tiny p-values say "not explained by
   per-window marginals", nothing more. It also preserves the *absence* of windows (no alignment
   at all) rather than permuting it; those absences are non-random.
3. **Null B conditions on the observed eligible-window positions**, so it inherits whatever
   produced them (including the mapping competition inside the paralogous family). Its
   candidate-level p = 0.09 / 0.32 must not be read backwards as "the tract is 90% likely to be
   false" — it only says donor contiguity adds nothing beyond the clustering.
4. **`estimated.identity` is not base-level identity.** It equals `dice.similarity` to ≤1e-7: it
   is coverage-weighted window-level estimated identity from alignment blocks. Margins quoted
   here are in that metric. No base-level percent-identity is reported or implied.
5. **The margin metric does not discriminate the retained from the withdrawn call.** §1 shows the
   two candidates overlap and that ~557 genome-wide windows exceed the retained tract's maximum
   margin. This note therefore cannot be used to argue the chr9q call on its margin — only the
   structural/coverage evidence (evidence-log rev 11/12) does that.
6. **The 20 kb strict run is not the tract length.** The chr3-derived terminal block is 54 kb by
   donor census (18 windows, 136,152,000–136,206,000) and ~54–60 kb by the FastGA analysis in the
   evidence log. The 20 kb is only what survives uninterrupted single-donor contiguity under the
   run filter; do not quote 20 kb as the tract length.
7. **The contig check covers the recipient only in full.** The donor contig was checked at `.fai`
   level (length identity with the assembly) but not gap-censused; and neither check excludes a
   correctly-recorded-but-wrong assembly join (the headers carry no join annotation).
8. **"Homolog absent" is an assembly statement, not a biological one** (see §2: 111 of 156 runs
   contain such windows, and the acrocentric p-arm runs are entirely of that kind because the
   father's p-arm copy is not on his `chrN` contig).
9. **Single parent–child comparison only.** Everything here is PAN027-paternal (SELF.p) vs PAN011
   (father). The maternal comparison (PAN027 maternal vs PAN010) and the child-side comparisons
   are not included even though `data/fig5_PAN027mat_vs_PAN010_joint.*` is present; the withdrawn
   child-side calls would need the same treatment and are out of scope for this item.
10. **Nothing here re-derives the Fig. 5 ribbon calls.** The scripts read the class-winner table
    as given by the sweepga→IMPG layer; they do not re-run sweepga/wfmash and cannot detect an
    upstream aligner error.

---

## Reproduce

```
cd /home/erikg/phrs
python3 scripts/pedigree/identity_margin_distribution.py            # 2 s, writes the .tsv
python3 scripts/pedigree/permuted_parent_child_fp.py --reps 10000   # 51 s, writes the .json
```

Both are deterministic given their seeds and take only repo-local inputs; no scratch directory
is touched and no text is materialized beyond the two declared artifacts.
