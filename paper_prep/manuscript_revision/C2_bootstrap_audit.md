# C2 Bootstrap Audit: Character Bootstrap and the 0.1% q-arm Result

Date: 2026-06-17
Task: `manuscript-revision-c2`

## Question

Audit `scripts/cladistics/char_bootstrap_d_m9.R` and its cached involvement
table. Determine what a table row represents, whether the bootstrap ever uses
the 15,668 x 15,668 sequence-level similarity matrix, and whether the reported
0.1% support for the six-q-arm grouping reflects structural blindness of this
bootstrap input or genuine instability of the full sequence-level signal.

## Files inspected

- Script: `scripts/cladistics/char_bootstrap_d_m9.R`
- Script default PHR involvement cache:
  `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv`
- 1 Mb involvement cache / repo-local snapshot:
  `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv`
  and `data/all-vs-all.1Mb.p95.id95.len.tsv`
- Arm-level distance cache:
  `data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv` and
  `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`
- Secondary prior analysis note:
  `paper_prep/synthesis/ANALYSIS_D_M9.md`

I did not load or scan the full sequence-level similarity matrix or RDS object.
All checks were limited to script text, `ls`, headers, row counts, and `head`.

## Code path

`char_bootstrap_d_m9.R` identifies the analysis as "character-level" PHR
resampling in its header: it says it resamples PHRs, recomputes an arm-level
Jaccard matrix per replicate from cached PHR-level cross-chromosome
involvement, and then runs `ape::nj()` plus UPGMA
(`scripts/cladistics/char_bootstrap_d_m9.R:7` to
`scripts/cladistics/char_bootstrap_d_m9.R:12`).

The input path is parameterized as `--phr-tsv`. The default in the script is:

```text
/home/guarracino/Desktop/Garrison/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv
```

with the canonical moosefs path documented in the header as:

```text
/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv
```

Relevant code lines:

- `PHR_TSV` is set from `--phr-tsv` with the default above
  (`scripts/cladistics/char_bootstrap_d_m9.R:47` to
  `scripts/cladistics/char_bootstrap_d_m9.R:48`).
- The table is loaded with `read.table(PHR_TSV, header = TRUE, sep = "\t", ...)`
  (`scripts/cladistics/char_bootstrap_d_m9.R:66` to
  `scripts/cladistics/char_bootstrap_d_m9.R:68`).
- Rows with empty or dot `chrs_involved` are removed
  (`scripts/cladistics/char_bootstrap_d_m9.R:71` to
  `scripts/cladistics/char_bootstrap_d_m9.R:75`).
- `arm == "parm"` and `arm == "qarm"` are converted to `p` and `q`, then
  combined with `self_chr` into `self_arm`
  (`scripts/cladistics/char_bootstrap_d_m9.R:77` to
  `scripts/cladistics/char_bootstrap_d_m9.R:83`).
- The reference arm labels are read from a small arm-level matrix, not from the
  sequence-level matrix (`scripts/cladistics/char_bootstrap_d_m9.R:98` to
  `scripts/cladistics/char_bootstrap_d_m9.R:104`).
- `make_phr_arm_long()` creates long-form `(phr, arm)` membership pairs by
  splitting `chrs_involved`, appending the same p/q suffix as the row's self
  arm, adding the self arm, and filtering to reference arms
  (`scripts/cladistics/char_bootstrap_d_m9.R:112` to
  `scripts/cladistics/char_bootstrap_d_m9.R:128`).
- A dense `N_phr x N_arm` matrix `M` is created only for PHR-by-arm membership,
  with comments stating this is `15k x 42`, about 5 MB as doubles
  (`scripts/cladistics/char_bootstrap_d_m9.R:135` to
  `scripts/cladistics/char_bootstrap_d_m9.R:139`).
- The per-replicate distance function builds a `N_arm x N_arm` count matrix `C`,
  row-normalizes it, and computes an arm-arm intersection distance
  (`scripts/cladistics/char_bootstrap_d_m9.R:149` to
  `scripts/cladistics/char_bootstrap_d_m9.R:191`).
- The bootstrap loop samples row indices with replacement, tabulates row
  multiplicities, recomputes the arm-level distance, and rebuilds trees
  (`scripts/cladistics/char_bootstrap_d_m9.R:308` to
  `scripts/cladistics/char_bootstrap_d_m9.R:316`).

The script therefore has two data scales:

1. `N_phr` rows from the involvement table.
2. `N_arm` reference arms, read from a small arm-level distance matrix.

There is no code path that reads `similarity.tsv.gz`,
`hprcv2.1Mb.subtelo.dist_matrix.rds`, or any 15,668 x 15,668 pairwise matrix.

## Cache paths and sizes

Lightweight metadata sampled on 2026-06-17:

```text
/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv
  size: 2.4M
  lines: 18,226 including header
  columns: 6
  data rows: 18,225
  rows with non-dot/non-empty chrs_involved by awk: 15,088

/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv
  size: 4.5M
  lines: 18,827 including header
  columns: 7
  data rows: 18,826
  rows with non-dot/non-empty chrs_involved by awk: 15,668

data/all-vs-all.1Mb.p95.id95.len.tsv
  size: 4.5M
  lines: 18,827 including header
  columns: 7
  data rows: 18,826
  rows with non-dot/non-empty chrs_involved by awk: 15,668

data/hprcv2.1Mb.subtelo.arm_dist_matrix.tsv
  size: 30K
  lines: 42 including header
  matrix labels: 41 arms plus row-name column
```

The 1 Mb table is the table matching the manuscript count of 15,668
signal-bearing rows. The script default, however, points to the shorter
non-`1Mb` path unless `--phr-tsv` is supplied. The prior D-M9 note reports the
bootstrap run as using `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.p95.id95.len.tsv`
and 15,089 signal-bearing PHR flanks
(`paper_prep/synthesis/ANALYSIS_D_M9.md:40` to
`paper_prep/synthesis/ANALYSIS_D_M9.md:47`). My direct `awk` count on that same
file found 15,088 rows with column 6 not empty or `"."`; I did not rerun R to
resolve this one-row discrepancy.

## Row schema and representative row

Script-default table header:

```text
seq    arm    self_chr    region_start    region_end    chrs_involved
```

Representative row from the script-default moosefs table:

```text
seq:            CHM13#0#chr10:134254995-134754994_chr10_qarm
arm:            qarm
self_chr:       chr10
region_start:   320000
region_end:     500000
chrs_involved:  chr1,chr2,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr13,chr16,chr17,chr18,chr19,chr20,chr21,chr22
```

The repo-local 1 Mb snapshot adds one column:

```text
seq    arm    self_chr    region_start    region_end    chrs_involved    arms_involved
```

Representative row from `data/all-vs-all.1Mb.p95.id95.len.tsv`:

```text
seq:            CHM13#0#chr10:134254995-134754994_chr10_qarm
arm:            qarm
self_chr:       chr10:134254995-134754994
region_start:   320000
region_end:     500000
chrs_involved:  chr1,chr2,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr13,chr14,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr4:191168094-191668093,chr4:191882350-192382349,chr4:193072741-193572740
arms_involved:  chr1p,chr1q,chr2q,chr4q,chr5q,chr6p,chr6q,chr7p,chr8p,chr9p,chr9q,chr10q,chr11p,chr13q,chr14p,chr16p,chr16q,chr17p,chr17q,chr18p,chr19p,chr19q,chr20p,chr20q,chr21q,chr22q
```

A row is not a pairwise PHR-PHR comparison and not one cell of a sequence
similarity matrix. A row represents one telomere-anchored PHR/flank record:
the row's own sequence/flank ID, its p/q side, its anchor chromosome or anchor
interval, the interval within the flank that carries inter-chromosomal signal,
and a collapsed list of chromosomes involved in that row's high-identity
inter-chromosomal sharing. The 1 Mb snapshot also includes an `arms_involved`
summary, but the audited script does not use that column. It uses only
`chrs_involved` plus the row's own `arm` to infer same-p/q arm labels.

## Does this bootstrap touch the 15,668 x 15,668 matrix?

No. The script never reads the full sequence-level matrix.

The evidence is direct:

- The only large-ish PHR input is `PHR_TSV`, read at
  `scripts/cladistics/char_bootstrap_d_m9.R:66` to
  `scripts/cladistics/char_bootstrap_d_m9.R:68`.
- The only distance matrix input is `REF_DIST`, read as a small arm-level
  table at `scripts/cladistics/char_bootstrap_d_m9.R:98` to
  `scripts/cladistics/char_bootstrap_d_m9.R:104`.
- The bootstrap resamples integer row indices from `1:N_phr`
  (`scripts/cladistics/char_bootstrap_d_m9.R:308` to
  `scripts/cladistics/char_bootstrap_d_m9.R:310`).
- Each replicate recomputes a small `N_arm x N_arm` distance matrix from row
  multiplicities and the precomputed PHR-to-arm long table
  (`scripts/cladistics/char_bootstrap_d_m9.R:311` to
  `scripts/cladistics/char_bootstrap_d_m9.R:316`).

The `15,668 x 15,668` sequence-level Jaccard matrix is documented elsewhere as
the output of `pggb + odgi similarity` and as the downstream sequence-level
community input, but it is not part of this character-bootstrap script.
`paper_prep/manuscript_revision/00_inventory.md` lists the heavy matrix paths
and explicitly warns against loading them on the head node; this audit did not
use them.

## Interpretation of the 0.1% q-arm support

The prior D-M9 summary reports `TIGHT_q` support as 0.1% for NJ and 0.0% for
UPGMA (`paper_prep/synthesis/ANALYSIS_D_M9.md:19` to
`paper_prep/synthesis/ANALYSIS_D_M9.md:27`). It also describes the method as
resampling PHR rows with replacement and rebuilding the arm-level tree for each
replicate (`paper_prep/synthesis/ANALYSIS_D_M9.md:38` to
`paper_prep/synthesis/ANALYSIS_D_M9.md:53`).

The code audit changes how that result should be read. The 0.1% support is
support under a chromosome-level involvement surrogate, not support under the
full 15,668 x 15,668 PGGB/ODGI sequence-level Jaccard matrix. The surrogate
collapses each PHR row to:

```text
self arm + {chromosome labels in chrs_involved, projected onto the same p/q side}
```

as implemented in `make_phr_arm_long()`
(`scripts/cladistics/char_bootstrap_d_m9.R:112` to
`scripts/cladistics/char_bootstrap_d_m9.R:128`). That representation does not
encode pairwise PHR-PHR similarity, per-base shared-segment counts, or
within-q-arm sequence distinctions from the full matrix. It therefore has
structural blindness to within-q splits that depend on information absent from
`chrs_involved`.

This supports the following narrow conclusion:

- The 0.1% q-arm value is genuine instability of the audited surrogate tree
  under row resampling.
- It is not evidence, by itself, that the six-q-arm grouping is genuinely
  unstable in the full 15,668 x 15,668 sequence-level similarity matrix.
- The script cannot adjudicate full-matrix q-arm stability because it never
  reads or resamples that matrix.
- The safe manuscript-facing phrasing is that the D-M9 character bootstrap was
  performed on a collapsed per-PHR involvement cache and exposes low support
  for the six-q-arm split in that surrogate. It should not be presented as a
  full sequence-matrix bootstrap or as a mechanistic explanation for why the
  q-arm grouping is biologically unstable.

## Bottom line

One row in the cache is one PHR/flank involvement record, not one pairwise
similarity observation. The bootstrap resamples those rows, reconstructs
PHR-to-arm memberships from `chrs_involved`, and rebuilds a small arm-level
surrogate distance matrix. It does not touch the 15,668 x 15,668 sequence-level
matrix. Therefore the reported 0.1% q-arm support should be treated as a
limitation/result of the collapsed involvement-table surrogate, especially its
loss of within-q resolution, rather than as a demonstrated instability of the
full sequence-level q-arm structure.
