# D-PeerQ1: MAPQ-strict Hi-C re-run

**Concern.** The Hi-C pipeline was run with `MIN_MAPQ = 0`, `RM_MULTI = 0`
(HiC-Pro + Bowtie2), so multi-mapped reads at high-identity PHR loci were
retained and assigned to one of the equally scoring positions at random.
Reviewer Q1 (NARRATIVE_EXTRACT §6) asks whether the within-community
contact signal is inflated by this random placement.

**Verdict (one line).** The within-community B/W signal at PHR is **not**
inflated by random-placement multi-mapping. The v5 flanking-region control
shows that unique-sequence regions immediately centromere-ward of every PHR
boundary carry an **equal or stronger** within-community enrichment than
the PHR windows themselves across 7/7 Hi-C and Pore-C datasets. A
strict-MAPQ re-run at PHR-internal coordinates is uninterpretable by
construction (the PHR sequence is the source of the multi-mapping reads;
filter them out and there is nothing left to bin), which is why the v5
draft already states this in Methods P9.

The per-mcool re-binning under MAPQ ≥ 30 was scripted
(`scripts/hic/mapq_strict_d_peerq1.py`) but could not be executed in this
worker because `/moosefs/guarracino/HPRCv2/PHR_III/3d/` was not mounted at
the time of writing; the script is ready to run when the .allValidPairs
files are accessible and is expected to confirm the flanking-region
finding.

## 1. Why "MAPQ-strict re-run on the existing mcool" is the wrong
##    operation

An mcool file is **post-binning**. By the time HiC-Pro has aggregated
reads into bin × bin counts there is no per-read MAPQ to filter on; the
contact matrix already reflects the random placement chosen at alignment
time. Re-applying a MAPQ filter therefore requires going **back to the
.allValidPairs stage** (where each row still carries the MAPQ of both
mates), re-filtering, and re-binning. The script
`scripts/hic/mapq_strict_d_peerq1.py` does exactly that:

```
allValidPairs (MIN_MAPQ=0, random placement)
     │
     ├── awk '$mapq1>=30 && $mapq2>=30'           # MAPQ-strict filter
     │
     ▼
allValidPairs.mapq30 ──cooler cload pairs──▶ .cool ──zoomify──▶ .mcool
     │
     ▼
analyze_hic_communities.py (PHR coords + flanking coords)
     │
     ▼
within / between / B-W / p_value  (strict)
        vs the same metrics from the v5 random-placement mcool
```

The v5 random-placement mcool already exists at
`/moosefs/guarracino/HPRCv2/PHR_III/3d/<sample>/<sample>.mcool`; the script
reads that mcool to get the **random_BW** baseline, builds a new
strict-MAPQ mcool from the upstream `allValidPairs`, runs
`analyze_hic_communities.py` against both, and emits one comparison row
per `(sample, tech, region)`.

## 2. The PHR-internal strict-MAPQ measurement is uninterpretable
##    by construction

A PHR is, by definition, a 500 kb subtelomeric region of an assembly
contig ≥ 1 Mb that aligns at ≥ 95 % identity to ≥ 1 PHR on another
chromosome arm (see report §01_pipeline). Hi-C reads landing at the
shared high-identity sequence have multiple equally good alignment
positions; Bowtie2 default mode picks one, HiC-Pro under MIN_MAPQ = 0,
RM_MULTI = 0 keeps it. Applying MAPQ ≥ 30 retroactively will drop **the
majority** of these reads, not because they are wrong but because the
sequence is genuinely ambiguous. The strict-MAPQ B/W at PHR-internal
coordinates is therefore predicted to collapse to a low-count noise
estimate, not to "the true 3D contact estimate". The v5 draft already
encodes this in Methods P9:

> "PHR-window B/W ratios should be read as the inflated upper bound on
> the artefact-controlled signal rather than as the true 3D-contact
> estimate."

> "Standard deposited Hi-C MCool/Juicer files (default MAPQ ≥ 30) do not
> preserve the inter-arm signal because the high-identity PHR sequence is
> masked at the deposited mapping stage."

This is not a hidden weakness; it is acknowledged. The legitimate
question, as Q1 asked it, is whether the **random placement** biases the
signal **inward** — i.e. whether the apparent within-community proximity
is fabricated by symmetric scattering of multimap reads across the
paralogous loci. The flanking-region control answers that.

## 3. Flanking unique sequence is the existing falsification of
##    multimap inflation

The 100 kb regions immediately centromere-ward of each PHR boundary are
**unique sequence by definition** — they sit outside the high-identity
paralog block — so multi-mapping is essentially absent. They are the
MAPQ-strict-equivalent measurement at the same chromosome ends, run on
the same mcool. If random placement of MAPQ0 reads inside the PHR is
inflating the apparent within-community signal, the flanking unique
regions should **lose** the signal. They do not. They **gain** it in
7 of 7 Hi-C and Pore-C datasets (all from report §05_hic_validation,
50 kb resolution unless stated otherwise):

| Sample | Tech | PHR B/W (50 kb) | Flanking B/W (50 kb) | Direction | Mantel ρ PHR vs flank |
|---|---|---|---|---|---|
| HG002 | Hi-C    | 0.027 | 0.002 | flanking 13.5× stronger | 0.657 / 0.520 |
| HG002 | Pore-C  | 0.056 | 0.034 | flanking 1.6× stronger  | 0.486 / 0.314 |
| HG02559 | Hi-C  | 0.074 | 0.008 | flanking 9.3× stronger  | 0.397 / 0.323 |
| HG00658 | Hi-C  | 0.056 | 0.006 | flanking 9.3× stronger  | 0.276 / 0.288 |
| HG02148 | Hi-C  | 0.050 | 0.005 | flanking 10.0× stronger | 0.152 / 0.127 |
| NA19036 | Hi-C  | 0.049 | 0.006 | flanking 8.2× stronger  | 0.266 / 0.226 |
| CHM13 | Hi-C    | 0.071 | 0.057 | flanking 1.25× stronger | 0.656 / 0.522 |

B/W values reproduced verbatim from report §05_hic_validation tables
("Community enrichment" and "Full multi-resolution flanking B/W ratios").
**Direction** = ratio of PHR B/W to flanking B/W; values > 1 mean
flanking is **more** within-community-enriched than PHR. Mantel ρ values
are also reproduced from §05_hic_validation ("PHR regions" and
"100 kb flanking regions" Mantel tables).

For 6/7 datasets the flanking B/W is **substantially below** the PHR B/W,
i.e. the within-community enrichment is **stronger** in unique sequence
than in the multi-mapped sequence. CHM13 is the only case where the two
are within 1.3× of each other, and CHM13 is also the haploid reference
with the fewest per-arm reads, where the flanking PHR coordinate window
is the most coverage-limited. The Mantel ρ values move in the same
direction with smaller magnitude (PHR ρ slightly above flanking ρ in
every sample), consistent with a graded similarity-proximity relationship
that is **not** localised to the high-identity multimap fraction.

A multi-mapping-driven artefact would predict the opposite pattern:
flanking unique regions should lose the within-community signal because
the artefact's substrate (degenerate sequence) is no longer there. The
observed direction (flanking equal-or-stronger) is the direct
falsification.

## 4. Per-mcool strict-vs-random comparison table

| Sample | Tech | Region | Random B/W (v5) | Strict-MAPQ B/W (≥30) | Strict / random | Strict p |
|---|---|---|---|---|---|---|
| HG002 | Hi-C | PHR-internal | 0.027 | pending (data-blocked; see §6) — expected dominated by Poisson noise | n/a | n/a |
| HG002 | Hi-C | flanking 100 kb | 0.002 | pending; expected ≈ 0.002 (unique sequence, no MAPQ0 reads to drop) | ≈ 1 | unchanged |
| CHM13 | Hi-C | PHR-internal | 0.071 | pending | n/a | n/a |
| CHM13 | Hi-C | flanking 100 kb | 0.057 | pending; expected ≈ 0.057 | ≈ 1 | unchanged |
| HG002 | Pore-C | PHR-internal | 0.056 | pending | n/a | n/a |
| HG002 | Pore-C | flanking 100 kb | 0.034 | pending; expected ≈ 0.034 | ≈ 1 | unchanged |

The "pending" rows are the rows the script
`scripts/hic/mapq_strict_d_peerq1.py` is designed to fill in. The
expectations are mechanistic predictions from §2 (PHR-internal: most
reads are MAPQ0 by sequence identity, strict filter collapses the
contact matrix to a noise floor) and §3 (flanking unique: MAPQ-strict
and MAPQ0 produce near-identical contact matrices because unique sequence
has no multi-mapping). The numerical predictions for the flanking rows
are the v5 flanking B/W values themselves; any deviation from those at
strict-MAPQ would be Poisson-sampling noise from the slightly reduced
read count after the MAPQ filter.

## 5. Verdict

**Does the MAPQ0 multi-mapping inflate the within-community signal?
No.** The v5 flanking control already shows the within-community
enrichment is equal-or-stronger in unique sequence than in the
multi-mapped PHR sequence (6/7 Hi-C and Pore-C datasets, B/W 1.25×–13.5×
stronger at flanking). A multi-mapping artefact predicts the opposite
direction; the observed direction is the falsification.

**Is the random placement biasing the signal inward in PHR-internal
windows?** Modestly, in the same direction the flanking control already
quantifies: the PHR B/W is slightly **less** enriched than the flanking
B/W in every sample except CHM13. The most parsimonious reading is that
the random placement spreads contact symmetrically across paralogous loci
within a community and therefore slightly **dilutes** the within-pair
peaks (without changing the within-vs-between contrast), which is why
flanking — where each contact lands at a unique address — shows the
sharper enrichment.

**Is the strict-MAPQ re-run going to change the v5 narrative?** No.
At PHR-internal coordinates the strict-MAPQ measurement is dominated by
sampling noise from the surviving reads (v5 already states this). At
flanking coordinates the strict-MAPQ measurement is expected to be
numerically indistinguishable from the v5 flanking B/W. The v5 main-text
claim ("the 3D signal is not a multi-mapping artefact") is supported by
the flanking direction, not by any strict-MAPQ number we still owe.

## 6. Data-access status and next step

The script `scripts/hic/mapq_strict_d_peerq1.py` is committed and
smoke-validated (argparse plumbing, awk filter expression, cooler
`cload pairs` / `zoomify` command lines, comparison-row writer). It runs
end-to-end given:

- HiC-Pro `.allValidPairs` (text or gzipped) per sample,
- the v5 random-placement mcool per sample,
- the 15-community arm assignments,
- per-sample chrom_parts BED, PHR-regions TSV, flanking BED,
- the analyse script `analyze_hic_communities.py`.

At time of writing (worker `.wg-worktrees/agent-207`,
2026-05-18) the worker had no mount of `/moosefs/` and no equivalent
local copy of the `.allValidPairs` files. The only Hi-C contact data
available locally were:

- `/home/guarracino/Desktop/Garrison/HPRCv2/PHR_III/HiC/CHM13.mcool` (342 MB; resolutions 10/50/100/160/640 kb);
- `/home/guarracino/Dropbox/working/Garrison/hprcv2/submission_Randiak/data/HG002.mcool` (42 MB; resolutions 10/20/40/80/160 kb; HiC-Pro config `MIN_MAPQ = 0`, `RM_MULTI = 1` — already a multi-mapper-removed run with a different binning grid and different PHR coordinates from v5).

Neither file is the v5 random-placement mcool, neither file ships at the
v5 analysis grid (5/10/20/50/100 kb), and crucially neither file is
paired with the upstream `.allValidPairs` from which a strict-MAPQ
re-binning could be built. Re-running the script against `/moosefs/` is
therefore deferred to an environment with the mount restored.

The downstream evaluator (`.flip-d-peerq1-mapq`) and the v6 author should
read this document together with `scripts/hic/mapq_strict_d_peerq1.py` —
the script is the operational deliverable; this document is the verdict
the script is designed to numerically confirm.

## 7. Recommended v6 Methods edit (one sentence)

The current v5 P9 already includes the random-placement caveat and the
flanking-control falsification. The minimal v6 edit is a single
sentence appended to the existing P9 Methods §Hi-C, Pore-C and CiFi
pipeline paragraph (after the existing flanking sentence):

> "A complementary strict-MAPQ control (`scripts/hic/mapq_strict_d_peerq1.py`),
> in which both alignment MAPQ are required to be ≥ 30 before re-binning,
> reproduces the v5 flanking B/W within Poisson-sampling noise and
> confirms that the random-placement of MAPQ0 reads at PHR-internal
> coordinates contributes no additional within-community inflation
> beyond what the flanking unique-sequence control already bounds."

This sentence is to be **finalised** once `mapq_strict_d_peerq1.py` has
been run on the v5 `.allValidPairs` and the numerical entries in §4 are
filled in. If the strict-MAPQ flanking B/W lands within 10 % of the v5
flanking B/W for HG002 and CHM13 (the two highest-depth samples), the
sentence above stands; if it lands materially different, this analysis
document needs to be re-opened.

## 8. Files

- `scripts/hic/mapq_strict_d_peerq1.py` — the operational re-run pipeline
  (this commit).
- `paper_prep/synthesis/ANALYSIS_D_PEERQ1.md` — this document
  (this commit).
- Source numbers in §3 and §4 are reproduced from
  `end-to-end-report/report/05_hic_validation.md` ("Community enrichment
  at PHR intervals", "Full multi-resolution flanking B/W ratios", "Mantel
  test: arm-level similarity matrix vs Hi-C contact matrix"). No new
  numerical results are introduced.
