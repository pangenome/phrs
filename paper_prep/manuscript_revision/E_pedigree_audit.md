# Cluster E Pedigree Audit: Circularity, Replication, and Claim Size

Date: 2026-06-17  
Task: `manuscript-revision-e1`  
Scope: judgment support only. This audit does **not** edit `submission/paper.tex`.

## Executive Recommendation

The Fig. 5 pedigree evidence is usable, but the claim should be framed as an
event-level consistency check, not as definitive proof of germline recombination
rates or independent recurrent exchange.

Recommended claim level:

> A T2T-resolved pedigree shows high-quality inter-chromosomal subtelomeric
> patches that are strongly enriched within the sequence-similarity communities
> inferred from the population graph.

Avoid:

> Inter-chromosomal recombination caught in the act.

The current evidence supports "pedigree-resolved patches consistent with recent
exchange", with strongest support in WashU T2T assemblies and weaker,
artifact-bounded support in fragmented CEPH1463 assemblies.

## Evidence Read

Primary manuscript text currently says the WashU pedigree has 538 high-quality
patches, 494 within the same sequence community, and a permutation null mean of
77.0% (`submission/paper.tex:239`, `submission/paper.tex:246`,
`submission/paper.tex:247`). The Fig. 5 caption repeats 494/538, Wilson CI
89.2--93.9%, and the Monte Carlo null (`submission/paper.tex:412`,
`submission/paper.tex:417`, `submission/paper.tex:419`). The Methods state
`odgi untangle nth-best=1`, the high-quality filter, and a "within-Leiden filter
as credibility constraint" (`submission/paper.tex:618`,
`submission/paper.tex:620`, `submission/paper.tex:622`), followed by the
within-Leiden permutation test (`submission/paper.tex:628`).

The source report defines a patch as a contiguous run of untangle segments with
the same parent chromosome and haplotype, applies the high-quality filter
`is_interchr=True`, `min_score >= 0.8`, and `500 bp <= size <= 100 kb`, and says
Leiden membership is used for validation/reporting (`end-to-end-report/report/14_pedigree_recombination.md:19`,
`end-to-end-report/report/14_pedigree_recombination.md:20`,
`end-to-end-report/report/14_pedigree_recombination.md:30`). It also states
the full table contains all 5,984 high-quality patches, including cross-community
and unknown calls (`end-to-end-report/report/14_pedigree_recombination.md:32`).

The survey makes the same status explicit: WashU has 538 high-quality
inter-chromosomal patches, 494 within Leiden, while fragmented CEPH1463 has
2,775 hifiasm patches with 324 within Leiden and 2,671 verkko patches with 359
within Leiden (`paper_prep/surveys/SURVEY_14_pedigree_recombination.md:20`,
`paper_prep/surveys/SURVEY_14_pedigree_recombination.md:22`,
`paper_prep/surveys/SURVEY_14_pedigree_recombination.md:23`,
`paper_prep/surveys/SURVEY_14_pedigree_recombination.md:24`). It identifies
`all_pedigrees_patches.tsv` and the `community_status` column as the downstream
filtering field (`paper_prep/surveys/SURVEY_14_pedigree_recombination.md:127`,
`paper_prep/surveys/SURVEY_14_pedigree_recombination.md:128`), and says a patch
is `within_community` iff the query and reference arms have the same Leiden
label (`paper_prep/surveys/SURVEY_14_pedigree_recombination.md:177`,
`paper_prep/surveys/SURVEY_14_pedigree_recombination.md:178`).

The upstream script confirms this order of operations. Patch construction and
pattern classification are based on untangle state runs and immediate neighbors
(`/moosefs/guarracino/HPRCv2/scripts/pedigree/analyze-pedigree-recombination.py:147`,
`/moosefs/guarracino/HPRCv2/scripts/pedigree/analyze-pedigree-recombination.py:223`,
`/moosefs/guarracino/HPRCv2/scripts/pedigree/analyze-pedigree-recombination.py:236`).
Community status is assigned afterward by comparing query and reference arm
Leiden labels (`/moosefs/guarracino/HPRCv2/scripts/pedigree/analyze-pedigree-recombination.py:441`,
`/moosefs/guarracino/HPRCv2/scripts/pedigree/analyze-pedigree-recombination.py:452`).
The high-quality inter-chromosomal set is defined from `is_interchr`,
`min_score`, and patch size, not from Leiden status
(`/moosefs/guarracino/HPRCv2/scripts/pedigree/analyze-pedigree-recombination.py:579`,
`/moosefs/guarracino/HPRCv2/scripts/pedigree/analyze-pedigree-recombination.py:580`).

I also counted the upstream table directly:

| Dataset | HQ inter-chr patches | Within | Cross | Unknown | Within fraction |
|---|---:|---:|---:|---:|---:|
| WashU | 538 | 494 | 44 | 0 | 91.8% |
| CEPH1463 hifiasm | 2,775 | 324 | 1,735 | 716 | 11.7% |
| CEPH1463 verkko | 2,671 | 359 | 1,525 | 787 | 13.4% |
| All | 5,984 | 1,177 | 3,304 | 1,503 | 19.7% |

Within WashU, the 538 patches occur in only three assayed child-haplotype
transmissions: PAN027 maternal from PAN010 (167 total, 165 within), PAN027
paternal from PAN011 (61 total, 57 within), and PAN028 maternal from PAN027
(310 total, 272 within). Pattern counts over all 538 are: 262 `acros_like`,
136 `gene_conversion_like`, 120 `sandwich_same_hap`, 18 `crossover_like`, and
2 `complex`. The within-Leiden subset is 229, 133, 115, 16, and 1 respectively,
matching the survey's filtered counts (`paper_prep/surveys/SURVEY_14_pedigree_recombination.md:28`,
`paper_prep/surveys/SURVEY_14_pedigree_recombination.md:38`).

## E1. Circularity: Is Leiden Used to Define the Tested Set?

Answer: for WashU, not for the tested 538-patch denominator. The high-quality
WashU set is defined by untangle inter-chromosomal state, patch score, and patch
size. Leiden membership is assigned afterward and the statistic is 494/538
within-community. That makes the 494/538 test a valid concordance/enrichment
test against the independently derived HPRCv2 community labels, as long as the
Methods clearly state that cross-community WashU patches remain in the
denominator.

However, the current Methods wording is ambiguous because "within-Leiden filter
as credibility constraint" appears before the permutation-test sentence. A
reader could infer that the tested set was pre-filtered to within-Leiden calls,
which would make 494/538 circular or undefined. The source report and survey
also say "Only within-community patches are reported", which is true for the
curated report tables but not for the 494/538 denominator.

Recommended correction is documentary, not analytical: separate "scoring the
tested set" from "curating examples".

### Methods Sentence Options

Preferred:

> For WashU, all high-quality inter-chromosomal patches were retained for the
> primary concordance test; Leiden community labels, inferred independently from
> the HPRCv2 population graph, were assigned only after patch calling, and the
> reported statistic is the fraction of these patches whose query and donor arms
> share a Leiden community.

Stronger anti-circularity wording:

> The within-Leiden label was not used to call WashU patches or to define the
> 538-patch denominator; it was used only as an orthogonal annotation for the
> enrichment test and for selecting high-confidence examples shown in the text.

If keeping "credibility constraint":

> For display and mechanistic examples, we emphasized patches whose query and
> donor arms fell in the same independently inferred Leiden community, but the
> enrichment statistic was computed over all 538 WashU high-quality
> inter-chromosomal patches.

Caption option:

> Of 538 high-quality inter-chromosomal patches called without using community
> labels, 494 were subsequently annotated as within a sequence-similarity
> community.

Avoid:

> within-Leiden filter as credibility constraint; the within-Leiden fraction was
> tested...

That sequence invites the circularity objection even if the implementation did
not do that for WashU.

## E2. Unit of Replication

The true replication unit for the 494/538 statistic is the **patch**, but those
patches are not independent biological replicates. They are nested within three
child-haplotype transmissions in one T2T pedigree, and concentrated in a smaller
number of chromosome arms and communities. The most defensible wording is:
"538 high-quality patch calls from three assayed transmissions in the WashU
T2T pedigree."

The true biological replication unit is not 538 meioses or 538 people. For
WashU, it is at most three child-haplotype transmissions:

| Transmission label | HQ patches | Within-Leiden |
|---|---:|---:|
| PAN027 maternal hap1 vs PAN010 mother | 167 | 165 |
| PAN027 paternal hap2 vs PAN011 father | 61 | 57 |
| PAN028 maternal hap1 vs PAN027 mother | 310 | 272 |

For mechanism subclasses, the replication is even smaller. The 16
within-community `crossover_like` patches are patch-level observations within
these transmissions, not 16 independent meioses. The survey also cautions that
T2T pedigrees are rare, WashU has `n = 4` samples, `sandwich_same_hap`
interpretation is open, and no de-novo-vs-inherited split is available
(`paper_prep/surveys/SURVEY_14_pedigree_recombination.md:273`).

### Replication Sentence Options

Preferred:

> Because multiple patches can occur within the same inherited haplotype and
> the WashU analysis samples three child-haplotype transmissions from one
> pedigree, we treat patches as the counting unit for the concordance statistic,
> not as independent meiotic replicates.

Short:

> The denominator is patch calls, nested within three WashU child-haplotype
> transmissions.

For Results:

> Across three assayed WashU transmissions, 494 of 538 high-quality
> inter-chromosomal patch calls mapped to the same independently inferred
> sequence community as their donor arm.

Avoid:

> 538 recombination events...

Use "patch calls" unless the sentence is limited to a pattern subclass and
explicitly hedged as "crossover-like" or "gene-conversion-like".

## E3. Is Lower-Quality Assembly Artifact Bounded?

Answer: partially bounded, but not fully quantified for the WashU calls.

The fragmented CEPH1463 assemblies provide a useful empirical bound on how much
untangle can degrade in non-T2T assemblies: only 12--13% of their high-quality
inter-chromosomal patch calls are within Leiden, compared with 92% in WashU.
This strongly indicates that fragmented assemblies produce many artifactual or
uninterpretable inter-chromosomal calls. It also supports the current decision
to use CEPH1463 only through cross-assembler parent-by-chromosome-pair features.

But the artifact fraction within WashU is not directly bounded by an orthogonal
truth set. The 44/538 WashU cross-community calls provide an upper-bound-like
count of discordant patch calls under the community model, but they are not
necessarily all artifacts. Conversely, some within-community WashU calls could
still be assembly, graph, or untangle artifacts because the Leiden filter checks
biological plausibility, not molecular truth. The source text explicitly says
lower-quality assemblies have many more likely artifacts, but does not validate
each WashU event with independent reads or a second T2T assembler.

A defensible statement is:

> Assembly-fragmentation artifacts are empirically bounded for CEPH1463 by the
> low within-community fraction and by the hifiasm/verkko intersection, but the
> per-event artifact rate in WashU is not directly estimated; WashU's T2T status
> and high within-community enrichment make the signal credible at aggregate
> level.

Avoid claiming:

> The artifact fraction is known to be 8%.

The 44/538 discordant WashU calls are not a calibrated false-positive rate.

### Artifact Sentence Options

Preferred:

> Fragmented CEPH1463 assemblies show that untangle can produce many
> inter-chromosomal patch calls outside the population-defined communities
> (12--13% within-Leiden), so we use those assemblies only as a
> cross-assembler validation set; the WashU T2T analysis is interpreted at the
> aggregate enrichment level rather than as a fully validated catalog of
> individual recombination events.

Short Results caveat:

> The signal is aggregate: individual patch calls remain candidates until
> confirmed by orthogonal long-read or read-backed recombination evidence.

Methods caveat:

> Cross-community and unknown-community calls were retained in the WashU
> denominator but not used as mechanistic examples.

Avoid:

> the same signal increases in lower-quality assemblies

This is imprecise. The signal does not "increase"; the number of
inter-chromosomal patch calls increases while within-community concordance
collapses to 12--13%.

## Claim-Title Recommendations

Best section title:

> Pedigree-resolved patches support recent subtelomeric exchange

More conservative:

> T2T pedigree untangling links recent subtelomeric patches to sequence
> communities

If the paper wants a stronger narrative title:

> A complete pedigree captures candidate exchange tracts within subtelomeric
> communities

Avoid:

> Recombination in a complete pedigree

This is acceptable as a broad section label but underspecifies that the
evidence is patch-based and candidate-level.

Avoid:

> Inter-chromosomal recombination caught in the act

This overstates the current evidence. It implies direct observation of meiotic
exchange events, while the analysis observes inherited sequence-state patches
from assemblies and classifies them by untangle context.

Recommended Fig. 5 title:

> Pedigree-resolved candidate subtelomeric exchange

Recommended caption lead:

> Candidate inter-chromosomal subtelomeric patches in the WashU
> three-generation T2T pedigree.

## Manuscript-Ready Replacement Blocks

### Results Paragraph Option

> To look for recent exchange directly, we applied `odgi untangle` to the
> WashU three-generation T2T pedigree, comparing each inherited child haplotype
> with the two haplotypes of the transmitting parent. Patch calls were made
> without using the HPRCv2 community labels. Across three assayed
> child-haplotype transmissions, 494 of 538 high-quality inter-chromosomal patch
> calls (92%) subsequently mapped between arms assigned to the same
> sequence-similarity community, above a marginal-aware permutation null. These
> patch calls include candidate gene-conversion-like tracts and
> crossover-like state switches, but they remain assembly-derived candidates
> rather than a per-event recombination truth set.

### Methods Block Option

> Pedigree patches were called from `odgi untangle nth-best=1` output by
> merging contiguous child-flank intervals assigned to the same transmitting
> parent chromosome and haplotype. The high-quality WashU denominator comprised
> all inter-chromosomal patches with `min_score >= 0.8` and
> `500 bp <= size <= 100 kb`. Leiden labels from the independently inferred
> HPRCv2 arm-level community graph were assigned only after patch calling. The
> primary statistic was the fraction of high-quality WashU patches whose query
> and donor arms shared a Leiden community, tested against a permutation null
> preserving per-arm patch-count marginals. Within-community status was used to
> select mechanistic examples and, for CEPH1463, as part of a conservative
> cross-assembler validation criterion.

### Caption Block Option

> Candidate inter-chromosomal subtelomeric patches in the WashU
> three-generation T2T pedigree. Purple marks high-quality inter-chromosomal
> patch calls from `odgi untangle`; community labels were assigned after patch
> calling. Of 538 patch calls across three child-haplotype transmissions, 494
> (92%; Wilson 95% CI 89.2--93.9%) joined arms in the same independently inferred
> sequence-similarity community, above a marginal-aware Monte Carlo null
> (mean 77.0%, 95% CI 75.4--78.8%; `p < 1e-4`).

## Decision Record

- Circularity: **No for WashU 494/538 if the denominator remains all
  high-quality inter-chromosomal patches; yes-risk in wording if "within-Leiden
  filter" is described before the denominator/test.**
- Replication unit: **patch calls nested within three child-haplotype
  transmissions in one T2T pedigree; not independent meioses.**
- Artifact bound: **bounded qualitatively and for CEPH1463 by cross-assembler
  concordance plus 12--13% within-Leiden rates; not quantitatively bounded as a
  WashU per-event false-positive rate.**
- Claim size: **use "candidate patches", "consistent with recent exchange",
  "aggregate enrichment", and "pedigree-resolved"; avoid "caught in the act" and
  avoid treating 538 as independent recombination events.**
