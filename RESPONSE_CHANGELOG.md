# Reviewer-response change log

Manuscript text is **corrected**, never phrased as a reply to reviewers. This log
tracks what changed and why, so the response letter can be written from it.

Naming: `PEDIGREE_NAMING.md` (MOM=PAN010, DAD=PAN011, SELF=PAN027, KID=PAN028,
PARTNER=KID's father, unsampled).

## Decisions

- One ectopic event is retained: SELF.p.chr9q terminus derived from DAD.h2.chr3q.
- The chr1p/chr5q signal is ancestral homology, not an event.
- KID-side candidates are withdrawn: PARTNER is unsampled, and those scans had used
  the paternal chromosome copy anyway.
- Fig5 is simplified to the ectopic event: A pedigree, B identity-crossover zoom,
  C paternal scan, and a new **D: pairwise alignments of the sequences involved**
  (SELF.p.chr9q vs DAD.h2.chr3q / DAD.h2.chr9q / DAD.h1.chr9q). Panels D/E
  (maternal scan, grandchild scan) retired. 1p/5q relabelled ancestral; 21p/4p
  dropped.
- Population-frequency support for the standing arrangement: NOT pursued.

## Done

| # | change | where |
|---|---|---|
| 1 | Results: one documented exception, not "two putative autosomal exchanges"; chr1p/chr5q stated as ancestral homology | `paper.tex` §Ectopic subtelomeric recombination within a pedigree |
| 2 | Results: maternal-scan paragraph trimmed; granddaughter candidates removed | same section |
| 3 | Results: closing paragraph now says the family transmits an ectopic product rather than a de novo exchange with resolved breakpoints | same section |
| 4 | Fig5: panels D/E removed; legend 21p/4p dropped, 1p/5q relabelled ancestral; caption rewritten, title now "Pedigree view of an ectopic subtelomeric exchange" | `paper.tex` Fig5 block |
| 5 | Methods: "Inference rule" paragraph (best parental match = source; exceptions PAR1, acrocentric arms, this case) | `paper.tex` §Pedigree direct-alignment windows |
| 6 | Methods: "Haplotype labelling" paragraph (per-chromosome parent of origin; maternal copy resolved by identity; unsampled PARTNER ⇒ KID paternal scans not reported) | same |
| 7 | Methods: Fig5C--E references reduced to C | same |
| 8 | Evidence record with all numbers, retractions and revisions rev 1-12 | `REVIEWS.txt` |
| 9 | Canonical naming + per-chromosome maternal table | `PEDIGREE_NAMING.md` |
| 11 | **Fig5B regenerated** with the naming convention: axis `SELF.p.chr9q position (Mb)`, series `DAD.h2.chr3q (source)` and `DAD.h2.chr9q (drifted copy)`, tract shaded 20 kb; the retired KID 21p/4p panel removed from the script | `submission/fig/MainFigures/Fig5B_zoom9q3q.{pdf,png}`, `submission/scripts/figures/make_fig5_zoom.R` |
| 15 | Fig5A regenerated with DAD/MOM/SELF/KID labels; Fig5C paternal ribbon regenerated with `SELF.p` / `DAD.h1` / `DAD.h2`; retired maternal+KID panel assets deleted and their 3 supporting-material rows dropped, survivors relabelled | `submission/fig/MainFigures/`, `submission/scripts/figures/`, `submission/supporting_material/Fig5_pedigree_supporting_material.tsv` |
| 16 | Response letter: flagger-QC point and chr1p/chr5q standing-homology point added; aligner-disagreement argument removed | `RESPONSE_LETTER_DRAFT.md` |
| 14 | Margin/support table written and a Methods paragraph added: PAR1 = +0.8 pp margin (obligate control); chr9q/chr3q = +0.3 pp margin but heterozygous-donor structure; acrocentric = pattern-based (multi-source NAHR), margin not applicable | `submission/supporting_material/Fig5_margin_table.tsv`, `paper.tex` |
| 13 | Methods de-meta-ified: the two emphatic paragraphs ("Inference rule", "Haplotype labelling") replaced by one neutral factual paragraph ("Chromosome-copy assignment"); all reviewer-facing reasoning moved out of the manuscript into `RESPONSE_LETTER_DRAFT.md` | `submission/paper.tex` |
| 12 | Fig5C ribbon script relabelled (`SELF.p (PAN027 paternal)`, `DAD.h1/h2`) and the maternal/KID gen calls removed; **regeneration blocked here: inkscape is not installed** | `submission/scripts/figures/make_fig5_panels.sh` |
| 10 | New **Fig5D** pairwise-alignment panel built and added to the figure block + caption; asset, box coordinates and generator script committed | `submission/fig/MainFigures/Fig5D_pairwise_9q3q.{pdf,png}`, `Fig5D_highlights.tsv`, `submission/scripts/figures/make_fig5d_pairwise.sh` |

## Audit result: PAN028 label exposure (pending item 3, partially done)

`all_pedigrees_patches.tsv`, WashU dataset = 538 patches:

| dataset | patches | comment |
|---|---|---|
| `PAN028 maternal (hap1) vs PAN027` | **310 (58%)** | KID-side; hap1 is maternal only on chr2,6,8,10,11,13,16,17,19,20 |
| `PAN027 maternal (hap1) vs PAN010` | 167 | SELF side; labels consistent (`PAN027_mat`) |
| `PAN027 paternal (hap2) vs PAN011` | 61 | SELF side; labels consistent (`PAN027_pat`) |

Of the 310 KID patches, **228 (74%) sit on chromosomes where hap1 is the PATERNAL
copy**, i.e. they are calls about PARTNER's unsequenced chromosome and cannot support
a transmission claim; only 82 are label-correct. Wrong-label counts by chromosome:
chr15 94, chr21 38, chr14 25, chr22 18, chr5 15, chr9 10, chr1 10, chr3 7, chr7 7, chr4 4.

`crossover_like` patches: 18 total, 15 in KID. **7 of those 15 are wrong-label**
(chr1, chr3, chr4, chr5, chr22) — including the two the drafts quote explicitly:
chr3 q 366,221-368,310 → chr9q:2 (2,089 bp: the call already retracted in rev 10) and
chr3 q 262,953-290,922 → chr7p:2 (27,969 bp, the "largest" event). The 8 label-correct
crossover-like calls are all chr13p → chr15p/chr21p/chr22p, i.e. acrocentric p-arm sharing.

Consequences:
- The submitted manuscript (`submission/paper.tex`) contains **no** 494/538 or
  "13 of 16 crossover-like" claim — checked: clean. So the submitted pedigree text
  needed only the corrections already made (items 1-7).
- The claim still stands in working documents and must be fixed before reuse:
  `paper_prep/figures/fig4/caption.md` line 9 ("494/538 = 92 %"),
  `paper_prep/synthesis/NATURE_DRAFT_v4.md` line 48 and
  `NATURE_DRAFT_v5.md` line 53 ("Thirteen of sixteen crossover-like events are in
  PAN028, confirming that meiotic-resolution inter-chromosomal breakpoints transmit
  across generations"; "largest 27.97 kb on PAN028 maternal chr3q").
- The 92% within-community figure is a property of query/ref arm membership, so it is
  not automatically invalidated, but the 228 wrong-label calls must be re-run on the
  correct maternal copy before the count is reported as evidence of transmission.

### Headline KID events rechecked on the correct maternal copy

| claim (child-side) | on KID.p (what was queried) | on KID.m (correct copy) | verdict |
|---|---|---|---|
| 27,969 bp chr3q -> chr7p, "largest crossover-like event" | SELF chr7p 99.1-99.4% over 100% coverage; SELF chr3q 11-13% covered | **SELF.p.chr3q 99.94% over 100% coverage**; chr7p 96.6% over 88.9% | **retired** — PARTNER's chr3 carries a chr7p-derived segment; on the daughter's maternal copy the region is plain homolog transmission, chr7p losing by 3.3% |
| 2,089 bp chr3q -> chr9q | (already retracted in rev 10) | SELF.p.chr3q 99.70% vs chr9q 99.44% | **retired** |
| 8 label-correct crossover-like calls (chr13p -> chr15p/chr21p/chr22p, <=2.7 kb) | n/a (hap1 IS maternal on chr13) | same | **survive** as acrocentric p-arm sharing |

The 27.97 kb case is the clean illustration of the whole defect: the same coordinates give
opposite answers on the two copies, and only the maternal copy is interpretable.

### Recheck of the 228 wrong-label KID patches against the existing untangle output

The PAN028 hap2 flanks are **already present** in
`PAN028_vs_PAN027.e50000.m1000.bed.gz` (435,332 rows) — the patch table simply only
extracted the hap1 rows. No new untangle run was needed. Mapping each hap1 patch
interval to the hap2 flank by telomere-anchored coordinates (the two haplotypes differ
in arm length, so start-anchored offsets are invalid) and taking the rank-1 untangle row
at those coordinates:

| outcome on the **correct maternal copy** | count |
|---|---|
| homolog wins (rank-1 is same chromosome) — call does not survive | **168** |
| interval outside the maternal flank (haplotype length difference) — untestable | 54 |
| **survives** as rank-1 inter-chromosomal | **6** |

The six survivors all have **acrocentric p-arm donors** (chr13p x3, chr22p x2, chr21p x1),
i.e. the rDNA exchange system, not new subtelomeric events.

### Corrected community statistic — it gets stronger

| set | patches | within Leiden community |
|---|---|---|
| as published | 538 | 494 (92%) |
| SELF-side (labels consistent) | 228 | 222 |
| KID, label-correct (hap1 = maternal) | 82 | 81 |
| KID, rechecked survivors on the maternal copy | 6 | 6 (acrocentric) |
| **corrected total** | **~316** | **~309 (98%)** |

The wrong-copy calls were disproportionately *cross*-community (they describe the
unsampled father's chromosomes, which need not respect the mother's partition), so
removing them raises the within-community rate from 92% to **98%** on a cleaner set.

Caveats to carry into the response: 54 intervals are untestable without the paternal
genome or longer flank windows; the recheck uses rank-1 logic plus telomere-anchored
mapping rather than a re-run of the upstream merge pipeline, so the definitive numbers
should come from re-running that pipeline with the per-chromosome copy selection.

### The dichotomy that survives (sex-specific germline architecture)

SELF-side patches, split by transmitting parent (labels verified, so these are supported):

| transmission | patches | acro<->acro | non-acrocentric |
|---|---|---|---|
| SELF.m <- MOM | 167 | **162 (97%)** | 5 (2.0-4.2 kb gene-conversion/sandwich-like) |
| SELF.p <- DAD | 61 | 42 | **19** |

The 19 paternal non-acrocentric calls contain exactly the two well-characterised cases:
chrXp->chrYp (PAR1, obligate) and chr3q->chr9q 3,414 bp crossover_like (the event), plus
small gene-conversion-like calls and chr9q->chr15q/16q. The 5 maternal non-acrocentric
calls are 2-4 kb gene-conversion/sandwich-like, including the same chr3q->chr9q pair we
showed to be base-level ties.

So the pedigree reads as: **paternal line = obligate PAR1 exchange plus one subtelomeric
ectopic exchange; maternal line = acrocentric-acrocentric exchange** - and nothing else.
This matches prior expectations in two independent ways: male meiotic recombination is
telomere-biased, and Robertsonian translocations are overwhelmingly of maternal origin,
i.e. ectopic recombination between acrocentric short arms arises preferentially in the
female germline. The acrocentric signature is also the majority of the maternal calls
(162), not a single lucky event, and recurrent rDNA NAHR is by nature repeated rather
than one-off, so a single pedigree is a weak constraint here.

Caveat to keep: acrocentric short arms are the most error-prone regions, so the
cross-assembler (hifiasm + verkko) filter in CEPH1463 remains the right guard, alongside
the pangenome-scale acrocentric communities (C15, C7) as population support.

### Discussion bullet candidate: sex-specific germline architecture (with the literature numbers)

Observation in our data (SELF-side, labels verified): maternal transmission 162/167
acrocentric-acrocentric; paternal transmission carries PAR1 plus the one chr3q/chr9q
subtelomeric exchange, with 19 non-acrocentric patches vs 5 on the maternal side.

Literature, with the case counts so the claim is proportionate:
- De novo Robertsonian translocations form mostly during **oogenesis**; the parental-origin
  and timing analysis of de novo ROBs (Bandyopadhyay et al., AJHG 2002;71:1456-1462)
  reports a maternal-origin predominance, with the common rob(13;14) and rob(14;21)
  among the characterised cases - a cohort of ~30 de novo ROBs, i.e. small.
- Reciprocally, **de novo reciprocal translocations** arise predominantly in
  spermatogenesis (Hoeckner et al., Cytogenet Genome Res 2012;136:242-245).
- Large compendia give the population-scale context: the collaborative European
  series covers **52,965 amniocenteses** (Ferguson-Smith & Yates, Prenat Diagn
  1984;4:5-45) and establishes the maternal-age dependence of acrocentric
  rearrangements; Kovaleva's compendium of ROB carriers with known parental origin
  reports 37 maternal vs 30 paternal (67 carriers, female-biased carrier sex ratio
  0.64), so across all carriers the origin ratio is close to parity even though the
  *de novo* class is oogenic.
- Acrocentric short arms are the most error-prone regions of the assembly, so the
  claim should be phrased as an observation in these data plus cross-assembler support
  (hifiasm + verkko in CEPH1463), not as a definitive rate.

Suggested Discussion sentence: "In this pedigree the asymmetry is sex-specific: the
paternal transmission carries the obligate Xp/Yp exchange and a single subtelomeric
ectopic exchange, whereas the maternal transmission carries acrocentric-acrocentric
exchange - consistent with de novo Robertsonian translocations arising predominantly
during oogenesis and with the telomere bias of male meiotic recombination."

### Acrocentric caveat — the numbers for Methods

Joining all 538 WashU patches (chromosomal coordinates via the flank map) to the
**HiFi** flagger calls on the matching contig, at the patch interval:

| query class | patches | flagger overlap | rate |
|---|---|---|---|
| acrocentric | 459 | 1 | **0.2%** |
| non-acrocentric | 79 | 0 | 0.0% |
| all | 538 | 1 | 0.2% (1.55 Mb acrocentric + 0.46 Mb non-acrocentric clean) |

So at HiFi resolution the patch intervals are essentially untouched by flagger, including
the acrocentric ones (1/459). The earlier "66/538 touch problematic regions" figure came
from the ONT-based calls, which mark the terminal tens of kb of every chromosome end in
every sample and therefore say nothing patch-specific. Methods sentence candidate:
"Patches were checked against the HiFi-based flagger calls: 1 of 538 overlaps such a call
(the acrocentric subclass 1 of 459), and the remaining overlap in the ONT-based calls is a
generic end-of-chromosome effect present in every sample."

### D-PeerQ1 / D-PeerQ2 response drafted (subagent) + a number fix it surfaced

`paper_prep/synthesis/RESPONSE_Q1_Q2_DRAFT.md` (subagent, reviewed):
- Q1: the pipeline places each read pair at ONE position, so multi-mapping cannot spread
  mass; placement error is symmetric with respect to community membership, so it adds noise
  but no directional enrichment; falsification is the unique-sequence flanking control; the
  a-fortiori argument (few engaged copies can dilute a signal but cannot manufacture one) is
  included. Closure experiment correctly specified against the upstream `.allValidPairs`
  stage, not the binned mcool (mcool carries no per-read MAPQ).
- Q2: not softened. Framed as a Discussion hypothesis with direct precedent - Liang et al.,
  Nature 619 (2023), doi:10.1038/s41586-023-06323-x, "Complementary Alu sequences mediate
  enhancer-promoter selectivity" (repeat-derived RNA determining which enhancer contacts
  which promoter), plus the aggregate-vs-sparse-participation argument.

NUMBER FIX surfaced by the subagent's audit: the widely quoted "PHR B/W 0.027 -> flanking
0.0031, 9x" pairs the 50 kb PHR value with the 100 kb flanking value. Report table
`end-to-end-report/report/05_hic_validation.md` gives HG002 flanking B/W 0.001 (5/10/20 kb),
0.002 (50 kb), 0.003 (100 kb), and PHR B/W 0.029 (5 kb), 0.029, 0.028, 0.028, 0.027 (50 kb),
0.021 (100 kb). Resolution-consistent pairs:
  10 kb: 0.029 -> 0.001 = 29x    50 kb: 0.027 -> 0.002 = 13.5x    100 kb: 0.021 -> 0.003 = 7x
Corrected in: `RESPONSE_Q1_Q2_DRAFT.md`, `NATURE_DRAFT_v4.md`, `NATURE_DRAFT_v5.md`,
`figures/fig3/caption.md`, `NATURE_DRAFT_OUTLINE.md`. Left as-is (historical): `NATURE_DRAFT_v2.md`
(frozen revision) and `REVISION_LOG_v5.md` (log quoting the superseded wording).

### D-PeerQ3 resolved (arm-type asymmetry)

New script `scripts/cladistics/peerq3_arm_type_asymmetry.py` + table, draft in
`RESPONSE_Q3_DRAFT.md`. Within-community pairs across the 41 arms / 15 communities:
P-P 19 (exp 14.8, p=0.148), Q-Q 22 (exp 13.5, **p=0.019**), P-Q 17 (exp 29.7,
**p=0.005 depleted**; 20,000 arm-type permutations with communities fixed).
Two communities are arm-type pure: **C6 all Q-Q (15 pairs, zero P-Q, p=0.021)** and
**C7 all P-P (10 pairs, zero P-Q, p=0.048)**; mixing is confined to C3 and C11; C1
(4q/10q) is Q-Q and C2 (10p/18p) P-P as the concern anticipated.
Orientation half: **resolved at the correct scale.** Using the HPRCv2 all-vs-all flank
alignments (`all-vs-all.1Mb.p95`, 46 CHM13 arm-end PAFs), blocks >=5 kb, homolog-vs-homolog
excluded: **98.6% of 9.6 Gbp of inter-arm aligned sequence is in the same telomere-relative
orientation** (764 arm pairs; P-P 97.8%, Q-Q 100.0%, P-Q 98.4%); per-pair MEDIAN is 100% and
only 7/764 pairs are below 90%. Modules are reassorted in order but inserted in a consistent
telomere-to-centromere orientation. Exceptions are known inverted-repeat regions: chr8p vs
chr1p 0.0%, chr8p vs chr2q 5.7%. Direct-repeat bias is what gene conversion/NAHR predicts
(inverted-repeat exchange gives dicentrics and is selected against).
Script `scripts/cladistics/arm_module_orientation.py`, table `arm_module_orientation.tsv`.

Superseded: the earlier 100 kb-window pass (WITHDRAWN). The first pass (CHM13 100 kb arm ends,
sensitive all-vs-all) gave 8/17 P-Q pairs alignable, 6 same / 2 opposite, but every such
alignment is a 300-2,000 bp scattered repeat match at 78-95% identity and **none is
telomere-offset matched** (0 of 42 alignments). Strand on a scattered repeat reflects which
paralog copy was hit, not arm-end orientation, so the call is retracted. P-Q orientation is
undefined at 100 kb arm-end scale; it needs 1 Mb arm ends with offset-aware block selection.
The arm-type composition result (P-Q depleted p=0.005; C6 all-Q, C7 all-P) is unaffected.

### Module-order analysis: the architecture is 'same parts, consistent orientation, arbitrary order'

`scripts/cladistics/arm_module_order.py` + `arm_module_order.tsv` (HPRCv2 all-vs-all flank
alignments; per-arm partner midpoint profiles; Spearman across shared partners):
- repertoire Jaccard same-community 0.828 vs different 0.571 (p = 4.7e-10) - communities share
  module content;
- order concordance same-community +0.209 vs different +0.018 (p = 0.02) but **NOT robust**:
  stratified by number of shared partners the advantage vanishes and reverses (n=24: -0.27 vs
  -0.28; n=25: -0.94 vs +0.37). Reported as "order not conserved".
Combined with orientation (98.6% same): same parts, consistent orientation, arbitrary order.

### Full architectural pattern completed

`scripts/cladistics/arm_module_architecture.py` (+ tiling/blocks/inversions TSVs) and
`submission/scripts/figures/make_module_map_architecture.R` (module-map figure, 26 arms tiled):
- BLOCK SIZE: same-community partners have longer shared blocks (median 13.3 kb vs 9.9 kb) and
  more blocks per pair (1,919 vs 318 median; count confounded by sampled individuals).
- INVERSIONS: only 3 net-inverted pairs (chr8p/chr1p, chr8p/chr2q, chr2q/chr8p), each built
  from 2,645-3,832 small blocks with largest block 36-50 kb -> repeat-family inversion, not a
  single large inverted segment.
- MODULE MAPS: 25 kb bin tiling; the most frequent top partners are acrocentric (chr14p,
  chr22p), consistent with rDNA short arms as the most promiscuous module donors.
Architecture summary: same module repertoire, consistent telomere-relative orientation
(98.6%), arbitrary order; community membership predicted by repertoire (Jaccard 0.83 vs 0.57,
p=5e-10) and block length, not arrangement.

### Statistical controls and bibliography closed (D-M4, D-M5, D-M9, D-M12, D-Bib, D-PeerQ2b)

See `RESPONSE_STATS_AND_BIB.md` for the letter text. Summary of what was actually run:
- **D-M4 CLOSED**: 10,000-permutation null preserving per-arm marginals. Published 469/538 =
  0.872 vs null 0.677 [0.664-0.693]; corrected 293/310 = 0.945 vs null 0.840 [0.832-0.852];
  both at the 100th percentile, p < 1e-4. `ANALYSIS_D_M4.json`.
- **D-M5 CLOSED**: arm-level Mantel replaces the non-independent per-pair Spearman.
  `scripts/hic/mantel_mouse_d_m5.py`: leptotene +0.41/+0.29, zygotene +0.58/+0.34,
  pachytene +0.63/+0.31, diplotene +0.33/+0.30 (20 kb/50 kb), all p < 1e-4, null 95% +-0.11.
- **D-M9 DONE and consequential**: character bootstrap B=1000 (`char_bootstrap_d_m9.R`).
  PAR1 100%, PAR2 100%, P10_18 100%, ACRO_p 54.4% NJ / 87.8% UPGMA; **TIGHT_q and DUX4 are
  not monophyletic in the full tree** (0% by construction) - the old 100% "sensitivity"
  number measured partition recovery, not clade monophyly. Manuscript wording for those two
  groups must be community assignment, not tree support.
- **D-M12 CLOSED**: CIs in `scripts/ci/results_d_m12.json` for pedigree fraction (published and
  corrected), F_ST jackknife, Mantel bootstrap, and the new mouse arm-level Mantel.
- **D-Bib CLOSED**: Vollger2023 is a *different* paper (Science 376:eabj6965) from
  concerted_evolution_nahr_Vollger2023 (Nature 617:325); Smith1976 already present; Cech2004
  added to the synthesis bib (already in the submission bib as cech2004chromosomeend).
- **D-PeerQ2b**: roadmap statement with the required analysis and the numbers that would
  settle it; no analysis this round.

### Authoritative recount with the repo's own merge logic

CORRECTION TO MY EARLIER CLAIM: nothing about the *environment* blocks the recount. Tool
versions are pinned in Methods (wfmash v0.23.0; pggb v0.7.4 commit 498c5d7 with seqwish
v0.7.11, smoothxg v0.8.2, odgi v0.9.2, GFAffix v0.2.1, vg v1.62.0) and the untangle BEDs,
flank fasta and community table are all retained. What was missing is the *code*: the upstream
script that emitted `all_pedigrees_patches.tsv` was never committed, and the merging script
that IS in the repo, `scripts/pedigree/untangle_multimap_tracts.py`, hardcoded the same
hap1-everywhere assumption in its `PAIRINGS` table - so the labelling error had a second home.

FIXED: `--kid-correct-copies` selects PAN028's per-chromosome maternally inherited flank
(hap2 on chr1,3,4,5,7,9,12,14,15,18,21,22,X; hap1 otherwise). Re-running the repository's own
merge (both default and 1 kb-bridging settings):

| query copies | SELF-side tracts | KID-side tracts | total |
|---|---|---|---|
| published (hap1 everywhere) | 208 + 142 = 350 | **567** | 917 |
| **corrected per-chromosome** | 208 + 142 = 350 | **155** | 505 |

The KID-side call count falls 73% (567 -> 155) and the SELF side is unchanged, as it must be.
This independently reproduces the 74% wrong-label estimate from the earlier patch-level
recheck, now with the real merge logic rather than a per-interval approximation.
Outputs: `scripts/pedigree/recount/`.

### Threshold sensitivity: the retained call does not depend on the frequency threshold

FastGA `-f` sweep on the chr9q/chr3q terminal call (SELF.p.chr9q terminal 60 kb vs DAD's four
windows; coverage / mean identity over covered positions, `-s200 -c70 -i0.6`):

| f | chr3q h2 | chr3q h1 | chr9q h2 | chr9q **h1** |
|---|---|---|---|---|
| 10 | 99.0% / .9959 | 98.7% / .9956 | 94.4% / .9953 | **7.3%** / .9838 |
| 32 (published) | 100% / .9950 | 98.7% / .9960 | 97.7% / .9945 | **7.3%** / .9843 |
| 512 | 100% / .9950 | 98.7% / .9919 | 97.7% / .9944 | **7.3%** / .9843 |
| 8192 | 100% / .9979 | 98.7% / .9919 | 97.7% / .9945 | **7.5%** / .9868 |

The exclusion of the father's other chr9q copy and the full coverage of his chr3q copy hold
at every threshold, and identity is a dead heat throughout - which is why coverage, not
identity, is the discriminating statistic. Frequency 10 is noisier (94-99% coverage);
permissive seeding only cleans it up. Table: `scripts/pedigree/threshold_sensitivity_9q3q.tsv`.

Methods now also pins the aligner: sweepga v0.1.1 (commit 02e7ed0d) with its bundled FastGA,
and states the division of labour between the candidate-calling threshold (f32) and the
permissive validation settings (f512-8192, `-s200`).

### Module map promoted to Extended Data Fig. 3

`paper_prep/figures/ed3_module_maps/` (generator, vendored `module_tiling.tsv`, `sources.tsv`,
`caption.md`, README explaining the directory name), asset
`submission/fig/ExtendedDataFigures/ED_Fig3_module_maps.pdf`, referenced as `fig:ed3` from the
Discussion. NOTE: the legacy `paper_prep/figures/ed3/` directory holds an unrelated superseded
TAR1/telomere figure; it was accidentally overwritten during this work and restored from git
before continuing, hence the separate directory.

### Population-structure argument RESTORED (D-M6 addressed with a matched baseline)

It had been cut from the submission body in the V7 pass (AGENTS.md: popgen/FST removed from body
and Methods, "recover from git history if a reviewer asks"). D-M6 is that ask, and the analysis
machinery already existed: `scripts/popogen/matched_fst_d_m6.py` runs from local files and
compares subtelomeric Hudson F_ST with matched genome-wide continental values for the same five
superpopulation labels and the same metric (Bhatia 2013 Table 1; Patterson 2012):

| pair | subtelomeric F_ST (jk95) | matched genome-wide | ratio | verdict |
|---|---|---|---|---|
| AFR-AMR | +0.102 [+0.020,+0.184] | +0.071 | 1.44 | equivalent |
| AFR-EAS | +0.153 [+0.063,+0.242] | +0.144 | 1.06 | equivalent |
| AFR-EUR | +0.108 [-0.026,+0.242] | +0.150 | 0.72 | equivalent |
| AFR-SAS | +0.103 [-0.003,+0.210] | +0.110 | 0.94 | equivalent |

Non-AFR pairs are near zero or slightly negative, as genome-wide. Conclusion: subtelomeric
sharing carries ordinary continental structure, not a subtelomere-specific signature - the
compartment is population-stratified like the rest of the genome and exchange has not homogenised
it across populations. Restored as Results subsection "Communities carry ordinary continental
population structure" plus Methods "Population differentiation", with the remaining step stated
(cohort-internal window-sampled control from the phased VCF). Build: rev 16, exit 0, 0 undefined.

### Population differentiation is block-heterogeneous, not a gradient

Per-block AFR-vs-nonAFR F_ST spans -0.08 to +0.35 (s.d. 0.12 over 100 arm x pair estimates):
C11/chr6q homogenised below genome-wide (-0.08; AFR-EUR -0.25), C1/chr10q and C6/chr22q
undifferentiated (+0.02, +0.03), C3/chr16q and C15/chrXp as differentiated as any autosomal
region (+0.35, +0.28). The AFR-EUR/AFR-SAS direction is towards less differentiation but their
jackknife intervals include the matched baseline, so it is not read as a difference. And
differentiation does NOT track sequence sharing (r = +0.07 between block-mean within-community
Jaccard distance and AFR F_ST, n=10), so this is heterogeneity, not a gradient set by exchange
intensity. Added to the Results subsection; build rev 17, exit 0, 0 undefined.

## Pending

1. **Build + visual check of Fig5**: all four panels now exist (A pedigree with the naming,
   B zoom, C paternal ribbon, D pairwise). Run `bash compile.sh` and inspect the assembled
   figure; optional for B is a split y-axis so the 1.000 plateau and the 0.993 dip separate.
2. **Retire or supplement the D/E assets**: `Fig5D_maternal.png`,
   `Fig5E_pan028.png`, and the scripts that produce them
   (`submission/scripts/figures/make_fig5_panels.sh`, `make_fig5_zoom.R`,
   `make_fig5_pedigree_supporting_material.py`) — move to supporting material or
   delete, and drop their rows from
   `submission/supporting_material/Fig5_pedigree_supporting_material.tsv`.
3. **PAN028 label exposure (audited, see above)**: re-run the untangle for the 228
   wrong-label KID patches with the per-chromosome maternal copy; fix the 92% and
   "13 of 16 crossover-like" wording in `paper_prep/figures/fig4/caption.md`,
   `NATURE_DRAFT_v4.md`, `NATURE_DRAFT_v5.md`; check `scripts/pedigree/*.py` and
   report §14 for the same assumption.
4. **Acrocentric paragraph**: state the flagger breakdown for the acrocentric arms
   specifically (global: 66/538 patches overlap problematic regions, 14.2% of patch
   bp vs 6.7% flank background, 0/133 gene-conversion tracts) so the error-prone
   status is explicit.
5. **Margin table**: per-event identity margin + alternative-haplotype exclusion for
   every retained call, for Methods/supplement.
6. **flagger citation** missing from `bibliography.bib` (mentioned by name only) - add the reference.
7. **Build verification**: LaTeX is not installed on this machine; `paper.tex` changes
   are unverified by compilation. Run `bash compile.sh` on a machine with
   pdflatex+bibtex and check `grep -c undefined paper.log`.
7. **Discussion/abstract audit** for "recurrent exchange" phrasing that implies
   multiple pedigree events; keep the population-level recurrent-exchange claim
   (supported by the pangenome survey) separate from the single germline event.
