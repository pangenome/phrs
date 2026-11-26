# Response to reviewers — assembled letter

This is the single assembled document. Source material: `REVIEWS.txt` (evidence record, rev 12),
`RESPONSE_LETTER_DRAFT.md` (pedigree), `RESPONSE_Q1_Q2_DRAFT.md` (multi-mapping, RNA proximity),
`RESPONSE_Q3_DRAFT.md` (arm-type and module architecture), `RESPONSE_STATS_AND_BIB.md`
(statistical controls, bibliography), `PEDIGREE_NAMING.md` (naming), `RESPONSE_CHANGELOG.md`
(change log). Manuscript wording is corrected in the text, not phrased as replies.

---

## Overview of what changed

Three things in this revision are substantive rather than stylistic.

1. **We found and fixed a labelling error in our own pedigree analysis.** The PAN028 assemblies
   name their two copies `hap1`/`hap2`; parent of origin is chromosome-specific, and we had used
   `hap1` everywhere as "maternal". On chr1, 3, 4, 5, 7, 9, 12, 14, 15, 18, 21, 22 and X that is
   the **paternal** copy, from a father whose genome is not in the data. Re-running our own
   merge code with the correct copy per chromosome reduces the child-side call count by 73%
   (567 → 155 tracts). The reviewers who suspected these were artefacts were right.
2. **We have withdrawn the last exchange candidate, and no exchange remains in the pedigree.**
   The single event we retained is an assembly artefact in one paternal haplotype plus a
   seed-frequency artefact at the second candidate (Point 4). Offering every candidate at once to
   IMPG similarity puts the rebuilt homologue first (dice 1.0000) and the alleged donor second
   (0.9990). Outside the acrocentric arms and the pseudoautosomal regions we claim no detectable
   exchange between non-homologous ends. The within-community statistic that summarises the
   pedigree is unaffected by this and remains stronger after the copy correction (0.945 corrected
   vs 0.918 published; null 0.84, p < 1e-4).
3. **We added an architectural result** that the first submission only implied: shared
   subtelomeric modules are orientation-consistent relative to the telomere (98.6% of 9.6 Gbp)
   while their order along the arm is not conserved, and community membership tracks repertoire
   rather than arrangement.

---

## Point 1 (Q1) — is the Hi-C within-community signal inflated by MAPQ0 placement?

No, and the reason is mechanical. The pipeline places each read pair at **one** position; it does
not spread multi-mapping mass over equally-scoring loci. Misplacement is therefore symmetric with
respect to community membership: a read placed at the wrong paralogue is as likely to land on a
between-community partner as a within-community partner, so placement adds noise without adding
directional enrichment. The existing falsification is the unique-sequence flanking control: in
HG002 Hi-C the between/within ratio is 0.027 at PHR coordinates and 0.002 in the flanking
unique-sequence windows at the same 50 kb resolution, a 13.5-fold strengthening when the
paralogous sequence is removed (10 kb: 0.029 vs 0.001, 29-fold).

There is also an a-fortiori argument. A PHR contains many near-identical copies but only a few
can be physically engaged in any given contact; sparse participation can *dilute* a signal, it
cannot *manufacture* one, and only the latter would rescue the artefact hypothesis.

We do not think a MAPQ-strict re-derivation is a valid control for this question, and we would
rather say so than offer it. MAPQ measures how well sequence discriminates between copies, not
how reliable a read is: inside a paralogous region MAPQ is zero *because* the copies are
indistinguishable, so filtering on MAPQ does not clean those windows, it deletes them, leaving
an empty matrix over precisely the windows under test. The direction of the distortion is
moreover known and bounded. Placement among equally-scoring loci is uniform, so a contact
between two specific loci is attenuated by roughly $1/k$ per end ($1/k^2$ for PHR-PHR contacts);
the attenuation is symmetric across community classes, so misplacement predicts *weaker*
enrichment inside duplicated windows and makes PHR-window enrichment a lower bound on the true
within-community preference. That is the direction we observe: the between/within ratio is 0.027
in PHR windows and 0.002 in the flanking unique-sequence windows, so enrichment is strongest
where multi-mapping cannot occur at all. The artefact hypothesis predicts the reverse.

What we can add, and have specified, is the test that does bear on this: because placement is
random, the draw itself can be repeated. If within-community enrichment is invariant across
independent placement seeds, the arbitrary choice carries no information about community
membership. Orthogonal contact assays in the same individual (ChIA-PET or Micro-C) would close
the question without relying on placement at all.

## Point 2 — the pedigree "92% within Leiden community" statistic

It now has a null, and the null is close to what the reviewer predicted. Ten thousand
permutations preserving the per-arm patch-count marginals give an expected within-community
fraction of 0.68 for the published patch set and **0.84** for the corrected one — not 0.5.
Against those nulls:

| set | observed | null mean | null 95% | p |
|---|---|---|---|---|
| published, 538 patches | 0.872 (469/538) | 0.677 | [0.664, 0.693] | < 1e-4 |
| corrected copies, 310 patches | **0.945** (293/310) | 0.840 | [0.832, 0.852] | < 1e-4 |

Wilson 95% intervals: 0.918 [0.892, 0.939] published; **0.945 [0.914, 0.965]** corrected. The
observation sits at the 100th percentile of both nulls.

## Point 3 — the chr5q/chr1p candidate "recurs, so it must be an artefact"

The recurrence is inheritance of a standing arrangement, and the candidate is not a
margin-supported exchange at all. On the maternal side the chr5q homologue explains the sequence
cleanly: 99.58% identity over the full terminal 60 kb versus 97.9% for the chr1p copies. The
chr1p-derived module is present in **both grandparents' own chr5q**, so it predates the family.
We now describe it as ancestral homology rather than as an exchange candidate, and we removed
our earlier "a different aligner did not reproduce it" argument, which was itself an artefact of
that aligner's default chunking filter. With the candidate no longer claimed, the objection
applies to nothing we assert.

## Point 4 — we have retracted the exchange, and the region is the opposite finding

Since the preprint we have re-examined the one candidate we retained, and it does not survive. We
withdraw it. Two independent artefacts produced the appearance of an exchange:

**(i) An assembly error in one paternal haplotype.** The terminal ~150 kb of the father's chr9q
haplotype 2 (PAN011#2#chr9 133,594,275-133,744,561) is misassembled in both releases; the release QC
marks a MISJOIN at the start of the block, and the assembled segment disagrees with the father's own
long reads — IMPG graph dice 0.9940 inside that tract against 1.0000 immediately outside it, where a rebuilt
assembly gives 1.0000 throughout.
Because his own chr9q copy was partly missing from the assembly, the closest *remaining* sequence, his
chr3q, won by default. The apparent interchromosomal donor was created by the missing haplotype.

**(ii) A seed-frequency artefact at the second candidate.** For chr1p/chr5q the father's chr5q
candidate never enters the alignment at the frequency threshold used for the whole-genome scan; at
threshold 48 and above the chr5q chain reaches the last base of its contig at identity 1.0000 and
matches at 48 diagnostic sites against 1 for chr1p.

**The definitive test offers every candidate at once.** Rather than comparing candidates one alignment
setting at a time, we placed the transmitted haplotype and every competing paternal end in a single
alignment — both the released and the rebuilt assembly of the candidate locus, the alternative
homologue, the competing non-homologue and a no-signal control — aligned them all-against-all with
sweepga at a frequency threshold high enough to filter nothing, and scored every pair with IMPG
similarity in 2 kb windows:

| candidate | graph dice to the daughter's paternal chr9q, across the tract |
|---|---|
| rebuilt chr9q haplotype 2 (father's own HiFi reads only, no daughter data) | **1.0000** |
| father's chr3q haplotype 2 (uniform inside and outside the tract) | **0.9990** |
| released (misassembled) chr9q haplotype 2 | 0.9940 inside, 1.0000 outside |
| father's chr9q haplotype 1 | 0.9930 |
| father's chr1p | 0.9940 |

The homologue wins outright. And the chr3q-to-homologue difference is 0.10% in graph dice — within the
range we measure at these ends generally, where the median chromosome end's closest non-homologous
partner is only 0.52% divergent and 11 of 45 ends lie within twice allelic divergence.

**This makes the region the opposite kind of evidence, and we now present it that way.** The father's
chr3q terminus is a *closer* match to his daughter's chr9q than his own other chr9q haplotype is
(0.9990 against 0.9930). A non-homologous end lying inside the spread of the homologous ones is
concerted evolution stated as a measurement, and it is now Figure 5. The pedigree figure has moved to
Extended Data, where it serves as a worked demonstration that identity-based donor assignment fails at
subtelomeric ends.

**What we can and cannot claim about exchange in a pedigree.** Across the whole three-generation family
the only heterologous sharing that behaves like obligate exchange is acrocentric-acrocentric and Xp/Yp —
the two compartments where non-homologous copies are as similar as alleles (paralog/allelic ratio 1.00
and 0.99). We therefore no longer claim detectable exchange between non-homologous ends anywhere else.
That is a narrower claim than the preprint's, and it is consistent with the mixing our divergence data
imply once it is expressed correctly. We no longer quote a rate of events: in an obligate region a
crossover occurs every generation, yet two copies there still coalesce ~62,500 generations ago, because
crossing over makes the two loci a single pool rather than a shallow genealogy. The data constrain a
ratio of genealogical depths, $R$, which is 1.0 where mixing is obligate (acrocentric short arms, PARs),
1.8 at the telomeric end of the compartment, 5.9 for the average block and 44 in the weakest blocks, so
mixing relative to the obligate case, $1/(R-1)$, runs from 1.2 down to 0.02.

## Point 5 — mouse ρ = 0.715 is a non-independent-pair Spearman

The pairs share chromosome ends, so they are not independent. We keep pair resolution and handle the dependence in the inference. Inference is at pair resolution (1,135 inter-chromosomal PHR pairs, 1 Mb flanks, 20 kb windows): Spearman ρ = 0.42, 0.61, 0.58 and 0.25 at leptotene, zygotene, pachytene and diplotene, with two-way clustered standard errors over the two arms of each pair (95% CI 0.37-0.47, 0.57-0.66, 0.53-0.62, 0.19-0.30) and a permutation null that shuffles similarity only among PHR pairs of the same arm pair (10,000 permutations; p = 1e-4 at leptotene, zygotene and pachytene, 1.5e-3 at diplotene). The association is elevated at zygotene and pachytene and lowest at diplotene. Zygotene minus pachytene is +0.038 at pair resolution (95% CI 0.009 to 0.067, p = 0.017), but a bootstrap over the 27 chromosome ends does not separate the two stages (95% CI -0.043 to +0.131), so the manuscript does not say "strongest at zygotene". The arm-level Mantel test (344 arm pairs, 10,000 row-and-column permutations; ρ = 0.40, 0.58, 0.63, 0.33, all p <= 1e-4) is kept in Methods as a display only: it discards the variance within arm pairs and puts the peak at pachytene. At 50 kb windows the four stages do not separate.

## Point 6 — the "100% bootstrap" on the NJ tree

Replaced with a genuine character-level bootstrap: 1,000 replicates resampling the 15,089
signal-bearing PHRs with replacement, recomputing the arm-level Jaccard matrix and rebuilding NJ
and UPGMA trees per replicate.

| clade | monophyletic in full tree | NJ | UPGMA |
|---|---|---|---|
| Xp/Yp (PAR1) | yes | 100.0% | 100.0% |
| Xq/Yq (PAR2) | yes | 100.0% | 99.9% |
| 10p/18p | yes | 100.0% | 100.0% |
| acrocentric p-arms | yes | 54.4% | 87.8% |
| TIGHT_q, 4q/10q | **no** | 0% | 0% |

The last row is a correction to our own language: those two groups are not monophyletic in the
arm-level tree, so the earlier "100% support" measured recovery of a *partition*, not monophyly
of a *clade*. They are now described as Leiden community assignments, and the 4q/10q relationship
rests on its direct sequence evidence (100% same-orientation, 86 Mb of shared sequence).

## Point 7 — missing confidence intervals

Every headline statistic now carries one: pedigree within-community (Wilson, both sets), F_ST
(block jackknife over arms), Mantel ρ (bootstrap/permutation), mouse arm-level Mantel
(permutation p < 1e-4 against a ±0.11 null range). No pedigree or 3D point estimate is reported
without an interval.

## Point 8 (Q3) — p-arm/q-arm asymmetry, now quantified

Within-community arm pairs across 41 arms and 15 communities (58 pairs), against 20,000 arm-type
permutations with communities fixed: P–P 19 (expected 14.8), Q–Q 22 (13.5, **p = 0.019**),
P–Q 17 (29.7, **p = 0.005**, depleted to 0.57×). The asymmetry is structural: **C6 is entirely
Q–Q (15 pairs, zero P–Q, p = 0.021)** and **C7 entirely P–P (10 pairs, zero P–Q, p = 0.048)**,
with cross-type mixing confined to C3 and C11; the two-arm communities fall where the concern
anticipated (C1 4q/10q is Q–Q, C2 10p/18p is P–P).

Orientation, which the concern also raised: measured per alignment block on the HPRCv2
all-vs-all flank maps with homolog pairs excluded, **98.6% of 9.6 Gbp of inter-arm alignment is
in the same orientation relative to the telomere** (P–P 97.8%, Q–Q 100.0%, P–Q 98.4%; per-pair
median 100%, only 7 of 764 pairs below 90%). Most P–Q pairs therefore align in the same
telomere-relative direction, so "P–Q occurs at opposite orientation" is not a general property;
the few net-inverted pairs are repeat-family architecture built from thousands of small blocks
rather than single inverted segments.

## Point 9 (Q2) — RNA-mediated proximity

We are keeping this, as a Discussion hypothesis, because it now has direct experimental
precedent rather than being speculation: Alu elements in enhancers and promoters pair through
complementary Alu-derived RNA and thereby determine which enhancer contacts which promoter
(Liang et al., *Nature* 619:727–735, 2023, doi:10.1038/s41586-023-06323-x). Our observation is
the same class of phenomenon at a different scale — an aggregate tendency of similar-sequence
subtelomeric regions to sit closer. An aggregate signal is compatible with only a small fraction
of repeats being engaged at any moment: sparse participation dilutes an aggregate, it does not
create one. We state it as a hypothesis with cited precedent rather than removing a documented
repeat-RNA proximity mechanism from the discussion.

The associated architectural result is new in this revision: modules shared between arm ends are
orientation-consistent relative to the telomere (98.6%) and carry a community-specific repertoire
(Jaccard 0.83 within versus 0.57 between, p = 5e-10; shared blocks 13.3 versus 9.9 kb median),
but their order along the arm is not conserved, within or between communities. That is the
pattern expected from repeated duplication plus orientation-biased conversion — gene conversion
and NAHR favour direct repeats, whereas inverted-repeat exchange yields dicentrics — rather than
from a conserved ancestral arrangement.

## Point 10 (Q2b) — generality across all segmental duplications

Not yet measured, and we say so. Closing it requires applying the same arm-level Mantel machinery
to all non-subtelomeric HPRCv2 segmental duplications (~50–100 Mb) and comparing that ρ with the
subtelomeric values (0.66 human; 0.58–0.63 mouse). It is on the roadmap as a separate analysis
rather than a revision item.

## Point 11 (M6) — F_ST matched-region control

Done, and it supports the reviewer's reading. Subtelomeric Hudson F_ST compared with matched
genome-wide continental values for the same five superpopulation labels and the same estimator
(Bhatia 2013 Table 1; Patterson 2012), block-jackknife intervals over ten arm/community blocks:

| pair | subtelomeric F_ST (95% jk) | matched genome-wide | ratio | verdict |
|---|---|---|---|---|
| AFR–AMR | +0.102 [+0.020, +0.184] | 0.071 | 1.44 | equivalent |
| AFR–EAS | +0.153 [+0.063, +0.242] | 0.144 | 1.06 | equivalent |
| AFR–EUR | +0.108 [-0.026, +0.242] | 0.150 | 0.72 | equivalent |
| AFR–SAS | +0.103 [-0.003, +0.210] | 0.110 | 0.94 | equivalent |

Subtelomeric differentiation is therefore not elevated above the genome-wide baseline, so there
is no subtelomere-specific population signature — the claim the reviewer objected to has been
removed from the manuscript and replaced by this comparison. The direction for AFR–EUR and
AFR–SAS is towards *less* differentiation than genome-wide (ratios 0.72 and 0.94), but the
jackknife intervals include the matched baseline, so we do not read that as a difference.

We also report the spread the aggregate conceals: per-block AFR-versus-non-African values run
from -0.08 (C11/chr6q, homogenised below genome-wide) to +0.35 (C3/chr16q, as differentiated as
any autosomal region; s.d. 0.12 over 100 block x pair estimates). Differentiation does not track
sequence sharing (r = +0.07 between block-mean within-community Jaccard distance and AFR F_ST),
so this is heterogeneity, not a gradient set by exchange intensity.

One methodological limit we state explicitly: local-ancestry inference is not validated in
duplicated subtelomeric sequence, so we deliberately do not build the argument on
ancestry-deconvolved genotypes. The comparison above uses continental superpopulation labels and
the standard Hudson estimator at arm resolution, which is the level these regions can support.

## Point 12 — bibliography

`Vollger2023` is not a duplicate: it is *Segmental duplications and their variation in a complete
human genome* (Science 376:eabj6965, 2023), distinct from `concerted_evolution_nahr_Vollger2023`
(*Increased mutation and gene conversion within human segmental duplications*, Nature 617:325,
2023). Smith 1976 (*Evolution of repeated DNA sequences by unequal crossover*) is present. Cech
2004 (*Beginning to understand the end of the chromosome*) was missing from the synthesis
bibliography and has been added; it was already present in the submission bibliography.

---

## Files and figures changed in this revision

- `submission/paper.tex` — pedigree section rewritten around the single retained event; chr1p/chr5q
  as ancestral homology; child-side panels retired; Methods (chromosome-copy assignment,
  assembly QC, per-event support table, statistical controls, threshold sensitivity); Discussion
  (module architecture); mouse analysis with pair-level clustered inference, arm-level Mantel as display.
- `submission/fig/MainFigures/` — Fig5A pedigree with the project naming, Fig5B zoom with
  `SELF.p`/`DAD.h2`, Fig5C paternal ribbon, **new Fig5D pairwise-alignment panel**; the maternal
  and grandchild ribbon panels are removed.
- `submission/supporting_material/Fig5_margin_table.tsv` — per-event support.
- `scripts/pedigree/untangle_multimap_tracts.py` — `--kid-correct-copies` (per-chromosome maternal
  copy); `scripts/pedigree/recount/` — the 567 → 155 recount.
- `scripts/hic/mantel_mouse_d_m5.py`; `scripts/cladistics/{char_bootstrap_d_m9.R,d_m9/,peerq3_*,
  arm_module_*}`; `scripts/ci/results_d_m12.json`; `scripts/pedigree/monte_carlo_null_d_m4.py`.
- `PEDIGREE_NAMING.md` — canonical naming (SELF/MOM/DAD/KID/PARTNER; copy notation; per-chromosome
  maternal table) and `AGENTS.md`/`CLAUDE.md` point to it.
- `RESPONSE_CHANGELOG.md` — item-by-item log with the numbers behind each change.
