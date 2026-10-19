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
2. **What survives is narrower and better supported.** One ectopic exchange remains, on the
   chromosome copy where every required genome is sequenced, and the within-community statistic
   that summarises the pedigree is stronger after the correction (0.945 corrected vs 0.918
   published; null 0.84, p < 1e-4).
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
cannot *manufacture* one, and only the latter would rescue the artefact hypothesis. For closure
we have added the MAPQ-strict re-derivation as a specified experiment: contact matrices
re-derived from the upstream pair files behind the existing mcool objects with MAPQ0 pairs
dropped, reporting the strict-versus-random placement ratio (mcool alone cannot support this
filter, because it is post-binning and carries no per-read MAPQ).

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

## Point 4 — the one retained exchange, and why it is not an artefact

| evidence | value |
|---|---|
| recipient | SELF.p.chr9q (PAN027 paternal chr9q terminus, telomere-first 13–57 kb of the last 500 kb) |
| donor | DAD.h2.chr3q (PAN011 chr3q haplotype 2) — **99.79–99.84% identity, 100% of the terminal 60 kb covered** |
| excluded | DAD.h1.chr9q (**7.3–7.5% coverage at every threshold tested**) |
| structure | the father is heterozygous for the module (present on his chr9q hap2, absent from hap1), so the exchange predates his somatic genome and arose in the paternal lineage's germline |
| independence from thresholds | unchanged across FastGA frequency 10, 32, 512, 8192 |
| assembly QC | 1 of 538 patches overlaps a HiFi flagger call (acrocentric subclass 1 of 459); the larger ONT overlap is a generic chromosome-end effect present in every individual |

Two points deserve emphasis. The discriminating evidence is **coverage, not identity** — the
competing matches are near-identical paralogues, so identity is a dead heat at ~99.5% while the
excluded haplotype simply fails to match at all. And the frequency threshold we chose for
tractability (32) does not drive the call: we swept it from 10 to 8192 and the exclusion and the
coverage advantage are unchanged.

Within the pedigree, the child's chr9q terminus is then a 99.97% copy of that already-derived
paternal tip — an ordinary telomere-proximal crossover, not a second exchange.

## Point 5 — mouse ρ = 0.715 is a non-independent-pair Spearman

Replaced with an arm-level Mantel test with 10,000 row-and-column permutations: ρ = 0.41, 0.58,
0.63, 0.33 at leptotene, zygotene, pachytene and diplotene (20 kb windows, 344 inter-chromosomal
arm pairs, all p < 1e-4, permutation null centred on zero with a 95% range of ±0.11). The value is
attenuated relative to the per-pair statistic (0.715 → 0.58 at zygotene, 0.63 at pachytene),
which is the inflation the reviewer suspected. The stage profile is now reported as descriptive:
it rises into synapsis and peaks at pachytene, when the synaptonemal complex is fully assembled,
falling at diplotene; at 50 kb windows the profile is flat and we say so.

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

Acknowledged as not done. The "population-genetic signature" claim has been removed from the main
text, and the manuscript now states that cross-arm sequences carry population-genetic structure
consistent with genome-wide patterns rather than a subtelomere-specific one. The matched
non-subtelomeric control (equivalent-length autosomal windows per superpopulation pair, Hudson
F_ST, reported as a difference against that baseline) remains the analysis that would settle it.

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
  (module architecture); mouse analysis replaced by the arm-level Mantel.
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
