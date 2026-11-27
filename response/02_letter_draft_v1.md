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
2. **We retracted both pedigree exchange candidates.** The chr9q/chr3q candidate was an assembly error in one paternal
   haplotype, and the chr1p/chr5q candidate was an alignment artefact (Points 2-4). We corrected the assemblies of all four
   family members with each person's own reads (318 changes, Merqury QV up in all eight haplotypes) and redid the whole
   pedigree analysis on them. Offering every candidate at once to IMPG similarity puts the corrected homologue first
   (graph dice 1.0000) and the alleged donor second (0.9985), and the father's other chr9q haplotype does not carry the
   sequence at all. Across the three transmissions the only runs between non-homologous chromosome ends are PAR1 and
   acrocentric short arms, and 80 of the 82 runs join arms of the same Leiden community (null 0.929, p = 0.0003).
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

## Referee #1, major point 1: a base-level conversion scan (letter only; merge into the answer on concerted evolution)

The 2 kb scan in Methods cannot see conversion tracts shorter than a bin. We therefore also ran the same test at base resolution, and we report it here rather than in the manuscript.

**The test.** At each position we compared the haplotype with the other haplotypes of its own arm and with the copies of every other arm, using at least 20 copies of each.
- A donor-type site is a position where at least 90% of the copies of another arm share the haplotype's base and at most 10% of the copies of its own arm do.
- A tract is a run of at least two donor-type sites with no own-arm-type site between them.

**Three controls before counting.**
- **Clustering beyond chance.** Independent errors or recurrent mutations would scatter donor-type sites among the own-arm sites.
  - For each haplotype window and donor arm we computed the exact expected number of runs if the two kinds of site were in random order.
  - We kept a run only when fewer than 0.05 runs at least as long are expected by chance in that window.
  - Of 69,720 runs, 31,663 pass, with at most 554 expected by chance (false discovery rate at most 1.7%).
- **Assembly quality.** Before removal, 2.8% of these events overlapped an HMM-Flagger call of the HPRC release, against 0.66% of scored 2 kb bins. We removed the 205 events that overlap a call and the 24 events of samples without Flagger calls.
- **Aligner independence.** We realigned the copies of a random sample of events with minimap2 instead of wfmash and recomputed the sites with the same thresholds. In that sample, 42 of 48 events found in one individual and 26 of 29 shared events reappear (68 of 77, 88%, 95% CI 79-94%).

**Result.**
- **What remains.** After these filters 7,085 events remain, with median length 396 bp and 4 donor-type sites.
- **Individual and shared.** 1,508 are found in a single individual, and 5,518 are shared between individuals.
- **New relative to the 2 kb scan.** 4,739 lie outside the 2 kb tracts.
- **Read support.** HiFi reads support all 30 events that we inspected in a sample drawn before filtering: no read column inside any of them disagrees with the assembly.

**Limits.**
- **Sensitivity.** With tracts of known size planted into real haplotypes, the base-level scan found 7 of 18. The scan sees a conversion only where the donor differs from the replaced sequence at two or more sites, so these counts are a lower bound.
- **No dating of shared events.** For events shared by several individuals we cannot tell an old conversion from an ancestral state retained on some haplotypes without an outgroup, so we do not date them.

## Points 2-4: the pedigree (Referee #1, major points 2 and 3 and minor point 8; Referee #3, major point 8)

The referees were right about the pedigree. We retract both exchange candidates. The corrected analysis finds no exchange between non-homologous chromosome ends outside PAR1 and the acrocentric short arms. The genome-wide result, inter-chromosomal subtelomeric sharing organised into communities across 232 individuals, does not depend on the pedigree.

**What we did.**

1. **Assembly correction.** We corrected the version 1.1 assemblies of all four family members, each with that person's own reads only.
   - We reassembled each flagged region with hifiasm from the person's reads.
   - We applied a replacement only when more of the person's HiFi reads and more of their ONT reads supported it, and only when it did not raise the number of assembly 21-mers absent from the person's Illumina reads (Merqury).
   - We polished the replaced pieces with the person's HiFi reads.
   - We fixed two chromosome ends. We replaced the terminal 150,287 bp of the father's chr9q haplotype 2 with a contig that we assembled from his HiFi reads alone. We removed a 10,434 bp tail without telomeric repeats from the daughter's maternal chr22q.

   We made 318 changes. The Merqury quality value rose in all eight haplotypes (56.83-60.78 before, 57.18-60.85 after), and the error k-mers fell from 773,971 to 712,344. One script rebuilds the corrected assemblies from the released ones, base for base. The list of changes, with coordinates and read support, is in the code repository (Data availability).
2. **Transmission scan.** We aligned each transmitted haplotype to both haplotypes of the transmitting parent, on the corrected assemblies.
   - We cut the query into 2 Mb pieces, because FastGA stopped short of some chromosome ends with whole-chromosome queries.
   - We filtered the alignments to 18:18 mappings and classified 2 kb windows.
   - We drew runs of at least 10 kb at identity at least 0.95.
3. **Result.**
   - **Father to daughter:** the only runs between non-homologous chromosomes are PAR1 (two runs, 138 kb) and 47 runs between acrocentric short arms (830 kb). There are also two thin runs from the interior of chr1 and chr3 to acrocentric short arms.
   - **Mother to daughter:** 9 runs, all between acrocentric short arms (182 kb).
   - **Daughter to granddaughter:** 22 runs, all between acrocentric short arms (396 kb).
   - **Leiden communities:** 80 of the 82 runs join arms of the same Leiden community (0.976, Wilson 95% CI 0.915-0.993). A null that permutes the donor arms and keeps the per-arm run counts expects 0.929 (95% interval 0.927-0.951, p = 0.0003).

### Referee #1, major point 2: a new exchange or an inherited haplotype

All three of the referee's concerns were correct.

- **All three calls fell on one paternal haplotype.** The chr9q/chr3q call came from a misassembly of that haplotype.
  - In both releases the terminal ~150 kb of the father's chr9q haplotype 2 disagree with his own long reads, and the release QC marks a misjoin at the start of the block.
  - The released terminus carries a reverse-complement copy of a segment of his own chr1p: 100% of its canonical 21-mers fall within chr1p, collinear over 60 kb.
  - His own chr9q sequence was therefore missing at that terminus, and his chr3q won by default.
- **chr5q/chr1p was called in two transmissions.** This was an alignment artefact.
  - At the seed-frequency cutoff of the original scan, the father's chr5q homologue never entered the alignment, so chr1p won by default.
  - At threshold 48 the chr5q homologue reaches the last base of its contig at identity 1.0000. It matches the daughter's tract at 48 diagnostic sites against 1 for chr1p.
  - No chr1p/chr5q run remains on the corrected assemblies.
- **The margin was 1% or less.** We agree that a margin of this size cannot identify a donor at subtelomeric ends, where the competing sequences are near-identical paralogues. We no longer use the margin as evidence (Methods, "Why a margin alone cannot identify a donor here").

On the four requests:

- **Comparison against both parental haplotypes.** Every scan now aligns the transmitted haplotype to both haplotypes of the transmitting parent. At the former candidate we went further, because a call can be created or lost by which candidates reach the aligner.
  - We offered every competing paternal end at once, in one all-against-all alignment, with a seed-frequency threshold high enough to filter nothing.
  - We scored each candidate with IMPG graph dice over the 20 kb tract (Fig. 5d).

| candidate (corrected assemblies, plus the released chr9q) | graph dice to the daughter's paternal chr9q, across the tract |
|---|---|
| corrected chr9q haplotype 2 (father's own HiFi reads only, no data from the daughter) | **1.0000** |
| father's chr3q haplotype 2 | **0.9985** |
| released (misassembled) chr9q haplotype 2 | 0.9918 |
| father's chr9q haplotype 1 | no alignment (the 20 kb tract has no hit in the last 3 Mb of this haplotype) |

  The homologue wins. The difference between chr3q and the homologue is 0.15% in graph dice. That is inside the range we measure at chromosome ends generally, where the median end's closest non-homologous partner is 0.52% divergent (Extended Data Fig. 5e).

  The father's chr3q matches his daughter's chr9q at 0.9985, while his other chr9q haplotype does not carry the sequence at all. A non-homologous end that carries a block missing from one of the homologous ends is the concerted-evolution signal itself, and we now present it as Fig. 5d. The pedigree figure moved to Extended Data Fig. 7, as a worked example of an assembly error that mimics an exchange.
- **A genome-wide distribution of identity margins.** We no longer call a donor from a margin, and no candidate remains to locate on such a distribution.
  - We show the window-level best matches across the former tract, on the released and on the corrected assemblies (Extended Data Fig. 7B).
  - We list every run of the final scan with its mean same-chromosome and non-homologous identity (supporting_material/Fig5_pedigree_runs.tsv).
- **A false-positive estimate from permuted parent-child pairs.** [For the authors: not run. Possible control: the same scan on an unrelated pair, the mother against the father as if one were the other's parent, reporting its non-homologous runs.]
- **An unbroken contig for each tract.** No candidate tract remains. The corrected chr9q end is one contig assembled from the father's HiFi reads, joined to the assembly at the end of a 120,784 bp exact match.

  For the 82 runs of the final scan we also counted Merqury error k-mers (assembly 21-mers absent from the person's Illumina reads) on both sides of each run.
  - **The 80 runs outside PAR1:** 1 such k-mer on the transmitted side and 0 on the donor side, against 33.8 and 42.4 expected from each haplotype's genome-wide rate.
  - **The two PAR1 runs:** 26 on the transmitted side and 49 on the donor side, against 2.4 and 5.5 expected.

### Referee #1, major point 3: the exchange rate and Lemmers et al. (2010)

After the retraction the pedigree calls no exchange outside PAR1 and the acrocentric short arms. It therefore implies no rate that conflicts with the four discrete 4q/10q transfers that Lemmers et al. (2010) found across human populations.

On the request for a per-meiosis rate with its uncertainty: zero events in three transmissions give a Poisson 95% upper bound of 1.0 event per transmission (-ln 0.05 / 3). Three transmissions cannot constrain the rate more tightly than that.

The population data constrain exchange better, as a ratio of genealogical depths, R (Discussion and the Extended Data table on exchange regimes).
- In PAR1 a crossover is obligate in every male meiosis, yet two copies there still coalesce ~62,500 generations ago, because exchange makes the two loci a single pool.
- R is 1.00 where mixing is complete (acrocentric short arms, PAR1, PAR2). It is 1.82 at the telomeric end of the compartment, 5.93 for the average block and 44.0 in the weakest blocks.
- Mixing relative to complete mixing, 1/(R-1), therefore runs from 1.22 down to 0.023.

### Referee #1, minor point 8: the pedigree numbers

- **PAR1, 138 kb in the text against 92 kb in the figure.** PAR1 gives two runs: 92 kb (chrX 10,000-102,000) and 46 kb (chrX 106,000-152,000). The 138 kb is their sum, and the text, the Extended Data Fig. 7 legend and the run table now agree.
- **PAR1 against its ~2.7 Mb (major point 8).** We call a run only where the transmitted haplotype is closer to a non-homologous parental sequence than to the homologous one. In PAR1 the runs therefore mark where the daughter's paternal X matches her father's Yp better than his Xp, not the extent of PAR1.
- **The chr5q/chr1p lengths.** We withdrew the candidate and removed the panels that showed it.
- **The crossover-map argument and "the strongest non-acrocentric signal occurred in the paternal comparison".** We removed both.
- **The chr21p/chr4p run with a pericentromeric partner.** It disappeared with the old scan. The corrected scan has two runs of the same kind:
  - 18 kb from the chr3 centromere region to chr21p;
  - 10 kb from the chr1q pericentromeric region at 131 Mb to chr22p.

  We report them as interior runs, not as exchanges between chromosome ends. They are the satellite-type calls the referee anticipated: 2 of the 82 runs, 28 kb of 1,574 kb.

### Referee #3, major point 8: candidate exchanges

We agree, and the revision goes further. After the corrections no candidate remains. The section is now a direct test with a negative result ("Direct test for exchange between non-homologous ends"): this pedigree contains no evidence of exchange between non-homologous ends outside the acrocentric short arms and the pseudoautosomal regions. We no longer use the pedigree as evidence for recurrent genome-wide ectopic recombination. The evidence for exchange now comes from the population data (Results, "Non-homologous ends are homogenized, most strongly near the telomere").

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

Every headline statistic now carries one: pedigree within-community (Wilson interval and a permutation null), F_ST
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

- `submission/paper.tex`: pedigree section rewritten as a direct test with a negative result, both candidates retracted;
  Methods (assembly correction, transmission scan, competitive similarity test, runs and Leiden communities, Merqury
  counts in the runs, chromosome-copy assignment, statistical controls); Discussion (module architecture, exchange
  regimes as ratios of genealogical depths); mouse analysis with pair-level clustered inference.
- `submission/fig/MainFigures/`: Fig. 5 is now the concerted-evolution figure, with the competitive similarity test as
  Fig. 5d. The pedigree moved to Extended Data Fig. 7 (A pedigree, B released and corrected windows across the former
  tract, C paternal scan on the corrected assemblies). Extended Data Fig. 4 shows the pairwise alignments at the former
  candidate locus, with the released and the corrected chr9q.
- `submission/supporting_material/Fig5_pedigree_runs.tsv` — the 82 runs of the final pedigree scan (arms, communities, identities, Merqury error k-mers).
- `scripts/pedigree/fig5_v131_check/` — assembly correction and its rebuild (`g1_rebuild.sh`), the corrected-release tables and README (`results/corrected_release/`), the pedigree scan (`fig5_align_chunk.sbatch`, `fig5_similarity_s.sbatch`), the Fig. 5d competition (`fig5d_competition.sbatch`) and the run statistics (`pedigree_run_stats.py`).
- `scripts/pedigree/untangle_multimap_tracts.py` — `--kid-correct-copies` (per-chromosome maternal
  copy); `scripts/pedigree/recount/` — the 567 → 155 recount.
- `scripts/hic/mantel_mouse_d_m5.py`; `scripts/cladistics/{char_bootstrap_d_m9.R,d_m9/,peerq3_*,
  arm_module_*}`; `scripts/ci/results_d_m12.json`; `scripts/pedigree/monte_carlo_null_d_m4.py`.
- `PEDIGREE_NAMING.md` — canonical naming (SELF/MOM/DAD/KID/PARTNER; copy notation; per-chromosome
  maternal table) and `AGENTS.md`/`CLAUDE.md` point to it.
- `RESPONSE_CHANGELOG.md` — item-by-item log with the numbers behind each change.
