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

## Pending

1. **Fig5C regeneration** (blocked: needs `inkscape`, not installed on this machine) -
   run `submission/scripts/figures/make_fig5_panels.sh`. Fig5B and Fig5D are done;
   `Fig5A_pedigree.pdf` should still print the nickname/copy scheme.
   Optional for Fig5B: split the y-axis so the 1.000 plateau and the 0.993 dip are
   visually separated (measured: same-chromosome 1.000 outside / 0.9933 inside the
   tract; chr3q 0.99988 / 0.99937).
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
6. **Build verification**: LaTeX is not installed on this machine; `paper.tex` changes
   are unverified by compilation. Run `bash compile.sh` on a machine with
   pdflatex+bibtex and check `grep -c undefined paper.log`.
7. **Discussion/abstract audit** for "recurrent exchange" phrasing that implies
   multiple pedigree events; keep the population-level recurrent-exchange claim
   (supported by the pangenome survey) separate from the single germline event.
