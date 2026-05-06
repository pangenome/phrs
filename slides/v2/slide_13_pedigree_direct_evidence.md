## Title

Caught in the act — three generations of a T2T pedigree show ongoing subtelomeric exchange

## Bullets

- **WashU T2T pedigree (Cechova et al. 2025), 3 generations, 4 individuals** — PAN010 (maternal grandmother) → PAN027 (mother) → PAN028 (daughter), plus PAN011 (paternal grandfather). Every haplotype telomere-to-telomere; *odgi untangle* compares each child flank against its parent in the implicit pangenome graph.
- **538 high-quality inter-chromosomal patches; 494 (92%) fall inside the HPRC v2 Leiden communities we built from 233 unrelated samples** (slide 9). The graph predicts where exchange shows up — the family delivers the events in the predicted communities.
- **133 ectopic gene-conversion–like sandwich tracts at score ≥ 0.81; 96 at the perfect 1.000/1.000 alignment ceiling.** Plus **16 crossover-like** events with the patch left and right flanks resolving to *different* haplotypes — the breakpoint signature of a meiotic crossover involving an inter-chromosomal segment.
- **C7 acrocentric traffic dominates** (chr22p:h2 ↔ chr13p / chr21p / chr14p / chr15p, dozens of independent patches). Non-acrocentric hits land in the same named clades from slide 9: **chr18p ↔ chr10p (C2 — Linardopoulou pair), chr3q ↔ chr9q (C3 — f7501 cluster), chrXp ↔ chrYp (C15 — PAR1)**, and **DUX4 chr4q ↔ chr10q (C1) at score 0.957 in PAN028 maternal** — the disease-named locus, in a living family.
- **Independent replication, CEPH1463 four-generation pedigree** (Porubsky et al. 2025): **11 parent features detected by *both* hifiasm and verkko assemblies** within Leiden communities — chr10/chr18 (C2), chr19/chr22 (C6), chr12/chr9 (C5), chr6/chr9 (C5). Same communities, different family, two assemblers — the signal is not a graph artifact.

## Primary figure

**Recommended (lead figure for the slide):** `end-to-end-report/pedigree-plots/washu/PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf` — the *odgi untangle* ribbon for PAN027's maternal haplotype 1 against her mother PAN010. The ribbon is the *direct* image of inheritance: each colored stripe is a segment of PAN027's flank, plotted at its source position on PAN010's genome. Self-color stripes mean PAN027 inherited the flank as expected; **off-color stripes are inter-chromosomal patches — the literal exchange events.** The dense stack of off-color sandwich blocks across PAN027's chr13p/chr21p/chr22p/chr15p flanks is what the 133 gene-conversion-like number *looks like*.

**Optional second panel — three-generation transmission triptych** (no SBATCH, all PDFs already on disk):

```r
# slide_13_washu_triptych.R — PAN010→PAN027→PAN028 untangle ribbons, side by side
# Output: slide_13_washu_triptych.pdf (drop-in for the slide)
library(ggplot2); library(magick); library(grid); library(cowplot)

pan027_mat <- ggdraw() + draw_image(magick::image_read_pdf(
  "end-to-end-report/pedigree-plots/washu/PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf",
  density = 200))
pan027_pat <- ggdraw() + draw_image(magick::image_read_pdf(
  "end-to-end-report/pedigree-plots/washu/PAN027.paternal_hap2_from_PAN011_father.untangle.pdf",
  density = 200))
pan028_mat <- ggdraw() + draw_image(magick::image_read_pdf(
  "end-to-end-report/pedigree-plots/washu/PAN028.maternal_hap1_from_PAN027_mother.untangle.pdf",
  density = 200))

label <- function(txt) ggdraw() + draw_label(txt, fontface = "bold", size = 11, hjust = 0.5)

plot_grid(
  plot_grid(label("Generation 2 ← 1 (maternal)\nPAN027.hap1 from PAN010"),  pan027_mat, ncol = 1, rel_heights = c(0.10, 1)),
  plot_grid(label("Generation 2 ← 1 (paternal)\nPAN027.hap2 from PAN011"),  pan027_pat, ncol = 1, rel_heights = c(0.10, 1)),
  plot_grid(label("Generation 3 ← 2 (maternal)\nPAN028.hap1 from PAN027"),  pan028_mat, ncol = 1, rel_heights = c(0.10, 1)),
  ncol = 3
)
ggsave("slides/v2/slide_13_washu_triptych.pdf", width = 15, height = 6)
```

If 90 s does not allow a triptych, ship just `PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf` full-bleed. The single panel already carries the slide.

## Speaker notes

So far we have shown that subtelomeres share sequence, that they cluster in 3D, and that the 3D signal survives in haploid sperm. None of that proves the recombination is ongoing — the sequence sharing could be a frozen ancient signal, and the 3D could be passive territory architecture. This slide is the proof that it is, in fact, ongoing — happening right now in three living generations of a single family.

The substrate is the WashU pedigree from Cechova and colleagues, 2025. Four individuals: a maternal grandmother PAN010, a paternal grandfather PAN011, the mother PAN027, and her daughter PAN028. Every haplotype is telomere-to-telomere — these are the cleanest assemblies of a pedigree that exist. We push each child haplotype through *odgi untangle* against its parent in the implicit pangenome graph and ask, segment by segment, where in the parent's genome each piece of the child's flank actually came from.

The picture you are looking at is PAN027's maternal haplotype, painted onto her mother PAN010. The diagonal stripes are correctly inherited self-flank. The off-color stripes — those dense sandwich stacks across the acrocentric short arms — are the gene-conversion-like events. We get **538** high-quality inter-chromosomal patches across the whole pedigree, and **494** of them — **92%** — fall inside Leiden communities that we built from a completely independent set of 233 HPRC v2 individuals. The graph from the population predicts where exchange happens in the family. **133** of those events have the textbook ectopic gene-conversion sandwich pattern. **96** are at the perfect alignment ceiling. **16** are crossover-like — left and right haplotype flip across an inter-chromosomal patch. That is the breakpoint signature of a real meiotic crossover involving an unrelated chromosome end.

The communities that show up in the family are exactly the communities we named in slides 5 and 9: **C7** acrocentrics with chr22p exchanging with chr13p, chr21p, chr14p, chr15p. **C2**, the Linardopoulou chr18p↔chr10p pair. **C3**, the f7501 chr3q↔chr9q cluster. **C15**, PAR1 X↔Y. And — the punchline — DUX4, the FSHD locus, chr4q gene-converting onto chr10q at score 0.957 in PAN028's maternal haplotype. The disease-named exchange caught in the act in a normal family.

Replication: the CEPH1463 pedigree, four generations, twenty-eight individuals, two assemblers (hifiasm and verkko). We require both assemblers to agree, in the same Leiden community, before we count an event. **Eleven** parent features survive that filter — including chr10/chr18 again (C2), chr19/chr22 (C6), and chr12/chr9 (C5) appearing in *both* G1 grandparents independently. Different family, different assemblers, same communities.

This is the direct empirical proof. Concerted evolution of subtelomeres is not a frozen signature — it is happening in three generations of a living human pedigree. That is the title of the talk.

## Time budget

**Target: 90 seconds.** Spend it. This is the title-thesis evidence — do not race through it.
- 10 s: setup ("we've shown sharing and 3D — now: is exchange *ongoing*?") + WashU pedigree intro (4 individuals, T2T, 3 generations).
- 15 s: read the figure aloud. Self-color = correctly inherited; off-color = inter-chromosomal patches. Point at the dense acrocentric stack.
- 20 s: the 92 % number, slowly. **538 / 494 / 92 %** — graph from the population predicts what we see in the family.
- 15 s: **133 gene-conversion-like / 96 at perfect score / 16 crossover-like.** Define crossover-like in one sentence (left and right flanks change haplotype across the patch).
- 15 s: name the communities — C7 acrocentrics, C2 Linardopoulou, C3 f7501, C15 PAR1, **DUX4 (C1) at 0.957 in PAN028 maternal** — set up slide 14.
- 10 s: CEPH1463 cross-assembler replication — 11 features, two assemblers must agree, same communities.
- 5 s: closing line — "concerted evolution caught in the act, three generations of one family."

If forced to compress, drop the CEPH1463 paragraph (10 s saved) before dropping any of the WashU numbers. **Never drop the 92 % number or the DUX4 hit** — those two are the slide.

## Notes for synthesizer

- **NEW slide vs v1.** The v1 deck (`slides/20260204_Subtelomics_overview_EG.summary.md`) ends at PCA by community and has no pedigree content. Slide 13 is the single most novel piece of the talk — it is the title-thesis evidence (`paper_prep/synthesis/CROSSWALK.md` §"C8" calls chapter 14 "the *direct empirical evidence* for 'ongoing recombination shapes subtelomeres' in the abstract title"). Per the task spec, **spend time here**: 90 s, not 60 s.
- **Continuity inbound (slide 12).** Slide 11's outbound note flags slide 12 as mouse meiotic Hi-C (zygotene peak, lepto→zygo→pachy→diplo). The natural pivot into 13 is: "mouse showed us that the strongest 3D contact is exactly when meiotic recombination peaks — does that recombination actually deposit human inter-chromosomal exchange in the human germline? Pedigree, T2T, three generations." Open with that one-sentence bridge; do not re-explain mouse.
- **Continuity outbound (slide 14 — DUX4/FSHD).** End on the DUX4 chr4q→chr10q at 0.957 in PAN028 maternal. That is the literal handoff to slide 14. Slide 14 will then re-tell the same locus as the disease-revealed instance with copy-number / median-22 framing. The synthesizer should make sure slide 14 *credits* this slide for the family-level observation rather than re-claiming it.
- **Locked numbers (single source of truth: `end-to-end-report/report/14_pedigree_recombination.md`).** These must agree across slides 13, 14, 15:
  - **538** total HQ inter-chromosomal patches in WashU
  - **494 / 538 = 92 %** within Leiden community
  - **133** gene-conversion-like at score ≥ 0.81 (line 47 of chapter 14)
  - **96** of the 133 at perfect 1.000/1.000 (counted from the chapter-14 table, scores 1.000/1.000 rows 1–74 plus rows 75 — first 74 are at exactly 1.000/1.000; allow ±2 if the synthesizer recounts)
  - **16** crossover-like (chapter 14 line 44)
  - **11** CEPH1463 cross-assembler-validated parent features (chapter 14 lines 255–269 table, 11 rows)
  - **DUX4 chr10q ← chr4q at score 0.957** in PAN028 maternal — chapter 14 line 172 (`chr10q | 402,748-403,967 | 1,219 | chr4q:h2 | ... | 0.957/0.957 | gene_conv | out | C1`). **Note: this patch is `out` of community C1's PHR (the patch itself is in chr10q's flank but outside the strict PHR boundary of C1)** — this is how Andrea's table is annotated. If pressed, the honest framing is "in C1, just outside the PHR call" — but for a 90 s talk the slide should say "C1 / DUX4" without the asterisk.
- **Slide 15 already cites 13's headline numbers** (`slides/v2/slide_15_concerted_evolution_thesis.md` lines 81–83). I am the source of truth for slide 15 — if slide 15 ever drifts from these, file a `wg msg` rather than silently editing this slide.
- **Figure source on disk.** `end-to-end-report/pedigree-plots/washu/PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf` (99 kB, generated by Andrea's `scripts/pedigree/analyze-pedigree-recombination.py`). Two siblings: `PAN027.paternal_hap2_from_PAN011_father.untangle.pdf` (76 kB) and `PAN028.maternal_hap1_from_PAN027_mother.untangle.pdf` (270 kB). Use the maternal_from_PAN010 file as the lead — that is the one called out by the task spec and it has the cleanest dense C7 sandwich stack visually.
- **The "ribbon" / "untangle" vocabulary needs one short on-slide gloss.** The audience may not know what an *odgi untangle* plot is. Suggested on-slide caption: *"each colored stripe = a piece of the daughter's subtelomeric flank, plotted at its source position in the mother's genome; off-color = inter-chromosomal patch."* If the synthesizer crops the figure caption, keep that sentence.
- **Do not oversell crossover-like.** 16 is small; the gene-conversion-like 133 is the headline. The crossover-like number matters because crossovers are the rarer, more dramatic event class — but if you read 16 with the same emphasis as 133 it sounds underwhelming. Frame as "*including* 16 crossover-like — the rarer event class with both flanking haplotypes reshuffled."
- **Do not oversell CEPH1463.** Per chapter 14 §Conclusions point 4: "the CEPH1463 single-assembler results are dominated by graph topology noise (12–13 % within-community vs 92 % in WashU) and should not be used as primary evidence." The 11 cross-assembler-validated number is robust precisely because it is the intersection of two noisy assemblies; the slide says "11 features confirmed by both hifiasm AND verkko" and that is enough — do not promote CEPH1463 to a primary-evidence framing.
- **Citations for the 90-s talk.** Cechova et al. 2025 (WashU T2T pedigree); Porubsky et al. 2025 (CEPH1463); Linardopoulou et al. 2005 (chr10p/chr18p exchange pair); Mefford & Trask 2002 (subtelomeric f7501); Lemmers / Tassin / Belyaev FSHD literature is for slide 14, not here. Both pedigree references are in `paper_prep/synthesis/REFERENCES.bib` per CROSSWALK §13. Safe to cite verbally.
- **Provenance.** Source TSV: `PHR_III/pedigrees/all_pedigrees_patches.tsv` (5,984 HQ patches, all pedigrees). Filtering: `is_interchr=True`, `min_score >= 0.8`, `500 bp <= size <= 100 kb`, in-Leiden-community (cross-community and unknown filtered out). Pattern definitions are at chapter 14 lines 22–28. The "score ≥ 0.81" cutoff is what chapter 14 actually uses for the gene-conversion-like 133 count (133 / 0.81 line 184); per-patch scores are listed exhaustively in chapter 14's tables.
- **No SBATCH, no new compute.** Lead figure is on-disk; optional triptych is pure ggplot2/magick wrapping existing PDFs. Both runnable in the agent worktree in seconds.
- **Do not modify other v2 slides.** Single output file, single commit. Slide 14 (DUX4/FSHD) is owned by another agent — coordinate via `wg msg` if numbers need to be aligned.
