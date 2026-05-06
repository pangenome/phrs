# SLIDES v2 PLAN — Concerted evolution and unorthodox recombination of human subtelomeres

**Authors:** Andrea Guarracino, Erik Garrison
**Venue:** Biology of Genomes, Cold Spring Harbor — May 2026
**Status:** synthesised by `bog-v2-slides` (agent-809) on 2026-05-06 from 15 fanout-produced slide files in `slides/v2/`
**Anchor:** `paper_prep/synthesis/ABSTRACT.md`
**Substrate:** `paper_prep/synthesis/CROSSWALK.md`, `end-to-end-report/` (Andrea 14-chapter report)
**Companion deliverables (this synthesis):** `slides/v2/figure_manifest.md`, `slides/v2/coherence_check.md`

This file consolidates the 15 per-slide briefs into one talk-level plan. §1–§4 are the synthesis layer (time budget, narrative arc, open questions, risks). §5 reproduces each slide's content verbatim under one heading per slide so the speaker can read the deck end-to-end.

---

## 1. Time budget validation

Target: **≤ 900 s** (15-minute slot). Sum of per-slide budgets: **980 s** (16 min 20 s). **Overrun: 80 s.** Flagged.

| # | Slide | Budget (s) | Cumulative (s) | Mins:s | Notes |
|---|---|---:|---:|---|---|
| 01 | Title | 30 | 30 | 0:30 | hard floor — title + author + venue + thesis preview |
| 02 | Implicit interval tree | 50 | 80 | 1:20 | data-structure plumbing |
| 03 | IMPG workflow (ER connectivity) | 80 | 160 | 2:40 | methods anchor — "alignment IS the graph" + ER 230× argument |
| 04 | Genome-wide identity heatmap | 70 | 230 | 3:50 | first empirical result — PAR2-scale callout |
| 05 | Interchrom similarities (n-chrom-per-region) | 60 | 290 | 4:50 | PHR scale: 105 kb / 144 kb / n=15,668 |
| 06 | Length distributions per arm | 50 | 340 | 5:40 | clade outliers named on histograms |
| 07 | All-vs-all heatmap + NJ tree | 80 | 420 | 7:00 | the cladistic structure (NJ-tree-from artifact) |
| 08 | All-vs-all in 2D — chrom & superpop | 70 | 490 | 8:10 | population structure secondary |
| 09 | All-vs-all PCA — communities (keystone) | 80 | 570 | 9:30 | the keystone — every clade word maps to a cluster |
| 10 | Hi-C / Pore-C bulk + Mantel exclusions | 80 | 650 | 10:50 | sequence communities are 3D |
| 11 | Single-cell 3D — GM12878 + sperm | 60 | 710 | 11:50 | bulk → single-cell → haploid germline |
| 12 | Mouse meiotic — zygotene bouquet | 60 | 770 | 12:50 | meiotic-3D capstone |
| 13 | Pedigree direct evidence (WashU + CEPH1463) | 90 | 860 | 14:20 | **the proof** — 92% of patches in predicted communities |
| 14 | Gene biology aside (DUX4 / OR4F / TAR1) | 50 | 910 | 15:10 | **explicitly compressible to 0** |
| 15 | Concerted-evolution thesis (closer) | 70 | 980 | 16:20 | thesis pull-quote, locked text |

**Overrun triage** (per slide-author flags):

- **Drop slide 14 entirely** (slide-14 author calls this "compressible to 0"). Saves 50 s → 930 s total → still 30 s over.
- **Cut slide 13 CEPH1463 paragraph** (per slide-13 notes). Saves ~20 s → 910 s.
- **Cut slide 03 ER-callout discussion** (per slide-03 notes — though the 230× number must survive, since slide 15 calls it back). Saves ~15 s → 895 s.

**Recommended default:** drop slide 14 + soft-cut slide 13 CEPH1463 → 910 s ≈ 15:10. With a brisk delivery on slide 03 (~5 s saved) the deck fits inside the 15-min slot.

Numbers that **must not be dropped** under any compression (per-slide author flags + abstract anchoring):

- 230× Erdős-Rényi threshold (slide 03; called back on slide 15)
- 15,668 PHRs / 41 of 48 arms / median 105 kb (slide 05; called back on slide 15)
- PAR2 ≈ 334 kb (slide 04 anchor; called back across the deck)
- 92% of 538 inter-chromosomal patches within Leiden communities (slide 13; called back on slide 15)
- Thesis pull-quote on slide 15 (locked text from task spec)

---

## 2. Narrative arc — talk story flow

The talk has a four-act shape: **methods → empirical observation → mechanism → proof**, with biology-of-the-locus aside (slide 14) and synthesis closer (slide 15). The transition sentences below are the speaker's hand-offs between adjacent slides; they are derived from the "Continuity inbound / outbound" entries in each slide's notes.

### Act I — Setup (slides 01–03, 160 s)

**Slide 01 → 02** *"What we did is survey inter-chromosomal subtelomeric relationships at HPRC v2 scale — 466 near-complete haplotypes. The question is older than the data, but until now the data didn't exist. **Before the biology, one slide of plumbing.**"*

**Slide 02 → 03** *"Each pairwise alignment becomes an interval; one tree per sequence; the union is an interval forest, walked by transitive closure. **That forest *is* the implicit pangenome graph.** Now zoom out to the full pipeline — wfmash all-vs-all over 18,827 telomere-anchored flanks."*

**Slide 03 → 04** *"Twelve percent sampling, 230× above the Erdős-Rényi threshold for graph connectivity — closure from any subtelomere reaches the whole genome. **So what does that look like?**"*

### Act II — Empirical observation (slides 04–09, 410 s)

**Slide 04 → 05** *"PAR2-scale interchromosomal homology at nearly every subtelomere — 10s to 100s of kilobases. **How many chromosomes are mixing here?**"*

**Slide 05 → 06** *"15,668 PHRs across 41 of 48 arms, median 105 kilobases — a PAR2-class exchange landscape replicated at every chromosome end. **The distributions per arm tell us *which* exchange landscape it is.**"*

**Slide 06 → 07** *"The fat-right-tail outliers are the named clades — acrocentric short arms, PAR2, PAR1, DUX4, the Linardopoulou pair. **The shape already encodes the community partition you're about to see.**"*

**Slide 07 → 08** *"Three algorithms — Leiden, UPGMA, neighbor-joining — agree on the same six abstract-anchored clades. The clades are real. **Now zoom from the 41×41 arm picture to the 15,668 individual flanks.**"*

**Slide 08 → 09** *"On the left, points cluster by arm community. On the right, the same points by superpopulation — population structure is real but secondary. **Now the keystone — what *names* go on those clusters?**"*

**Slide 09 → 10** *"Every clade word in the abstract is a cluster on this plot. Everything that follows is evidence about the clusters you can already see here. **From sequence into the nucleus.**"*

### Act III — Mechanism (slides 10–12, 200 s)

**Slide 10 → 11** *"Sequence-defined communities are physical — Hi-C, Pore-C, CiFi all light up the diagonal at 50 kb, and the signal *strengthens* when known confounds are removed. **Bulk Hi-C is a population average — does the signal survive at single-cell resolution?**"*

**Slide 11 → 12** *"Same direction in GM12878 Dip-C, same direction (much stronger) in sperm. The negative control rules out chromosome-territory crowding. **If the haploid product still shows it, what about the cells where the recombination happens?**"*

**Slide 12 → 13** *"Mouse meiotic Hi-C across four prophase stages — Mantel ρ peaks at zygotene, exactly the bouquet stage when telomeres cluster at the nuclear envelope. **We see the contact at the bouquet stage; here is the recombination it produces, one generation later.**"*

### Act IV — Proof + biology aside + closer (slides 13–15, 210 s)

**Slide 13 → 14** *"494 of 538 inter-chromosomal patches in a three-generation T2T pedigree fall inside the Leiden communities we built from 233 unrelated samples — 92%. The graph predicts the family. DUX4 chr4q→chr10q at 0.957 in PAN028 maternal — caught in the act. **(Brief aside on the biology inside these regions.)**"*

**Slide 14 → 15** *"DUX4, OR4F, TAR1 — disease, decay, exchange machinery — all writing themselves into the same architecture. **Now the synthesis.**"*

**Slide 15 (closer):** Method (the implicit pangenome graph, 230×) → empirical (15,668 PHRs at PAR2 scale) → mechanism (Hi-C/Dip-C/mouse meiotic 3D) → proof (92% pedigree). Then the locked thesis pull-quote: *"Subtelomeres concertedly evolve through ongoing inter-chromosomal exchange — observable in pedigrees, predicted by 3D, recovered by an implicit pangenome graph across HPRC v2."*

### Callback discipline

Every slide in §5 has a "Notes for synthesizer" / "Cross-slide callbacks" block — these are the tendons that hold the talk together. The high-traffic callback nodes:

- **Slide 03's 230× ER number** is the methodological licence for "no chromosomal partitioning"; called back on slide 15 (thesis bullet 1).
- **Slide 05's PHR scale (105 kb / 144 kb / n=15,668)** is the central scale anchor; called back on slides 06, 09 (PCA layout), 15 (thesis bullet 2).
- **Slide 06–07–09's clade vocabulary** (C1, C2, C6, C7, C14, C15) is the noun set for the rest of the talk; slides 13 and 14 reuse it; slide 15 enumerates 4 of them in the thesis bullet 2.
- **Slide 10's bulk-3D headline** (B/W = 0.056 / Mantel ρ → 0.79 after exclusions) is the bulk anchor; refined per-cell on slide 11; refined meiotic on slide 12; recapped on slide 15.
- **Slide 13's 92% number** is the title-thesis evidence; called back on slide 15 (thesis bullet 4) and is the line the talk closes on before the pull-quote.

---

## 3. Open questions for Erik (aggregated from per-slide "Notes for synthesizer")

All items below are extracted verbatim or paraphrased from the 15 slide files. Each links back to its source slide. These are decisions only Erik (lead author) can make; the synthesizer cannot resolve them.

### 3.1 — Talk logistics

- **Date.** Slide 01 says "May 2026 — date: TBD". Confirm the actual BoG 2026 date before final render. (Slide 01.)
- **Affiliations.** Slide 01 author lists names without affiliations. If Erik wants UTHSC / pangenome.eu / etc. logos, add as a footer (not a bullet). (Slide 01.)
- **Title preference.** Slide 01 author **recommends the manuscript title** ("Concerted evolution and unorthodox recombination of human subtelomeres") over the v1 deck's "Inter-chromosomal subtelomeric relationships". The recommendation is to make the manuscript title the headline and demote the v1 phrasing to bullet 3. **If Erik prefers the v1 phrasing for accessibility**, swap bullet 3 ↔ headline. (Slide 01.)

### 3.2 — Methods / framing decisions

- **PCA vs MDS / PCoA terminology.** The artifact is `cmdscale(...)` (classical MDS, equivalent to PCoA on a Jaccard distance). The abstract and v1 deck call it "PCA". Two paths (slide 08 author): (a) **relabel only** — change abstract C6 + slides 08/09 to "MDS / PCoA"; lower friction. (b) **run a real PCA** on a haplotype × PHR-presence binary feature matrix (REWRITE_PLAN TASK-16); higher friction; biology unchanged. **Slide-08 recommends (a) for the talk.** Erik's call. See `coherence_check.md` §3 for the cascade to axes labels (slide 09). (Slides 08, 09.)
- **12% sampling rate provenance.** Slide 03 cites 12% as "fraction of C(n,2) pairs that pass wfmash's k-mer prefilter and reach full alignment." CROSSWALK §7b flags this as outstanding — the number must be computed from the actual on-disk PAFs before manuscript Methods. For the talk, the speaker can punt to Methods. (Slide 03.)
- **HG002 CiFi gap in ED5b.** Slide 10 author notes the no-acrocentric exclusion analysis covers 7 samples, not 8 (HG002 CiFi is missing). Slide says "7/7 sample × technology cells tested" — do not claim "8/8". Confirm with Erik whether to plug the gap or live with it. (Slide 10.)
- **"PCA" verbal use.** Slide 09 author says "for the talk this is fine — speaker should use 'PCA' verbally". Erik to confirm whether the speaker should commit to MDS / PCoA verbally instead, given §3.2 path (a). (Slide 09.)

### 3.3 — Visual / asset decisions

- **v1 deck PDF restoration.** Slides 02, 03, and 09 reference `slides/20260204_Subtelomics_overview_EG.pdf`, which is **not in this worktree**. Restore the file, render substitutes, or substitute alternates. See `coherence_check.md` §1. (Slides 02, 03, 09.)
- **Optional R-rendered overlays / composites.** Eight slides ship optional R scripts inline (slide 03 ER inset, slide 06 callout, slide 07 stitch, slide 09 legend, slide 10 composite, slide 11 composite, slide 12 trajectory, slide 13 triptych, slide 15 overlay). Two are **load-bearing** (slide 09 clade legend, slide 12 trajectory inset); six are polish. Erik's call: render all of them up-front, or rely on Keynote-level arrangement of the existing PDFs. See `figure_manifest.md` §5. (Slides 03, 06, 07, 09, 10, 11, 12, 13, 15.)
- **Slide 13 triptych vs single-panel.** Slide-13 author offers a 3-panel transmission triptych (`PAN027.maternal_hap1` + `PAN027.paternal_hap2` + `PAN028.maternal_hap1`) as optional second panel. If 90 s budget gets cut, ship the single-panel maternal-from-PAN010 figure. (Slide 13.)
- **Slide 14 explicit skip.** Slide-14 author flags the entire slide as "compressible to 0" if the talk runs long. **Default recommendation: drop slide 14** if any of the earlier slides over-runs. Erik to pre-commit. (Slide 14, see also §1 and §9 in `coherence_check.md`.)

### 3.4 — Numbers to lock down (single source of truth)

The following numbers appear on multiple slides with potential drift if anyone re-edits. Each has a designated source-of-truth file:

- **PHR scale:** 15,668 / 41 of 48 / median 105 kb / mean 144 kb / range 5–500 kb / PAR2 ≈ 334 kb. Source: `CROSSWALK.md` §C4 + `end-to-end-report/report/01_pipeline.md` §"Arm-level community detection". Used: slides 04, 05, 06, 09, 15.
- **Methods scale:** 18,827 telomere-anchored flanks; ~12% sampling; ER `p* ≈ 5.21×10⁻⁴`; ~230× threshold. Source: slide 03 (computed in-slide). Used: slides 03, 15.
- **Communities ↔ clades:** C1=DUX4 (4q,10q), C2=Linardopoulou (10p,18p), C6=q-arm clade (1q,13q,17q,19q,21q,22q), C7=acrocentric_p (13p,14p,15p,21p,22p), C14=PAR2 (Xq,Yq), C15=PAR1 (Xp,Yp; 18q n=1). Source: `CROSSWALK.md` §C5 + `end-to-end-report/report/01_pipeline.md`. Used: slides 06, 07, 09, 13, 14, 15.
- **Pedigree numbers:** 538 HQ patches / 494 within-community = 92% / 133 gene-conversion-like / 96 at perfect 1.000/1.000 / 16 crossover-like / 11 CEPH1463 cross-assembler-validated / DUX4 chr4q→chr10q at 0.957 in PAN028 maternal. Source: `end-to-end-report/report/14_pedigree_recombination.md`. Used: slides 13, 15.
- **Hi-C numbers:** B/W 0.056 / p=3.9e-85 (HG002 Pore-C); Mantel ρ before/after exclusions HG002 0.66→0.79, CHM13 0.66→0.85, HG02148 0.15→0.72; per-pair ρ=0.83 NA19036 / 0.81 HG02148 at 10 kb. Source: slide 10 + `end-to-end-report/report/05_hic_validation.md`. Used: slides 10, 15. **Note:** slide 15 currently cites a slightly different summary range (0.027–0.074) — see `coherence_check.md` §5.
- **Mouse meiotic:** Mantel ρ 0.687 / 0.718 / 0.683 / 0.577 (lepto/zygo/pachy/diplo, 50 kb, 1 Mb window). Per-PHR-pair Spearman ρ=0.715, p=4.4e-55, n=344 at zygotene. Source: `end-to-end-report/report/08_mouse.md` lines 87–113. Used: slide 12.

If any of these numbers drifts during a future edit, treat the source-of-truth file as canonical and update the slide(s) accordingly.

---

## 4. Risks and decisions for the lead author (Erik)

This section is the explicit "needs your call" list, distinct from §3 in that these are *risk* / *decision* items rather than open questions about content. Each names the slide(s), the risk, and the recommended action.

### 4.1 — High-priority blockers

| # | Risk / decision | Slide(s) | Recommendation |
|---|---|---|---|
| 1 | v1 deck PDF missing from worktree | 02, 03, 09 | Restore from backup before deck assembly; if not available, render substitutes per `coherence_check.md` §1. |
| 2 | Time budget overruns 15-min slot by 80 s | all | Pre-commit to dropping slide 14 + soft-cutting slide 13 CEPH1463 paragraph — see §1 above. |
| 3 | Slide 06 says "six clades" but lists only 5 (omits C6) | 06 | Re-word slide 06 bullet 3 — "round out the **five outliers visible on this panel**; the sixth clade C6 emerges in slide 09's PCA." See `coherence_check.md` §2. **Cannot be done in this synthesis pass — slide files are read-only per task spec.** Flagged for next round. |
| 4 | PCA vs MDS terminology inconsistent across slides 08/09 + abstract | 08, 09, abstract | Pick path A (relabel only) per slide 08 author; cascade to abstract C6 wording. See `coherence_check.md` §3 + §10. **Cannot be done in this synthesis pass.** Flagged for next round. |
| 5 | Slide 12 trajectory inset not yet rendered | 12 | Run `Rscript` from the slide-12 R block before deck assembly. Without it the zygotene-peak claim is a single number on faith. See `figure_manifest.md` row 12-trajectory. |

### 4.2 — Soft decisions / taste

| # | Decision | Recommendation |
|---|---|---|
| 6 | Whether to run all eight optional R composites or arrange in Keynote | Render all up-front (~5 min one-pass); Keynote arrangement is more flexible but accumulates manual error across 15 slides. |
| 7 | Affiliations on slide 01 | Footer logos, not bullets. |
| 8 | Slide 13 triptych or single panel | Single panel `PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf` is sufficient — only run the triptych composite if the deck has slack. |
| 9 | Title (manuscript vs v1 phrasing) | **Manuscript title** per slide-01 recommendation — the C8 thesis ("concerted evolution / unorthodox recombination") sets up the slide-15 callback. |
| 10 | CEPH1463 framing depth on slide 13 | Soft cut if compressing — primary evidence is WashU 92%; CEPH1463 is replication. |

### 4.3 — Manuscript-level (not blocking the talk)

| # | Decision | Recommendation |
|---|---|---|
| 11 | Whether to run a real PCA artifact for the manuscript | TASK-16 in REWRITE_PLAN. Talk version uses MDS/PCoA per §4.1 #4. |
| 12 | Compute exact 12% sampling rate from on-disk PAFs | TASK in CROSSWALK §7b. Talk version is fine with the round number. |
| 13 | Slide 14 in the SI vs not? | Talk-only; per slide-14 notes, do not promote any of the three vignettes (DUX4 / OR4F / TAR1) to load-bearing manuscript claims. |

### 4.4 — Things to confirm with Erik before render

- All numbers in §3.4 — re-verify against source-of-truth files if any chapter / TSV has been updated since the slide fanout.
- Pull-quote text on slide 15 ("Subtelomeres concertedly evolve through ongoing inter-chromosomal exchange — observable in pedigrees, predicted by 3D, recovered by an implicit pangenome graph across HPRC v2.") — slide-15 author flags this as locked from the task spec. Confirm.
- Speaker-delivery tone for slide 14 if kept — frame as "and the biology is interesting too" aside, skippable.

---

## 5. Per-slide content

The 15 slide files are reproduced verbatim below, one heading per slide. The slide files themselves remain authoritative; this section is a single-pane-of-glass copy for the speaker / synthesizer. Any divergence between this section and the file in `slides/v2/slide_NN_*.md` should be treated as a synthesis-time copy error — the file wins.

---

### slide_01_title

Source: `slides/v2/slide_01_title.md`

## Title
Concerted evolution and unorthodox recombination of human subtelomeres

## Bullets
- Andrea Guarracino, Erik Garrison
- Companion to HPRC v2 (Nature, in submission)
- Inter-chromosomal subtelomeric relationships at HPRC v2 scale — 466 near-complete haplotypes
- Biology of Genomes, Cold Spring Harbor — May 2026 (date: TBD)

## Primary figure
None — title slide. (No figure required by task spec.)

## Speaker notes
Title slide. Land the full title once, then translate it for a BoG audience that includes folks who don't work on subtelomeres: what we actually did is survey *inter-chromosomal* relationships between subtelomeres at population scale, and what we found is that this looks like concerted evolution — gene-conversion-like and crossover-like exchange between non-homologous chromosomes — happening more broadly than previously appreciated.

Frame the talk as the companion to HPRC v2: same data, new question. The HPRC v2 main paper hands us 466 near-complete haplotypes; we use them to re-examine a question genetics has parked since the 1990s — how related are the ends of different chromosomes? — and answer it at scale.

Co-authored with Andrea Guarracino. Manuscript in submission to Nature. Fifteen-minute talk; thirty-second slide; move quickly to motivation.

## Time budget
30s

## Notes for synthesizer
- **Title recommendation: use the manuscript title** ("Concerted evolution and unorthodox recombination of human subtelomeres") as the primary headline on slide 01, with the v1 phrasing ("inter-chromosomal subtelomeric relationships") demoted to a bullet that translates the manuscript title for the BoG audience. **Rationale:** (i) single source of truth — the brief identifies the manuscript title as canonical, and slide 02's notes already commit the deck to threading C2/C8 phrasing from `paper_prep/synthesis/ABSTRACT.md`; using the same title on slide 01 makes that thread legible from the start; (ii) the manuscript title carries the C8 thesis ("concerted evolution", "unorthodox recombination") which slide 14/15 (closing) should land — opening with it sets up a callback; (iii) the v1 title is not lost — it appears in bullet 3 as the speaker's own translation. If the synthesizer disagrees and prefers v1 title for accessibility, swap headline ↔ bullet 3 — the slide still works.
- **Date placeholder.** Bullet 4 says "May 2026 — date: TBD". Confirm the actual BoG 2026 date with Erik before final render. Project context (`paper_prep/synthesis/ABSTRACT.md`, `AUDIT_REPORT.md`) consistently says "BoG this week" as of 2026-05-06, so the talk is likely 2026-05-05 to 2026-05-09; substitute the exact date once known.
- **Cross-slide threading.** This slide sets up: (a) HPRC v2 framing (slide 02 follow-up), (b) the "concerted evolution / ongoing recombination" thesis the closing slide must land (per C8). Slide 02 introduces the implicit pangenome graph / IMPG; this slide should *not* introduce technical machinery.
- **Affiliations.** Bullets list the authors but do not include affiliations (UTHSC for Erik, presumably). If the synthesizer wants institutional logos / affiliation lines, add as a footer rather than a bullet — keep the bullet count ≤ 5 and the slide light.
- **No callbacks needed from earlier slides** (this is slide 01).
- **DO NOT** modify v1 deck files (`slides/20260204_Subtelomics_overview_EG.*`) — task spec is explicit. They are the historical baseline.

---

### slide_02_implicit_interval_tree

Source: `slides/v2/slide_02_implicit_interval_tree.md`

## Title
Implicit interval tree — the data structure under the implicit pangenome graph

## Bullets
- Each pairwise alignment becomes an interval `[start, end)` on its target sequence (CIGAR + target range)
- Li & Rong (2020) "implicit interval tree": store intervals in a sorted array; tree is implied by index, not pointers — O(log n) overlap queries with near-zero memory overhead
- Build one interval tree per target sequence; the union over all sequences is an **interval forest**
- That forest, queried by transitive closure, **is** our implicit pangenome graph — no graph is ever materialized explicitly
- Foundation for the rest of the talk: every plot you see comes from queries against this structure

## Primary figure
`slides/20260204_Subtelomics_overview_EG.pdf` (page 2 — kept verbatim from v1):
horizontal interval plot of 10 example alignment intervals; explicit-tree diagram on the left showing nodes labelled `[Start, End) Index, MaxEnd` across levels 0–3, with a highlighted node `[300, 320)` linked back to its alignment record; compact ordered-array form on the right showing the same intervals as the implicit representation (no pointers, just a sorted run with implied tree structure).

## Speaker notes
One slide of plumbing before the biology — the rest of the talk is just queries against this object.

An alignment is an interval on a target sequence: start, end, CIGAR. Li & Rong (2020) showed you can index a set of intervals as an *implicit* tree — a sorted array where tree topology is recovered from the index, not stored as pointers. O(log n) overlap queries, essentially zero memory beyond the intervals.

We build one tree per target sequence. The union is an interval forest. Walked by transitive closure, that forest **is** the implicit pangenome graph — the same object a `pggb` graph represents, never materialized. The implementation is **IMPG** (https://github.com/pangenome/impg, locally `~/impg`); `impg query -x` does the transitive lookup that powers everything downstream.

So when a later slide says "we queried the pangenome at chr18 q-arm," what's literally happening is `impg query` walking these trees. Keep this picture in your head — it's the substrate for every plot that follows.

## Time budget
50s

## Notes for synthesizer
- Slide 01 (title) sets up "we surveyed subtelomeres at HPRC v2 scale"; this slide answers *how is that even tractable* before slide 03 walks through the IMPG workflow figure. The transition from 01 → 02 → 03 is: motivation → data structure → pipeline.
- Slide 03 already shows the per-sequence interval-tree stack as part of the IMPG workflow; this slide deliberately *introduces* that visual vocabulary so 03 reads as "now zoom out to the full pipeline" rather than re-explaining trees.
- The "implicit pangenome graph" phrasing here matches **C2** of the canonical abstract (`paper_prep/synthesis/ABSTRACT.md`) — the speaker should land that exact phrase so it threads through the rest of the deck.
- Callback opportunity: when slide 04 (HPRC query) or slide 07 (all-vs-all heatmap) appears, the speaker can gesture back to "every cell is an `impg query`."
- Figure is unchanged from v1 — no new R/ggplot2 work needed. If a future revision wants a clearer ordered-array highlight, that is a minor v3 polish, not blocking for BoG.
- IMPG citation: Garrison et al., `pangenome/impg` (https://github.com/pangenome/impg). Underlying algorithm: Cordes, Li & Rong, "cgranges: a C/C++ library for fast interval overlap queries", 2020 (the implicit-interval-tree formulation).

---

### slide_03_impg_workflow

Source: `slides/v2/slide_03_impg_workflow.md`

## Title

**The implicit pangenome graph — wfmash all-vs-all *is* the graph (no GFA, no construction step)**

Subtitle: every PAF edge is a pairwise mapping; an interval forest per sequence is the index. IMPG (https://github.com/pangenome/impg) does transitive closure over it.

## Bullets

- **All-vs-all wfmash (`-p 95`) over n = 18,827 telomere-anchored 500 kb flanks → 18,827 PAFs.** Each PAF edge is one pairwise alignment; the union of edges *is* the implicit pangenome graph. No pggb, no GFA, no explicit graph construction step.
- **Index = interval forest, one implicit interval tree per sequence (Li & Rong 2020).** Built directly from PAF target intervals — `O(n log n)` build, `O(log n + k)` interval queries, no extra data structures on disk.
- **Query = `impg query -x` (transitive closure)** — chase chains of overlapping pairwise mappings outward from any seed interval to recover every reachable region across haplotypes. This is the operation; the graph is never materialized.
- **Sampling is dense, not sparse.** wfmash's k-mer prefilter evaluates ~12% of the C(n,2) ≈ 1.77×10⁸ pair space at full alignment cost. The Erdős-Rényi connectivity threshold for G(n, p) on n = 18,827 is **p\* = log(n)/n ≈ 5.21×10⁻⁴**; **12% is ~230× above p\***, so the random sub-graph is densely connected w.h.p. — transitive closure from any subtelomere reaches genome-wide. This is what licenses "no chromosomal partitioning."
- **Cite IMPG**: https://github.com/pangenome/impg (Garrison et al.). The implicit-graph framing is canonical and predates this work — we use it; we do not invent it.

## Primary figure

Re-use the v1 slide 3 visual stack: **wave (wfmash) → dotplot (all-vs-all alignment) → interval-tree forest, one tree per sequence**. Source: `slides/20260204_Subtelomics_overview_EG.pdf` page 3 (do NOT modify; extract / re-photograph for the v2 deck during synthesis).

Add one small callout panel at lower right — the Erdős-Rényi sanity check. R/ggplot2 generation script (light, no SBATCH; runs in the agent worktree):

```r
# slide_03_er_callout.R — ER connectivity threshold sanity check
# Output: slide_03_er_callout.pdf  (~3 in × 2 in, intended as inset on slide 3)
library(ggplot2)
n <- 18827
p_star  <- log(n) / n            # ER connectivity threshold ≈ 5.21e-4
p_obs   <- 0.12                  # wfmash k-mer-prefilter evaluation rate
ratio   <- p_obs / p_star        # ≈ 230×

df <- data.frame(
  label = c("ER threshold\np* = log(n)/n", "wfmash sampling\np ≈ 12%"),
  p     = c(p_star, p_obs),
  fill  = c("threshold", "observed")
)
df$label <- factor(df$label, levels = df$label)

ggplot(df, aes(label, p, fill = fill)) +
  geom_col(width = 0.55) +
  scale_y_log10(
    breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1),
    labels = c("1e-4","1e-3","1e-2","1e-1","1")
  ) +
  scale_fill_manual(values = c(threshold = "#888888", observed = "#1f77b4")) +
  geom_text(aes(label = sprintf("%.2g", p)), vjust = -0.5, size = 3.5) +
  annotate("text", x = 1.5, y = 0.4,
           label = sprintf("~%.0f× above threshold\n→ densely connected\n→ closure reaches\n   genome-wide", ratio),
           size = 3.2, hjust = 0.5) +
  labs(x = NULL, y = "edge probability  p   (log scale)",
       title = sprintf("n = %s flanks", format(n, big.mark = ","))) +
  guides(fill = "none") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"))

ggsave("slide_03_er_callout.pdf", width = 3.0, height = 2.0)
ggsave("slide_03_er_callout.png", width = 3.0, height = 2.0, dpi = 300)
```

## Speaker notes

This is the methods anchor — slow down here. When we say *implicit pangenome graph*, we mean it literally: we do **not** build a GFA. We run wfmash all-vs-all on 18,827 telomere-anchored 500 kb flanks; the resulting PAF set **is** the graph. Each PAF edge is one pairwise alignment. We index it as an interval forest — one implicit interval tree per sequence, in the Li-and-Rong 2020 sense from the previous slide. To query we use IMPG (github.com/pangenome/impg), `impg query -x`, which walks transitive closure: start anywhere, chase overlapping mappings outward, collect what's reachable. No graph constructed, nothing to break.

Is 12% sampling enough? The Erdős-Rényi connectivity threshold for n = 18,827 is p\* = log(n)/n ≈ 5.2×10⁻⁴. Twelve percent is **230× above that**, so the random sub-graph is densely connected w.h.p. — closure from any subtelomere reaches every other one. *That* is what licenses "no chromosomal partitioning." Everything downstream rides on it.

## Time budget

**Target: 80 seconds.** This is the methods anchor — do not rush. Roughly 25 s on "the alignment IS the graph, IMPG is the query interface" (frame-shift from v1), 35 s on the Erdős-Rényi argument with the n = 18,827, p\* ≈ 5.21e-4, 230× numbers spoken out loud, 20 s on cite-IMPG and segue into the next slide.

## Notes for synthesizer

- **Slide 02 sets up** the implicit interval tree data structure (Li & Rong 2020). Slide 03 reuses it — call back explicitly: "the data structure from the previous slide, one tree per sequence, gives us the index." Do not re-explain it.
- **Slide 03 sets up slide 04** (querying the pangenome → genome-wide identity heatmap). The link is: "now that we have transitive closure, here's what we see when we ask the graph what each subtelomere matches." Make sure slide 04's lead bullet starts from `impg query -x`, not from a fresh introduction of the alignment.
- **Frame-shift vs v1 is the load-bearing change.** The v1 slide 3 title was "IMPG: IMplicit Pangenome Graphs" and just walked the pipeline. The v2 title must say *this is the graph — there is no GFA*. If the synthesizer compresses, keep that clause; drop the algorithmic detail before dropping the framing.
- **The 12% / Erdős-Rényi callout is new in v2 and must survive compression.** It is the methodological justification for "no chromosomal partitioning" in the abstract (CROSSWALK §7b). Without it, the abstract's wording is asserted but not defended in the talk.
- **Cite line:** `https://github.com/pangenome/impg` should appear on the slide, not just in the speaker notes — credibility cue for an audience that will recognize Garrison et al. tooling.
- **If a referee in Q&A pushes on "12% — derived how?":** the answer is "from the on-disk PAFs — fraction of C(n,2) pairs that pass wfmash's k-mer prefilter and reach full alignment." Methods writer must compute this from the actual PAF set (CROSSWALK §7b explicitly flags this as outstanding); the speaker can punt to Methods.
- **Do NOT** swap in a pggb GFA visualization here — that belongs to a downstream similarity / community-detection slide if at all. Mixing the two undoes the entire frame-shift.
- **No SBATCH needed** for the inset; the R script is ~3 in × 2 in, runs in seconds locally inside the agent worktree.

---

### slide_04_genome_wide_identity

Source: `slides/v2/slide_04_genome_wide_identity.md`

## Title
Genome-wide identity heatmap — interchromosomal homology at PAR2 scale

## Bullets
- 466 HPRCv2 haplotypes vs CHM13, per-position max identity to any matching chromosome (100 kb windows).
- Most of each chromosome is silent; dense red bands appear where assemblies reach the telomeres.
- chr18 q-arm inset: tight subtelomeric block of >98% identity matches from many chromosomes.
- These blocks span **10s–100s of kb — comparable in scale to PAR2 (~334 kb on Xq/Yq)**.
- PAR2-scale pseudo-homology at nearly every subtelomere → motivates the all-vs-all view that follows.

## Primary figure
**Recommended:** `paper_prep/figures/fig1/figure_fig1.pdf`, **panel (a)** — genome-wide stacked identity heatmap with the chr18 q-arm subtelomeric inset already embedded. This is the polished version of the v1 slide-4 visual concept (1:1 substitute, no rebuild needed).

**Alternative if a single-chromosome zoom is preferred:** `identity_heatmap_chr18.zoom_last1mb.pdf` (the chr18q subtelomeric panel alone, larger and more legible than the inset).

**PAR2 callout overlay (new; for synthesizer):** add a small callout on top of panel (a) — preferred placement is over the chrX/chrY rows of the main heatmap (or in white space at the right edge of the inset). Suggested text:

> Extended interchrom homology spans 10s–100s kb at nearly all subtelomeres — *comparable in scale to PAR2* (~334 kb, chrXq/chrYq).

If the synthesizer wants to render the overlay programmatically rather than annotate in Keynote/PowerPoint, the following ggplot2 snippet wraps any imported panel-1a PNG with the callout. It needs no data and no SBATCH:

```r
# slide_04_par2_callout.R — overlay a PAR2-scale callout on panel 1a
library(ggplot2); library(png); library(grid); library(cowplot)

panel_a <- readPNG("paper_prep/figures/fig1/figure_fig1.png")  # or a cropped 1a export
bg      <- ggdraw() + draw_image(panel_a)

callout <- ggdraw() +
  draw_label(
    "Extended interchrom homology: 10s–100s kb\n(comparable to PAR2, ~334 kb on Xq/Yq)",
    x = 0.5, y = 0.5, hjust = 0.5, size = 11, fontface = "bold"
  ) +
  theme(plot.background = element_rect(fill = "#FFF7E6", colour = "black", linewidth = 0.6))

# place callout in upper-right of the panel; tune x/y/width/height to taste
ggdraw(bg) +
  draw_plot(callout, x = 0.62, y = 0.82, width = 0.36, height = 0.10)

ggsave("slides/v2/slide_04_panel_with_par2_callout.pdf", width = 12, height = 7)
```

## Speaker notes
This is the empirical foundation of the talk. For every 100 kb window across 466 HPRCv2 haplotypes we plot the maximum alignment identity to any *other* chromosome. Most of each chromosome is silent — alignments only hit the same chromosome, as expected. But at the telomeres, dense red bands appear wherever assemblies reach the chromosome end. These are the inter-chromosomal exchange blocks. The chr18 q-arm inset zooms in on one: a tight band at the very end where many other chromosomes match at over 98% identity, extending tens to hundreds of kilobases inward from the telomere. Here is the reframe — this scale is what PAR2 looks like on the sex chromosomes. PAR2 is about 334 kb. We are seeing PAR2-scale pseudo-homology at nearly every subtelomere. That dramatically expands the known scope of pseudohomologous regions in the human genome, which is the thesis the next slides will quantify.

## Time budget
70 seconds.

## Notes for synthesizer
- **Concept preserved from v1 slide 4** (genome-wide 100kb identity heatmap + chr18q inset). Visual is a clean 1:1 substitute with the publication panel.
- **One new asset needed:** the PAR2 callout text. It is the *only* new content vs v1; it lands the abstract's central reframing ("comparable in scale to PAR2"). If you skip the overlay, at minimum say it on screen as a subtitle bullet.
- **Continuity inbound (slide 03):** the previous slide ends on IMPG / all-vs-all alignment as the *method*. This slide is the first *result* — open with "what does that look like?" The dense bands at telomeres are the answer.
- **Continuity outbound (slide 05):** sets up the unique-chromosomes-per-region view (v1 slide 5 / Fig 1b). The natural pivot at the end is "*how many* chromosomes are mixing here?" — which is what slide 05 quantifies.
- **PAR2 number provenance:** 334 kb is the canonical PAR2 length on chrXq/chrYq (T2T-CHM13). It is referenced repeatedly in `paper_prep/synthesis/AUDIT_REPORT.md` (C4) and `paper_prep/synthesis/CROSSWALK.md` (C4 entry, "PAR2 anchor 334 kb"). Safe to cite verbally.
- **Do not relitigate the inset.** v1 used chr18q, abstract / Fig 1a both keep chr18q — keep chr18q.

---

### slide_05_interchrom_similarities

Source: `slides/v2/slide_05_interchrom_similarities.md`

## Title

Interchromosomal similarities — n-chromosomes per region (HPRCv2)

## Bullets

- Orange traces: number of unique chromosomes matching each 100 kb window across CHM13; CEN / PAR / PHR / XTR painted as background reference.
- Spikes are not confined to PARs — every subtelomere lights up, plus centromeres and the acrocentric short arms.
- **PHR scale: median 105 kb, mean 144 kb, range 5–500 kb** (15,668 PHR sequences across 41/48 chromosome arms).
- That places typical subtelomeric pseudohomology on the same length scale as PAR2 (~334 kb) — and present at nearly every chromosome end, not just X/Y.
- Read the plot as: a PAR2-class exchange landscape replicated 41 times across the human genome.

## Primary figure

`p_num_chromosomes_wide.pdf` (worktree root) — produced by `plot-impg-coverage.R` (`p_num_chromosomes_wide`, line ~321; saved line ~552). This is the v1 slide 5 figure preserved verbatim.

Cross-reference for synthesizer: `paper_prep/figures/fig1/figure_fig1.pdf` panel **1b** is the publication-quality version of the same view (genome-wide num-chromosome heatmap, 100 kb windows). v1 traces are easier to read at conference distance, so keep the v1 plot here; reserve fig1b for the manuscript.

Optional annotation overlay (if synthesizer wants the PHR scale rendered into the figure rather than spoken — drop into the existing `plot-impg-coverage.R` after `p_num_chromosomes_wide` is built, ~line 334):

```r
# PHR scale annotation — median 105 kb, mean 144 kb (CROSSWALK C4 / Andrea ch. 01)
p_num_chromosomes_wide_annot <- p_num_chromosomes_wide +
  annotate("segment", x = 0, xend = 144000,
           y = -2, yend = -2,
           colour = "#1b7a3a", linewidth = 1.2) +
  annotate("text", x = 144000, y = -2,
           label = "PHR scale: median 105 kb, mean 144 kb (n=15,668)",
           hjust = -0.05, vjust = 0.5, size = 3, colour = "#1b7a3a")
ggsave("p_num_chromosomes_wide_phr_annot.pdf",
       plot = p_num_chromosomes_wide_annot,
       width = 16, height = 9)
```

## Speaker notes

This is the same view I showed last time — number of unique chromosomes per 100 kb window across CHM13, with CEN, PAR, PHR, and XTR painted as background. Last time I asked you to notice that the orange spikes pile up at the chromosome ends; this time I want to put numbers on it. Across our 466 HPRCv2 haplotypes we recovered 15,668 pseudohomologous regions — PHRs — spanning 41 of 48 chromosome arms. Median length 105 kb, mean 144 kb, range 5 to 500 kb. To anchor that scale: PAR2 is about 334 kb. So a typical subtelomeric PHR is on the same order as PAR2 — and unlike PAR2, this is happening at nearly every chromosome end. The takeaway is one sentence: **extended pseudohomology at nearly all subtelomeres, comparable in scale to canonical pseudoautosomal regions, but replicated dozens of times across the genome.** That is the central observation the rest of the talk explains — what these communities look like, who is in them, and why they exist.

## Time budget

60 seconds.

## Notes for synthesizer

- **Callbacks/setup.** Slide 04 is the genome-wide identity heatmap (avg identity per matching chromosome, 100 kb windows) with the chr18q-arm inset. This slide flips the same data from "how identical?" to "how many chromosomes share it?" — keep the visual continuity (same 100 kb window grid, same chromosome ordering on the y-axis, same Mbp x-axis). The phrase "subtelomeric patterns are known qualitatively but need more precise quantification" is set up on slide 04; slide 05 delivers the first quantification (PHR scale).
- **Forward setup.** This slide hands off to slide 06 (length distributions of inter-chromosomal matches) and slide 07 (all-vs-all heatmap). The "median 105 kb, mean 144 kb" number stated here should match the histogram summary on slide 06 — please verify the synthesis pass keeps these consistent (single source of truth: CROSSWALK §C4, citing Andrea end-to-end-report ch. 01).
- **PAR2 framing is a deliberate abstract callback.** The ABSTRACT explicitly compares PHR scale to PAR2; CROSSWALK flags this as a framing gap (REWRITE_PLAN TASK-03 — a PAR2-vs-typical-subtelomere length panel for Fig 1 is planned). The talk version of that comparison lives here as a one-liner ("PAR2 ~334 kb"). If the synthesizer ends up adding a dedicated PAR2-comparison slide, this bullet can be slimmed.
- **Figure provenance.** Primary figure is the v1 deck figure (`p_num_chromosomes_wide.pdf`), not the manuscript Fig 1b. Per task: keep v1. Fig 1b is noted only as the manuscript companion so the synthesizer doesn't accidentally swap them.
- **Annotation overlay is optional.** I provided the R snippet so the scale bar can be rendered into the figure if the synthesizer prefers a visual annotation over the bullet/speech version. Default: keep the figure clean and let the speaker say the numbers.
- **Numbers to lock down.** 15,668 PHRs / 41 of 48 arms / median 105 kb / mean 144 kb / range 5–500 kb / PAR2 ≈ 334 kb. All five appear consistently in ABSTRACT, CROSSWALK §C4, and `framing_synthesis.md` (table comparing PARs vs PHR communities, line ~50).

---

### slide_06_length_distributions

Source: `slides/v2/slide_06_length_distributions.md`

## Title

Length distributions of inter-chromosomal matches per arm — the outliers are named clades

## Bullets

- Faceted histograms (kept verbatim from v1 slide 6): one panel per chromosome arm; **p-arm = blue, q-arm = orange**; pink fill marks the five introvert arms with no inter-chromosomal hits (2p, 3p, 5p, 8q, 11q, 14q).
- The bulk of arms cluster around the population scale stated on slide 05 — **median 105 kb, mean 144 kb, range 5–500 kb** (CROSSWALK §C4; n=15,668 PHR sequences across 41/48 arms).
- **The fat-right-tail outliers are not noise — they are the abstract's clades.** Acrocentric short arms (13p, 14p, 15p, 21p, 22p = **C7**) carry the heaviest, most uniform long-length distributions: fully homogenized rDNA-adjacent arms. PAR2 (Xq, Yq = **C14**) and PAR1 (Xp, Yp = **C15**) sit at the canonical pseudoautosomal scale (~334 kb on Xq/Yq). 4q and 10q (= **C1**, D4Z4 / DUX4) and the 10p–18p Linardopoulou pair (= **C2**) round out the six abstract-anchored clades.
- Pink panels are a biological signal, not missing data: arms without enough cross-chromosomal sharing to enter the 41×41 matrix.
- One-sentence reframe: **the shape of the per-arm length distribution recapitulates the community partition you are about to see on slide 09**.

## Primary figure

**Reuse v1 verbatim:** `/moosefs/guarracino/HPRCv2/PHR_III/plots/all-vs-all.1Mb.p95.id95.len_length_dist_by_chr_arm.pdf` (the exact figure on v1 slide 6 — faceted histograms by chromosome arm, p blue / q orange, pink for missing introvert arms). The underlying data is `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv` (18,827 sequences; columns `seq, arm, self_chr, region_start, region_end, chrs_involved, arms_involved`; PHR length per row = `region_end - region_start`).

**One new asset (this slide's deliverable):** a clade-callout overlay that names the four outlier groups directly on the facet grid. Renders to PDF / PNG locally with vanilla ggplot2 — the synthesizer can either drop the rendered overlay box onto the v1 figure in Keynote, or rebuild the histogram panel and let the overlay annotate live (data file is small enough to read directly with `readr::read_tsv`).

```r
# slide_06_clade_callouts.R
# Renders a callout legend that names the outlier facets on the v1 length-
# distribution facet grid. Pure ggplot2; no SBATCH; runs in seconds.
library(ggplot2)

cal <- data.frame(
  row   = 5:1,                                        # C1 lands at top of table
  C     = c("C7", "C14", "C15", "C1", "C2"),
  arms  = c("13p, 14p, 15p, 21p, 22p",
            "Xq, Yq",
            "Xp, Yp  (+18q, n=1)",
            "4q, 10q",
            "10p, 18p"),
  clade = c("acrocentric short arms — fully homogenized (rDNA-adjacent)",
            "PAR2 — pseudoautosomal q-end (~334 kb scale)",
            "PAR1 — pseudoautosomal p-end",
            "4q–10q DUX4 / D4Z4 — long-tail, copy-number diverse",
            "10p–18p — Linardopoulou 2005 pair"),
  fill  = c("#E5D8FA",   # C7 acrocentric (matches slide_09 palette)
            "#CDEAD3",   # C14 PAR2
            "#CDEAD3",   # C15 PAR1
            "#FDE2C8",   # C1 DUX4
            "#FFF1AA"),  # C2 Linardopoulou
  stringsAsFactors = FALSE
)

cols <- data.frame(
  x     = c(0.05, 0.22, 0.55),
  width = c(0.14, 0.30, 1.20),
  field = c("C", "arms", "clade"),
  stringsAsFactors = FALSE
)

cells <- do.call(rbind, lapply(seq_len(nrow(cols)), function(j) data.frame(
  x        = cols$x[j] + cols$width[j] / 2,
  width    = cols$width[j],
  y        = cal$row,
  text     = cal[[ cols$field[j] ]],
  fill     = cal$fill,
  fontface = ifelse(j == 1, "bold", "plain"),
  stringsAsFactors = FALSE
)))

hdr <- data.frame(
  x        = cols$x + cols$width / 2,
  width    = cols$width,
  y        = 6,
  text     = c("C", "Arms (outlier facets)", "Named clade — why the tail is long"),
  fill     = "#EEEEEE",
  fontface = "bold",
  stringsAsFactors = FALSE
)

title <- data.frame(
  x        = 0.85,
  width    = 1.6,
  y        = 7,
  text     = "The fat-right-tail facets are the abstract's clades",
  fill     = "#FFFFFF",
  fontface = "bold",
  stringsAsFactors = FALSE
)

note <- data.frame(
  x        = 0.85,
  width    = 1.6,
  y        = 0,
  text     = "Pink facets (2p, 3p, 5p, 8q, 11q, 14q) = introvert arms, no cross-chrom hits.",
  fill     = "#FFFFFF",
  fontface = "italic",
  stringsAsFactors = FALSE
)

p <- ggplot() +
  geom_tile(data = rbind(cells, hdr),
            aes(x = x, y = y, width = width, height = 0.95, fill = I(fill)),
            colour = "grey75") +
  geom_text(data = cells, aes(x = x - width/2 + 0.01, y = y, label = text,
                              fontface = fontface),
            hjust = 0, size = 3.0) +
  geom_text(data = hdr, aes(x = x - width/2 + 0.01, y = y, label = text,
                            fontface = fontface),
            hjust = 0, size = 3.2) +
  geom_text(data = title, aes(x = x - width/2 + 0.01, y = y, label = text,
                              fontface = fontface),
            hjust = 0, size = 3.6) +
  geom_text(data = note, aes(x = x - width/2 + 0.01, y = y, label = text,
                             fontface = fontface),
            hjust = 0, size = 2.6, colour = "grey35") +
  coord_cartesian(xlim = c(0, 1.7), ylim = c(-0.5, 7.5), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(6, 6, 6, 6))

ggsave("slide_06_clade_callouts.pdf", p, width = 8.0, height = 3.2)
ggsave("slide_06_clade_callouts.png", p, width = 8.0, height = 3.2, dpi = 200)
```

The five colored rows are exactly the abstract-anchored outlier clades (C1, C2, C7, C14, C15); palette matches slide_09 so the eye carries the same color-coding from histograms → community PCA. Place the rendered callout in the upper-right white space of the v1 facet grid (or beneath it as a banner) — this is the only on-slide change vs v1.

## Speaker notes

This is the same panel I showed last time — one histogram per chromosome arm, p-arms in blue, q-arms in orange, pink for the introvert arms with no cross-chromosomal hits at all. Most arms cluster around the population scale you saw on the previous slide: median around 105 kilobases, mean 144, range 5 to 500. What I want you to notice now is the **shape** of the outlier facets — the ones with the heaviest right tails, the most uniformly long distributions. Those are not noise, those are the clades the abstract names. The five acrocentric short arms — 13p, 14p, 15p, 21p, 22p — community 7 — fully homogenized, rDNA-adjacent. Xq with Yq — community 14 — PAR2, around 334 kb, the canonical pseudoautosomal scale. Xp with Yp — community 15 — PAR1. 4q with 10q — community 1 — the DUX4 / D4Z4 pair, long-tail and copy-number diverse. 10p with 18p — community 2 — Linardopoulou's 2005 pair. The pink facets — 2p, 3p, 5p, 8q, 11q, 14q — that absence is also a signal: those arms simply do not participate in the inter-chromosomal exchange landscape. So a single sentence carries this slide: **the shape of these per-arm distributions already recapitulates the community partition I will show you on slide 9.** The histograms know the clades before we even cluster.

## Time budget

50 seconds. Roughly: 10 s recap of the v1 visual ("same panel as last time, p blue / q orange, pink = introvert arms, scale anchored to the previous slide"); 30 s naming the five outlier clades on the callout (C7 acrocentric → C14 PAR2 → C15 PAR1 → C1 D4Z4 → C2 Linardopoulou, ~6 s each); 10 s landing the reframe ("the shape already recapitulates the community partition you'll see on slide 9").

## Notes for synthesizer

- **Concept preserved from v1 slide 6.** Faceted histograms by chromosome arm (p blue / q orange / pink for missing introvert arms). Do not rebuild the panel — the v1 PDF at `/moosefs/guarracino/HPRCv2/PHR_III/plots/all-vs-all.1Mb.p95.id95.len_length_dist_by_chr_arm.pdf` is what goes on the slide.
- **One new asset:** the clade-callout legend (R script above). It is the only addition vs v1. The script needs no data and no SBATCH; the rendered PDF / PNG sits over white space in the v1 facet grid (upper-right) or below it as a banner.
- **Continuity inbound (slide 05).** Slide 05 lands the population PHR scale (median 105 kb, mean 144 kb, range 5–500 kb, comparable to PAR2 ~334 kb) and the 15,668-PHRs-across-41-arms count. Slide 06 takes that same scale and shows it varies by arm in a *biologically structured* way — the outlier arms are the named clades. Lock these numbers to a single source of truth: CROSSWALK §C4, citing Andrea ch. 01.
- **Continuity outbound (slide 07 / 08 / 09).** This slide's reframe ("the shape already encodes the community partition") sets up the all-vs-all heatmap (slide 07) and the keystone PCA-by-community slide (slide 09). The five clade names spoken here (C1, C2, C7, C14, C15) MUST be the same five spoken on slide 09 — same colors, same arm lists. Synthesizer: keep the slide_06 and slide_09 callout palettes consistent (the R scripts already share the C1/C2/C7/C14/C15 fill colors).
- **Abstract-clade ↔ outlier-facet map (load-bearing — copy verbatim if compressing):**
  - C1 (4q, 10q) — DUX4 / D4Z4 long-tail; chapter 03 (28 DUX4L genes), chapter 04 (43.4 % type-discordance).
  - C2 (10p, 18p) — Linardopoulou 2005 Fig. 5 pair.
  - C7 (13p, 14p, 15p, 21p, 22p) — acrocentric short arms, fully homogenized (rDNA-adjacent). Heaviest, most uniform long-length distributions on the panel.
  - C14 (Xq, Yq) — PAR2, ~334 kb canonical scale.
  - C15 (Xp, Yp; 18q n=1) — PAR1.
  - Source: `paper_prep/synthesis/CROSSWALK.md` §C5; `end-to-end-report/report/01_pipeline.md` §"Arm-level community detection"; v1 slide 10 community list.
- **Missing-introvert callout (2p, 3p, 5p, 8q, 11q, 14q).** Same six arms appear on v1 slides 7, 8, 9, 10 in the upper-right corner. Keep the v1 phrasing here too — "introvert arms" is Erik's term and is consistent with "no inter-chromosomal hits" (their PHR length distribution is empty, hence the pink fill in v1).
- **Numbers to lock down across slides 05 / 06.** 15,668 PHRs / 41 of 48 arms / median 105 kb / mean 144 kb / range 5–500 kb / PAR2 ≈ 334 kb. All sourced from CROSSWALK §C4 (citing Andrea ch. 01) and `paper_prep/synthesis/ABSTRACT.md`. If a synthesis pass adjusts any of these, update slide 05 and slide 06 together.
- **Why the outliers are outliers (one-line each, for Q&A).** C7: rDNA recombination homogenizes the entire short arm; C14/C15: pseudoautosomal obligate crossover homogenizes both arms across X/Y; C1: D4Z4 / DUX4 macrosatellite copy-number polymorphism makes the right tail extend; C2: 10p–18p is the historical Linardopoulou exchange pair (long-known sequence-level homology block).
- **Do NOT touch the panel rendering.** The v1 PDF is on `/moosefs/guarracino/...` — out-of-worktree, but accessible read-only. Per task: keep v1. The slide deliverable is the *callout* on top of v1, not a rebuild.
- **If 50 s is cut to 35 s** (deck pacing): drop C2 (Linardopoulou) from the spoken walk — keep C7 acrocentric, C14 PAR2, C15 PAR1, C1 D4Z4. C2 is in the abstract clade list but is the smallest-effect outlier on the histograms and is recoverable on slide 09 if needed.
- **Cross-slide concern.** The v1 figure has small per-facet `n=` annotations (e.g., n=428, n=446); per the v1 summary these are sample counts per arm. Speaker should NOT verbalize per-facet counts at 50 s — they are visible on the figure for anyone who wants to read them. Save the per-arm n discussion for Q&A.

---

### slide_07_allvsall_heatmap_nj_clades

Source: `slides/v2/slide_07_allvsall_heatmap_nj_clades.md`

## Title

All-vs-all at the arm level — clustered heatmap + NJ tree with named clades (the cladistic structure)

## Bullets

- **Two-panel layout.** Left: 41×41 arm-level Jaccard distance heatmap (Fig 1c — `hprcv2.1Mb.subtelo.arm_dist_matrix.tsv`), arms ordered by Leiden k=15 community, cyan rectangles delimit the 15 communities, the original v1 UPGMA k=14 dendrogram on top recovers 14 of 15 communities (12/15 exact). Right: neighbor-joining tree on the **same** 41×41 distance matrix, rooted at the acrocentric MRCA, tip labels colored by the **six named clades from the abstract**.
- **The cladistic structure recovers expected pseudohomology and reveals novel clades** — every clade word in the abstract is a monophyletic block on the NJ tree, with **100% perturbation-bootstrap support** at each named clade's MRCA (1000 reps, σ = 25% of off-diagonal IQR):
  - **PAR1** (Xp/Yp) — red
  - **PAR2** (Xq/Yq) — blue (Xq–Yq edge length essentially zero, the shortest pair on the tree)
  - **acrocentric short arms** (13p, 14p, 15p, 21p, 22p) — green; chosen as the root, monophyletic
  - **10p–18p** (Linardopoulou 2005) — orange
  - **tight q-arm clade** 22q, 21q, 19q, 1q, 13q, 17q — purple (exact six-arm match to the abstract wording)
  - **4q–10q DUX4** — brown
- **NJ ↔ Leiden 1:1 mapping** (the partition is not an algorithm artefact): PAR1 = C15, PAR2 = C14, acrocentric_p = C7, 10p–18p = C2, tight q-arm = C6, DUX4 = C1.
- **Missing introvert arms** (no inter-chromosomal PHR detected, absent from this 41×41 matrix and from the NJ tree): **2p, 3p, 5p, 8q, 11q, 14q**. Same six arms drop out on every downstream similarity view (slides 08–09).
- **This is the highest-value v1→v2 swap.** v1 showed the heatmap with an unlabeled UPGMA dendrogram and a single "?" overlay; v2 adds the NJ tree the abstract names, with bold labels on every clade the audience will hear about for the rest of the talk.

## Primary figure

Two panels side-by-side.

- **Left panel (heatmap):** `paper_prep/figures/fig1/figure_fig1.pdf` — extract panel **(c)** only (arm-level 41×41 Jaccard distance heatmap with UPGMA k=14 dendrogram on top and the 15 Leiden community boxes overlaid). This is the published-quality version of v1 slide 7's heatmap; preserve the cyan community rectangles, the dendrogram, and the p/q-arm side annotation.
- **Right panel (NJ tree):** `paper_prep/figures/nj_tree_arms/nj_tree_annotated.pdf` (PNG companion alongside) — produced by the upstream `nj-tree-from` task in this worktree. 41 tips, rooted at the acrocentric MRCA, tip labels in bold for the six named clades and color-keyed to a legend with the clade names spelled out. Bootstrap support printed at every internal node.

If a side-by-side composite is wanted in the deck (no SBATCH; runs in agent worktree, ~2 s):

```r
# slide_07_heatmap_nj_combined.R — stitch the Fig 1c arm-level heatmap and the
# annotated NJ tree into a single landscape composite for the slide.
library(magick)

heatmap <- image_read_pdf("paper_prep/figures/fig1/figure_fig1.pdf",
                          pages = 1, density = 200)
# Crop panel (c) only — adjust to your fig1 layout if it changes.
hi <- image_info(heatmap)
panel_c <- image_crop(heatmap,
                      sprintf("%dx%d+%d+%d",
                              as.integer(hi$width  * 0.55),  # right ~55% width
                              as.integer(hi$height * 0.55),  # lower ~55% height
                              as.integer(hi$width  * 0.45),
                              as.integer(hi$height * 0.45)))

njtree <- image_read_pdf("paper_prep/figures/nj_tree_arms/nj_tree_annotated.pdf",
                         density = 200)

# Pad to equal height, then concatenate horizontally.
h <- max(image_info(panel_c)$height, image_info(njtree)$height)
panel_c <- image_extent(panel_c,
                        geometry_size_pixels(width = image_info(panel_c)$width,
                                             height = h),
                        gravity = "center", color = "white")
njtree  <- image_extent(njtree,
                        geometry_size_pixels(width = image_info(njtree)$width,
                                             height = h),
                        gravity = "center", color = "white")
combo <- image_append(c(panel_c, njtree))
image_write(combo, "slide_07_heatmap_nj_combined.png", format = "png")
```

Crop coords are heuristic; if the synthesizer prefers, render the heatmap panel directly from `paper_prep/figures/fig1/figure_fig1.R` (the panel-c block) and skip the PDF crop.

**Annotation overlay (recommended, on the NJ tree only)** — the upstream task already paints the six abstract clades in bold colored text with a top-right legend; no new ggplot work needed for the figure itself. The only on-slide annotation to add is a one-line caption beneath the pair: *"Same 41×41 arm-level Jaccard distance matrix; left = clustered heatmap (Leiden k=15 + UPGMA), right = NJ tree (rooted at acrocentric MRCA, 1000-rep bootstrap)."*

## Speaker notes

This is the all-vs-all picture at the arm level. The matrix is forty-one by forty-one — every chromosome arm we have a signal on, against every other arm. Cell color is Jaccard distance on the pangenome graph; the cyan boxes are the fifteen Leiden communities; the dendrogram on top is the original UPGMA, which already pulled out fourteen of the fifteen blocks.

What I want you to look at is the right panel. This is a neighbor-joining tree on the **same** distance matrix, rooted at the acrocentric short-arm clade because that clade is monophyletic and gives a stable orientation. Every clade name you read in the abstract — PAR1, PAR2, acrocentric short arms, ten-p with eighteen-p, the four-q ten-q DUX4 pair, and the tight q-arm clade of twenty-two-q, twenty-one-q, nineteen-q, one-q, thirteen-q, and seventeen-q — every one of those is a monophyletic block on this tree, in bold color, with one hundred percent bootstrap support at the named-clade root. We perturbed the distance matrix a thousand times at twenty-five percent of its off-diagonal IQR, rebuilt the tree, and these six clade roots never broke.

So **the cladistic structure recovers expected pseudohomology — the pseudoautosomals, the acrocentric short arms — and it reveals novel clades**: the ten-p eighteen-p pair which Linardopoulou drew in 2005, the four-q ten-q DUX4 pair, and most strikingly the six-arm q-arm clade that the abstract names. Same answer from Leiden, from UPGMA, and from neighbor-joining. The next several slides are evidence about clusters you can already see here.

## Time budget

**80 seconds.** Roughly: 15 s setting up the two-panel framing ("same forty-one by forty-one matrix, two clustering views"); 10 s on the heatmap (cyan boxes = 15 communities, UPGMA dendrogram agrees on 14 of 15); 35 s walking the six abstract-anchored clades on the NJ tree (PAR1, PAR2, acrocentric, tight q-arm, 10p–18p, 4q–10q DUX4 — about 5 s each, name and color); 10 s on bootstrap robustness ("100% support at every named clade, three algorithms agree"); 10 s landing the "expected + novel" closing line and segueing to slide 08.

## Notes for synthesizer

- **Highest-value v1→v2 swap.** The v1 deck had the heatmap with an unlabeled UPGMA dendrogram and a single hand-drawn "?" overlay on what we now call C1 / DUX4. v2 keeps the heatmap (it works) but **adds the NJ tree the abstract explicitly names** (`Cladistic analysis based on neighbor-joining trees…`), and labels the six abstract clades in bold colored type. This single swap closes the C5 framing gap from CROSSWALK §C5 ("abstract names NJ but no NJ exists") and resolves REWRITE_PLAN TASK-01 + TASK-13 for the talk.
- **Layout decision.** Two panels side-by-side. If only one panel fits at conference distance, prefer the **NJ tree** alone (it carries the abstract terminology) and demote the heatmap to a small inset. The heatmap is the v1 figure and is recognizable; the NJ tree is the new artifact and the named-clade evidence.
- **Figure provenance — strict.** Heatmap panel = Fig 1c from `paper_prep/figures/fig1/figure_fig1.pdf` (already in repo, manuscript-quality, do not rebuild). NJ tree = `paper_prep/figures/nj_tree_arms/nj_tree_annotated.pdf` (produced by the upstream `nj-tree-from` task in this branch series; commit `602a9d3`). Both PNG companions are present alongside the PDFs. **Do not** re-render either — the only optional artifact is the `magick` stitch above, which is a 2-second compositing step.
- **Cross-slide callbacks (load-bearing).**
  - **Inbound (slide 06)** — slide 06 ends on PHR length distributions; this slide picks up the same 15,668-PHR object and asks "how do the *arms* group?". Pre-cue verbally: "we just summarized PHRs by length; now we summarize arms by who they share with."
  - **Outbound (slide 08)** — slide 08 takes the **same** 15,668 × 15,668 *sequence-level* Jaccard matrix and projects it into 2D, colored by chromosome and superpop. Hand-off line for the speaker: "we just saw the forty-one by forty-one **arm** picture — next slide unfolds it to the fifteen-thousand-flank picture."
  - **Outbound (slide 09)** — slide 09 is the **keystone** that walks the same six abstract clades on a 2D PCA layout with arm-level points. The clade vocabulary (PAR1=C15, PAR2=C14, ACRO_p=C7, 10p–18p=C2, TIGHT_q=C6, DUX4=C1) **must match** between this slide and slide 09. The Leiden ↔ NJ ↔ abstract-clade map is the single source of truth for the talk.
- **NJ ↔ Leiden ↔ abstract-clade map (single source of truth — copy verbatim):**
  - **C15 ↔ NJ block PAR1 ↔ "Xp/Yp via the pseudoautosomal regions"** (red on tree)
  - **C14 ↔ NJ block PAR2 ↔ "Xq/Yq via the pseudoautosomal regions"** (blue on tree)
  - **C7 ↔ NJ block ACRO_p ↔ "acrocentric short arms"** (green on tree; root)
  - **C2 ↔ NJ block 10p–18p ↔ "10p–18p homology"** (orange on tree)
  - **C6 ↔ NJ block TIGHT_q ↔ "tightly linked clade involving 22q, 21q, 19q, 1q, 13q, and 17q"** (purple on tree; exact six-arm match to abstract wording)
  - **C1 ↔ NJ block DUX4 ↔ "DUX4-containing homology between 4q and 10q with wide copy number diversity"** (brown on tree)
  - Source: `paper_prep/synthesis/CROSSWALK.md` §C5; `paper_prep/figures/nj_tree_arms/README.md`; upstream `nj-tree-from` task log: "all 6 abstract clades recovered as monophyletic … 100% support at every named-clade MRCA."
- **Bootstrap caveat — explain only if asked.** Support is from a perturbation bootstrap (Gaussian noise on the distance matrix, σ = 25% of off-diagonal IQR, 1000 reps), not a Felsenstein column-bootstrap, because the input is a derived distance summary rather than an alignment. Deeper backbone edges show 32–90% support — the named clades are robust; the relative ordering of clades is not. Mention only if reviewer-style Q&A pushes; the abstract uses the result, not the support method.
- **Missing introvert arms (2p, 3p, 5p, 8q, 11q, 14q)** — same six arms as slides 08/09. Their absence is a biological signal (no inter-chromosomal PHR detected) and absent from the 41×41 matrix by construction. One-line corner annotation; do not spend speaker time unless asked. (Note: chr18q has n=1 sequence and is typically grouped into C15/PAR1 on the heatmap; if the v1 heatmap shows it, leave the v1 layout alone.)
- **Three-algorithm agreement is the robustness story.** Leiden k=15 (community detection, modularity-based), UPGMA k=14 (agglomerative, on the same matrix, Fig 1c dendrogram), and NJ (distance-based phylogenetics) all recover the same six abstract-anchored clades. Speaker should land this in one sentence — "three algorithms, same answer, the clades are real." This is the answer to "is the clustering an artefact of Leiden?".
- **Numbers to lock down (single source of truth).** 41 × 41 arm-level matrix / 15 Leiden communities / UPGMA k=14 (12/15 exact, 14/15 partial agreement) / NJ 41 tips / 1000 perturbation reps σ=0.0163 / 100% support at all six named-clade MRCAs / 6 missing introvert arms (2p, 3p, 5p, 8q, 11q, 14q). All consistent with `paper_prep/synthesis/CROSSWALK.md`, `paper_prep/figures/fig1/caption.md`, and `paper_prep/figures/nj_tree_arms/README.md`.
- **Do NOT** re-color the heatmap or recompute the partition. The Leiden partition and the UPGMA dendrogram in Fig 1c are the published Andrea version; the NJ tree is the *additional* view that closes the abstract-naming gap. Two views of one matrix is the point.

---

### slide_08_pca_chromosome_superpop

Source: `slides/v2/slide_08_pca_chromosome_superpop.md`

## Title

**All-vs-all in 2D — colored by chromosome (left) vs superpopulation (right)**

Subtitle: same projection of the 15,668-sequence Jaccard distance matrix; left panel asks *what arm is it from?*, right panel asks *what population?* — only the left answer holds up.

## Bullets

- Two-panel layout, **identical points**, only the coloring changes. Left: each point is one of 15,668 subtelomeric flanks colored by source chromosome (chr1–22, X, Y); shape encodes p / q arm. Right: same scatter recolored by 1KGP superpopulation (AFR / AMR / EAS / EUR / SAS).
- **Left panel: clusters are arm-community shaped.** Points group by Leiden community (15 communities over 41 arms — D4Z4 chr4_q+chr10_q, acrocentric p-arms chr13/14/15/21/22 p, PAR1 chrX_p+chrY_p, PAR2 chrX_q+chrY_q, the 6-arm q-arm community chr1_q/13_q/17_q/19_q/21_q/22_q, etc.). Visible structure ≈ arm-community structure.
- **Right panel: superpopulations are mixed across all clusters.** The arm-community clusters do not split by AFR/AMR/EAS/EUR/SAS — population structure is **real but secondary** to arm-community structure.
- **Population structure *is* there, just smaller.** Hudson Fst on cross-arm affinity (chapter 04): mean **Fst = 0.044**; AFR vs non-AFR pairs **0.10–0.15**, non-AFR–non-AFR pairs **−0.05 to 0.02**. The AFR-deepest split mirrors the human out-of-Africa tree (chapter 12 novel contribution #19) — a population signal exists, but it lives at finer scale than the dominant arm-community signal driving the global 2D layout.
- **Missing introvert arms** (no inter-chromosomal PHR detected, excluded from this projection): **2p, 3p, 5p, 8q, 11q, 14q** — same six arms across both panels.

## Primary figure

Two-panel side-by-side. **Both panels are existing artifacts in Andrea's pipeline output:**

- Left panel — by chromosome: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-chromosome.pdf` (PNG companion: `hprcv2.1Mb.subtelo.mds.color-by-chromosome.png`).
- Right panel — by superpopulation: `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.color-by-superpopulation.pdf` (PNG companion: `hprcv2.1Mb.subtelo.mds.color-by-superpopulation.png`).

**Important — what these are.** Both produced by `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R` line 556 — `cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)` on the 15,668 × 15,668 Jaccard distance matrix from `pggb -p 95` + `odgi similarity --all -P`. Variance-explained labels on the axes come from `fit_full$eig / sum(abs(fit_full$eig)) * 100` (script line 579). This is **classical MDS (= PCoA on a Jaccard distance matrix)**, not PCA on a feature matrix — see Notes for synthesizer.

If a side-by-side composite is needed for the slide deck (no SBATCH; runs in agent worktree):

```r
# slide_08_pca_combined.R — side-by-side composite of color-by-chromosome
# and color-by-superpopulation MDS panels. Inputs are Andrea's existing
# panel PDFs; this script just stitches them.
library(magick)
library(grid)
library(gridExtra)

base <- "/moosefs/guarracino/HPRCv2/PHR_III/similarity"
left  <- image_read_pdf(file.path(base, "hprcv2.1Mb.subtelo.mds.color-by-chromosome.pdf"),       density = 200)
right <- image_read_pdf(file.path(base, "hprcv2.1Mb.subtelo.mds.color-by-superpopulation.pdf"),  density = 200)

# pad to equal height, then concatenate horizontally
h     <- max(image_info(left)$height, image_info(right)$height)
left  <- image_extent(left,  geometry_size_pixels(width = image_info(left)$width,  height = h),
                     gravity = "center", color = "white")
right <- image_extent(right, geometry_size_pixels(width = image_info(right)$width, height = h),
                     gravity = "center", color = "white")
combo <- image_append(c(left, right))
image_write(combo, "slide_08_pca_combined.png", format = "png")
```

**Annotation overlay (recommended)** — drop a single shared caption beneath the two panels: *"Same MDS / PCoA projection (cmdscale on 1 − Jaccard, k = 5); n = 15,668 flanks across 41 arms; missing introvert arms: 2p, 3p, 5p, 8q, 11q, 14q."* Speaker delivers the "left = arm-community, right = mixed → population structure is secondary" beat verbally.

## Speaker notes

This is the same all-vs-all object as the heatmap two slides back, just projected into two dimensions so we can ask two questions of one picture. Each point is one of 15,668 subtelomeric flanks. The projection is classical MDS — PCoA on the Jaccard distance matrix from the pangenome graph. **One scatter, two colorings.**

On the left I've colored by source chromosome. The clusters you see are arm-community-shaped: D4Z4 pulls chr4q and chr10q together; the acrocentric p-arms — chr13p, 14p, 15p, 21p, 22p — collapse into one cluster; PAR1 puts Xp on top of Yp; the six-arm q-arm clade with chr1q, 13q, 17q, 19q, 21q, 22q sits as its own cloud. The structure on the left panel is the arm-community structure from slide 07's heatmap, replotted.

On the right I've taken the *same* points and recolored them by 1000 Genomes superpopulation — AFR, AMR, EAS, EUR, SAS. The colors mix across every cluster. The clusters do **not** resolve by population. So the headline is: in this 2D view, what dominates is *which chromosome arm a sequence comes from*, not *which population the haplotype is from*.

That doesn't mean population structure is absent — it means it's secondary. Chapter 4 of our analysis tests it directly with Hudson's Fst on cross-arm affinity: mean Fst is 0.044, AFR-versus-non-AFR pairs sit at 0.10 to 0.15, and non-AFR pairs sit at zero. Same pattern you'd recognize from any human-population genetics study — AFR is the deepest split, mirrors out-of-Africa. So the population signal is present and quantifiable; it just doesn't dominate the global 2D layout, because the arm-community signal is much stronger. **Two-tier hierarchy: arm-community first, population structure within.**

## Time budget

**Target: 70 seconds.** Roughly 15 s setting up the two-panel framing ("same plot, two colorings, different questions"); 25 s on the left panel — call out 3–4 specific clusters by name (D4Z4 4q/10q, acrocentric p, PAR1, q-arm clade); 20 s on the right panel mixing → "population structure is real but secondary"; 10 s landing the Fst numbers (0.044 mean, 0.10–0.15 AFR vs non-AFR, out-of-Africa shape).

## Notes for synthesizer

- **CRITICAL — PCA vs MDS labeling (CROSSWALK §3 C6).** The v1 deck slides 8 and 9 both said "PCA". The artifact is **MDS / PCoA**, computed by `cmdscale(as.dist(jaccard_dist_df), eig = TRUE, k = 5)` at `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R` line 556, output cached to `hprcv2.1Mb.subtelo.full_mds.rds`. There is **no `*pca*.rds` in the HPRCv2 tree** (verified by `find /moosefs/guarracino/HPRCv2 -name "*pca*"` — only `HG002.GRCh38_no_alt.dipcall.bed` matches, unrelated). I have re-titled this slide accordingly. Two paths going forward; the synthesizer must pick one and propagate:
    - **(a) Relabel only.** Change abstract wording from "Principal component and community detection analyses" → "MDS / PCoA and community detection analyses on the Jaccard distance matrix". Erik's clarification per CROSSWALK is "abstract is not locked"; this is the lower-friction option and it is what the existing artifact actually shows.
    - **(b) Run a proper PCA** on a haplotype × PHR-presence-binary feature matrix (REWRITE_PLAN TASK-16). This produces a strict PCA but adds work and a new artifact; the *biological* conclusion (arm-community first, population secondary) does not change.
    - **My recommendation: (a)** for the talk this week, with TASK-16 as the manuscript follow-up. The talk slide says "MDS / PCoA"; if the abstract is later rewritten to PCA, the slide can be relabeled.
- **The v1 axis labels were "Dimension 1 (16.05%)" and "Dimension 2 (11.2%)".** Those numbers come from `var_explained_full[1]` / `var_explained_full[2]` in the script (line 579–581) — i.e., from `fit_full$eig / sum(abs(fit_full$eig)) * 100`. They are the **MDS** variance-explained, not PCA variance-explained. They are real numbers from the actual run; keep the values, just relabel the projection.
- **Combining v1 slides 8 + 9 into one combined slide is the right move** — they share axes, share points, share the missing-introvert-arms callout in the corner, and the comparison is the *content*: arm coloring "explains" the layout while superpop coloring does not. Showing them as separate slides loses 30+ seconds of redundant setup.
- **Callbacks.** Slide 07 (all-vs-all heatmap with arm-level dendrogram) sets up the same Jaccard distance matrix at the arm × arm level; this slide unfolds the underlying *sequence × sequence* matrix into 2D. Pre-cue verbally: "same matrix, zoomed out from 41 × 41 arm averages to 15,668 individual flanks."
- **Forward setup for slide 09 (community coloring) and slide 10 (DUX4 / acrocentric / PAR mechanism slides).** This slide deliberately uses the **chromosome** and **superpopulation** colorings, not the **community** coloring, so that slide 09 can introduce the Leiden communities as the *labeled* version of what the audience just saw cluster naturally. Do NOT pre-empt the community coloring here.
- **Hudson Fst footnote.** The numbers (mean 0.044, AFR vs non-AFR 0.10–0.15, non-AFR pairs −0.05 to 0.02, AMR-EUR closest at Fst 0.22 in the out-of-Africa tree) are CROSSWALK §1 ch.04 + ch.12 / end-to-end-report 04_heterogeneity.md lines 105–113. They support the abstract's "resolve subtelomeric clustering across human populations" claim *independent* of whether the projection is PCA or MDS. If the synthesizer rewrites the abstract per (a) above, the C6 Results paragraph should still cite these Fst numbers — they are the substance of the population-clustering claim.
- **Missing introvert arms (2p, 3p, 5p, 8q, 11q, 14q).** Six arms with no detected inter-chromosomal PHR. These are absent from the projection by construction (no Jaccard rows). Mention as a one-line corner annotation; do NOT spend speaker time explaining unless asked in Q&A — the explanation belongs in slide 11/12 (community detection caveats).
- **Figure provenance.** Both panel PDFs already exist in Andrea's PHR_III/similarity/ directory; no new alignment, no SBATCH, no graph rebuild needed. The optional R stitching script is light I/O only and runs in seconds inside the agent worktree.
- **Numbers to lock down (single source of truth).** 15,668 flanks / 41 of 48 arms / 6 missing introvert arms / Dim1 16.05% / Dim2 11.2% / Hudson Fst mean 0.044, AFR vs non-AFR 0.10–0.15. All consistent with CROSSWALK §3 C6 + end-to-end-report 01_pipeline.md and 04_heterogeneity.md.

---

### slide_09_pca_communities_clades

Source: `slides/v2/slide_09_pca_communities_clades.md`

## Title
All-vs-all PCA — 15 arm-level communities, named for the abstract's clades (the keystone slide)

## Bullets
- Same PCA layout as v1 slide 10: each point is one chromosome arm in PHR-similarity space; clusters = Leiden k=15 communities (C1–C15) on the 41×41 arm-level Jaccard distance matrix (chapter 01 §"Community detection"; arms `hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv`).
- **Every named clade in the abstract is a community on this plot.** PAR1 = **C15** (Xp/Yp + 18q outlier). PAR2 = **C14** (Xq/Yq). Acrocentric short arms = **C7** (13p,14p,15p,21p,22p). The abstract's 22q/21q/19q/1q/13q/17q clade = **C6** (exact match). 10p–18p Linardopoulou pair = **C2**. 4q–10q DUX4 = **C1**.
- Three interpretive zones replace v1's question-marked overlays: **PAR-driven** (lower-left, C14+C15), **concerted-exchange PHR core** (center, C6 + C7 + C3), **DUX4 / D4Z4** (right, C1 — confidence upgraded; CROSSWALK §C5 ties it to chapter 03's 28 DUX4L genes and chapter 04's 43.4 % type-discordance).
- Same PC1 / PC2 axes as v1 (16.05 % / 11.2 %); same five missing introvert arms (2p, 3p, 5p, 8q, 11q, 14q). UPGMA k=14 recovers 12/15 of these communities — independent confirmation that the partition is not a Leiden artefact.
- **One slide that ties methods to abstract.** Every clade word the audience just read in the abstract has a colored cluster here, with an arm list, with a recombinational interpretation. After this slide the rest of the talk is *evidence for* what is grouped on this plot.

## Primary figure
**Reuse:** `slides/20260204_Subtelomics_overview_EG.pdf` page 10 — keep the PCA scatter, the chromosome-arm point labels, the PC1 / PC2 axes, and the three interpretive arrows verbatim. Do **not** re-render the scatter (the data lives outside the worktree at `/moosefs/guarracino/HPRCv2/PHR_III/similarity/…`).

**One new asset (this slide's deliverable):** replace v1's right-margin community list with an abstract-anchored legend table. Renders to PDF / PNG locally with vanilla ggplot2 (no `gridExtra`, no `cowplot` — confirmed working in this worktree's R install). Place to the right of (or beneath) the imported v1 PCA.

```r
# slide_09_clade_legend.R
# Renders the abstract-anchored legend that replaces the v1 community list on
# the right margin of the PCA scatter. Runs in seconds, no data dependencies,
# pure ggplot2.
library(ggplot2)

leg <- data.frame(
  row   = 15:1,                       # so C1 lands at top of the rendered table
  C     = sprintf("C%d", 1:15),
  arms  = c("4q, 10q",
            "10p, 18p",
            "3q, 7p, 9q, 11p, 16q, 19p",
            "7q, 12q",
            "6p, 9p, 12p, 20q",
            "1q, 13q, 17q, 19q, 21q, 22q",
            "13p, 14p, 15p, 21p, 22p",
            "15q",
            "16p",
            "17p",
            "1p, 5q, 6q, 8p",
            "2q, 20p",
            "4p",
            "Xq, Yq",
            "18q (n=1), Xp, Yp"),
  clade = c("4q–10q DUX4 / D4Z4",
            "10p–18p (Linardopoulou 2005)",
            "f7501 duplicons (fixed + AFR-enriched)",
            "private 7q/12q pair",
            "RPL23A / WASH duplicons",
            "concerted q-arm clade (22/21/19/1/13/17q)",
            "acrocentric short arms (rDNA-adjacent)",
            "chr15_q (single arm)",
            "chr16_p (single arm)",
            "chr17_p (single arm)",
            "OR4F21 sharing (Linardopoulou block 5)",
            "2q/20p pair",
            "chr4_p (single arm)",
            "PAR2 (Xq/Yq)",
            "PAR1 (Xp/Yp)"),
  fill  = c("#FDE2C8","#FFF1AA","#FFFFFF","#FFFFFF","#FFFFFF",   # C1 DUX4, C2 10p18p
            "#D6E8FF","#E5D8FA","#FFFFFF","#FFFFFF","#FFFFFF",   # C6 q-arm, C7 acro
            "#FFFFFF","#FFFFFF","#FFFFFF","#CDEAD3","#CDEAD3"),  # C14 PAR2, C15 PAR1
  stringsAsFactors = FALSE
)

cols <- data.frame(
  x      = c(0.05, 0.30, 0.70),
  width  = c(0.22, 0.38, 1.05),
  field  = c("C","arms","clade"),
  stringsAsFactors = FALSE
)

cells <- do.call(rbind, lapply(seq_len(nrow(cols)), function(j) data.frame(
  x        = cols$x[j] + cols$width[j]/2,
  width    = cols$width[j],
  y        = leg$row,
  text     = leg[[ cols$field[j] ]],
  fill     = leg$fill,
  fontface = ifelse(j == 1, "bold", "plain"),
  stringsAsFactors = FALSE
)))

hdr <- data.frame(
  x        = cols$x + cols$width/2,
  width    = cols$width,
  y        = 16,
  text     = c("Community","Arms","Abstract clade / interpretation"),
  fill     = "#EEEEEE",
  fontface = "bold",
  stringsAsFactors = FALSE
)

p <- ggplot() +
  geom_tile(data = rbind(cells, hdr),
            aes(x = x, y = y, width = width, height = 0.95, fill = I(fill)),
            colour = "grey75") +
  geom_text(data = cells, aes(x = x - width/2 + 0.01, y = y, label = text,
                              fontface = fontface),
            hjust = 0, size = 3.0) +
  geom_text(data = hdr, aes(x = x - width/2 + 0.01, y = y, label = text,
                            fontface = fontface),
            hjust = 0, size = 3.3) +
  coord_cartesian(xlim = c(0, 1.7), ylim = c(0.5, 16.5), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(6, 6, 6, 6))

ggsave("slide_09_clade_legend.pdf", p, width = 8.0, height = 5.5)
ggsave("slide_09_clade_legend.png", p, width = 8.0, height = 5.5, dpi = 200)
```

The colored rows are exactly the abstract-anchored six clades (C1 / C2 / C6 / C7 / C14 / C15); everything else stays neutral so the eye lands on the rows that match the abstract's wording. The synthesizer can dock this PNG to the right of the v1 PCA scatter (the shape — 8.0 × 5.5 inches — matches the v1 right margin).

## Speaker notes
This is the keystone of the talk — every clade word in the abstract is a colored cluster on this plot. The PCA is on the 41×41 arm-level Jaccard distance matrix; clusters are Leiden communities, k=15 chosen by silhouette. UPGMA on the same matrix recovers twelve of these fifteen — the partition is not an algorithm artefact.

Walk the audience around the plot in the order of the abstract. The pseudoautosomals — PAR2 is **C14**, Xq with Yq; PAR1 is **C15**, Xp with Yp. Acrocentric short arms — **C7** — thirteen-p, fourteen-p, fifteen-p, twenty-one-p, twenty-two-p, all five of them, the rDNA-adjacent homogenization clade. The novel q-arm clade the abstract calls "twenty-two-q, twenty-one-q, nineteen-q, one-q, thirteen-q, seventeen-q" is **C6** — an exact six-arm match, the largest non-acrocentric inter-chromosomal community we see. Ten-p with eighteen-p — Linardopoulou's 2005 pair — is **C2**. And four-q with ten-q, the DUX4/D4Z4 pair, is **C1** — what v1 marked with a question mark is now an established clade in the abstract, twenty-eight DUX4L genes, copy-number diversity, type-discordance forty-three percent.

So this slide says: the abstract's vocabulary maps onto a single empirical structure. Everything that follows is evidence *about* clusters you can already see here.

## Time budget
80 seconds. Roughly: 10 s on "every point is an arm, clusters are communities, k=15 by silhouette, UPGMA agrees on 12/15"; 50 s walking the six abstract-anchored communities (PAR2, PAR1, acrocentric, q-arm clade C6, 10p–18p, 4q–10q DUX4) — about 8 s each, name the clade, point at the cluster, name the arms; 20 s closing line "everything that follows is evidence about clusters you can already see here" + segue to the next slide.

## Notes for synthesizer
- **Layout is preserved from v1 slide 10.** Do not re-render the PCA scatter — keep the v1 panel as-is (PC axes, point shapes, point labels, the three interpretive arrows). The only on-slide change is the right-margin legend, which the R script above generates.
- **Update the three interpretive arrows from v1**: "PARs-driven" → "PAR1 / PAR2 (C15 / C14)"; "PHRs-driven" → "concerted-exchange core (C6 + C7 + C3)"; "DUX4-driven?" → "4q–10q DUX4 (C1)" — *drop the question mark*; the abstract treats it as established.
- **Cross-slide callbacks:**
  - **Slide 07 (v1 all-vs-all heatmap, v2 slide 07/08)** sets up the same arm-level distance matrix; this slide is the projection of that matrix into 2-D.
  - **Slide 08 (PCA by superpop, v1 slide 9)** sets up the projection coordinates without color-coding by community; this slide reuses those coordinates with the community color-coding.
  - **Outbound to slide 10 (within-community heterogeneity / Fig 2a)** — the speaker should land "C7 is the only community where paralog distance is less than allelic distance — the acrocentric p-arms are *fully* homogenized; we will see that on the next slide." That is the cleanest hand-off.
  - **Outbound to slide 11 (population structure / Fig 2c-d)** — communities here become the units for cross-arm Fst and the out-of-Africa tree.
- **Abstract-clade ↔ community map (the load-bearing block — copy verbatim if compressing):**
  - C1 (4q, 10q) ↔ "4q–10q DUX4 with copy-number diversity"
  - C2 (10p, 18p) ↔ "10p–18p" (Linardopoulou 2005 Fig. 5)
  - C6 (1q, 13q, 17q, 19q, 21q, 22q) ↔ "tightly linked clade involving 22q, 21q, 19q, 1q, 13q, and 17q" (exact six-arm match)
  - C7 (13p, 14p, 15p, 21p, 22p) ↔ "acrocentric short arms"
  - C14 (Xq, Yq) ↔ "Xq/Yq via the pseudoautosomal regions" → PAR2
  - C15 (Xp, Yp; 18q n=1) ↔ "Xp/Yp via the pseudoautosomal regions" → PAR1
  - Source: `paper_prep/synthesis/CROSSWALK.md` §C5; `end-to-end-report/report/01_pipeline.md` §"Arm-level community detection".
- **MDS-vs-PCA framing caveat (CROSSWALK §C6 / §3 row):** Andrea's report uses MDS / PCoA on the Jaccard distance, not strict PCA on a feature matrix. The v1 slide called it "PCA" and the abstract uses "principal component … analyses". For the BoG talk this is fine — MDS on a Jaccard distance is metric-equivalent to PCoA, and the audience reads "PCA" as "the 2-D projection of a similarity structure". Speaker should use "PCA" verbally; the methods writer (REWRITE_PLAN TASK-16) is responsible for resolving this in the manuscript.
- **Missing introvert arms (2p, 3p, 5p, 8q, 11q, 14q)** — keep the upper-right callout from v1; do NOT reframe. Their absence is not a quality issue, it is a biological signal (these arms have ~no inter-chromosomal hits, so they have no row in the 41×41 matrix). One sentence in the speaker notes if a Q&A pushes on it.
- **Do not** re-color the points by superpopulation here — that is the previous slide's job. Mixing collapses both slides into one and erases the keystone framing.
- **If the time budget is cut** (e.g., 60 s total for the deck): drop C3 and C5 from the spoken walk (they are duplicon-sharing groups, not in the abstract clade list); never drop PAR2 / PAR1 / C7 / C6 / C2 / C1 — those are the abstract-anchored six.
- **R script is light** (renders a 7.5×5 inch table grob in seconds; no data files, no SBATCH). If the script fails for missing `gridExtra`, fall back to a hand-typed table block in Keynote — content is what matters, rendering pathway is interchangeable.

---

### slide_10_hic_bulk_mantel_exclusions

Source: `slides/v2/slide_10_hic_bulk_mantel_exclusions.md`

## Title

**Hi-C / Pore-C confirm sequence communities are 3D — and the signal *strengthens* when known confounds are removed**

Subtitle: arm-level Mantel of subtelomeric Jaccard similarity vs inter-chromosomal contact, before and after acrocentric / sex / strong-community exclusions.

## Bullets

- **Lead visual: HG002 Pore-C inter-arm contact matrix, 50 kb, 77 arm-haplotypes ordered by sequence community.** Diagonal blocks light up; **B/W = 0.056, p = 3.9 × 10⁻⁸⁵** (within-community contacts vastly outweigh between).
- **Bulk Mantel (similarity × contact) is positive in 7/8 datasets** at 50 kb — CHM13 ρ = 0.66, HG002 Hi-C ρ = 0.66, HG002 Pore-C ρ = 0.49, NA19036 ρ = 0.27, HG02148 ρ = 0.15 (the only borderline sample).
- **Exclude acrocentric p+q + chrX/Y + the four strong communities (D4Z4, acro p, PAR1, PAR2) and every sample's ρ goes up:** HG002 0.66 → **0.79**, CHM13 0.66 → **0.85**, **HG02148 0.15 → 0.72** (the marginal sample becomes one of the strongest).
- **At 10 kb, community-free per-sequence-pair correlations reach** ρ = 0.83 in NA19036 and ρ = 0.81 in HG02148 (`05_hic_validation.md` §"Individual sequence-pair similarity vs Hi-C contact") — finer resolution + no aggregation **strengthens** the signal further.
- **One-line takeaway: sequence similarity predicts 3D contact, and the signal strengthens after confound exclusions.** It is therefore not driven by nucleolar acrocentric clustering, by PAR contact, or by the few largest communities — it is a generic property of subtelomeric homology.

## Primary figure

Two-panel layout, both panels already exist as published assets — no rebuild needed:

- **Left (bulk):** `paper_prep/figures/fig3/figure_fig3.pdf`, **panel (a)** — HG002 Pore-C contact matrix, 50 kb, 77 arm-haplotypes, ordered by Leiden sequence community. Annotation: keep the **B/W = 0.056, p = 3.9 × 10⁻⁸⁵** label visible in the corner.
- **Right (exclusions):** `paper_prep/figures/ed5/figure_ed5.pdf`, **panel (b)** — Mantel ρ before vs after acrocentric + sex exclusion (50 kb), one point per HPRC sample, identity diagonal drawn. Y > X for 7/7. Source TSVs are in `community_based/50000bp/` and `no_acrocentric/50000bp/` (`<sample>_global_test.tsv`).

If the synthesizer wants a single composite for the slide instead of two pasted panels, the script below does it from the published PNGs (no SBATCH, no re-running the analysis):

```r
# slide_10_bulk_mantel_composite.R — bulk Hi-C panel + Mantel exclusion panel side-by-side
# Output: slides/v2/slide_10_bulk_mantel_composite.pdf
suppressPackageStartupMessages({
  library(ggplot2); library(png); library(grid); library(cowplot)
})

bulk      <- readPNG("paper_prep/figures/fig3/figure_fig3.png")    # crop to panel (a) before use
exclusion <- readPNG("paper_prep/figures/ed5/figure_ed5.png")      # crop to panel (b) before use

p_left  <- ggdraw() + draw_image(bulk)      + draw_label(
  "HG002 Pore-C contacts ordered by sequence community\nB/W = 0.056, p = 3.9e-85",
  x = 0.02, y = 0.97, hjust = 0, vjust = 1, size = 10, fontface = "bold")

p_right <- ggdraw() + draw_image(exclusion) + draw_label(
  "Mantel ρ: full vs (no acrocentric + sex)\n7/7 above identity — signal strengthens",
  x = 0.02, y = 0.97, hjust = 0, vjust = 1, size = 10, fontface = "bold")

annotation <- ggdraw() + draw_label(
  "HG002 0.66 → 0.79   CHM13 0.66 → 0.85   HG02148 0.15 → 0.72",
  x = 0.5, y = 0.5, size = 12, fontface = "bold"
) + theme(plot.background = element_rect(fill = "#FFF7E6", colour = "black", linewidth = 0.5))

top  <- plot_grid(p_left, p_right, nrow = 1, rel_widths = c(1, 1), labels = c("a", "b"))
fig  <- plot_grid(top, annotation, ncol = 1, rel_heights = c(1, 0.10))

ggsave("slides/v2/slide_10_bulk_mantel_composite.pdf", fig, width = 13, height = 6.5)
```

The script depends only on the two existing PNGs in the repo and `ggplot2` + `cowplot` + `png` — same toolchain used by neighboring v2 slides.

## Speaker notes

This is the moment we move from sequence to nucleus. Left panel: the HG002 Pore-C inter-arm contact matrix at 50 kilobases, 77 arm-haplotypes, ordered along both axes by their Leiden sequence community. The diagonal blocks are precisely the communities we built from pangenome graph similarity, and they light up — within-community contacts are eighteen-fold higher than between, p ≈ 10⁻⁸⁵. Sequence-defined communities are physical. Right panel: the bulk Mantel test asks whether arms with more similar subtelomeres also contact each other more often. Across eight HPRC datasets the answer is yes in seven — HG002 Hi-C and CHM13 each at ρ = 0.66, Pore-C at 0.49, with HG02148 the marginal exception at 0.15. Now the robustness check. The skeptic's worry is that the signal is just acrocentric nucleolar clustering or pseudoautosomal contact. Strip out the chr13–22 p- and q-arms, chrX, chrY, and the four strongest communities — D4Z4, acrocentric p, PAR1, PAR2 — and every sample's ρ goes up: HG002 to 0.79, CHM13 to 0.85, and HG02148, the marginal sample, jumps to 0.72. Drop to ten-kilobase resolution and treat individual sequence pairs without any community labels and you reach 0.83 in NA19036, 0.81 in HG02148. The signal is not a nucleolar artifact; it is a generic property of subtelomeric homology that gets cleaner the more carefully you look.

## Time budget

80 seconds.

## Notes for synthesizer

- **NEW slide vs v1.** v1 covered Hi-C only briefly. This slide consolidates two figures from the manuscript (Fig 3a + ED5b) into one slot so the talk lands the bulk-Mantel + exclusion-robustness story in 80 s without spawning a separate ED5 slide.
- **Continuity inbound (slide 09):** the previous slide should have set up Hi-C / Pore-C / CiFi as the orthogonal validation modality and introduced the per-haplotype 3D pipeline (`analyze_hic_communities.py`). Open this slide with "the matrix" — the contact matrix in panel a is the picture, not just a number.
- **Continuity outbound (slide 11):** subsequent slides (per-community detail, RPE-1 cell-type/cell-cycle modulation, Dip-C / sperm radial, mouse meiotic, etc.) all build on the *bulk* result that this slide establishes. End on the takeaway sentence so the next slide can refine into per-community / per-cell-type behavior without re-justifying the bulk effect.
- **Source tables to cite verbally if asked:**
  - HG002 Pore-C bulk B/W: `community_based/50000bp/hg002_porec_global_test.tsv` → 0.056, p = 3.9e-85.
  - Mantel full vs exclusions: `05_hic_validation.md` lines 327–351 (full table), with the **no acro pq + sex** column being the headline ("HG002 0.79, CHM13 0.85, HG02148 0.72"). The "no strong" column gives the same direction (HG002 0.765, CHM13 0.837).
  - 10 kb community-free peaks: `05_hic_validation.md` lines 157–166 — NA19036 ρ = 0.827 at 10 kb; HG02148 ρ = 0.809 at 10 kb.
- **Do not invert the convention.** B/W < 1 means *within*-community contacts are higher (between/within ratio). Mantel ρ uses similarity × contact, so positive ρ means more-similar arms contact more. Both go in the "communities are physical" direction; do not flip the signs verbally.
- **CHM13 caveat (only mention if asked):** CHM13 is haploid → ARI is high (0.54) but per-community W/B power is limited (singletons). Bulk Mantel is fine; do not lean on CHM13 per-community numbers.
- **HG002 CiFi was not in the no-acrocentric run** (`05_hic_validation.md` no-acro table covers 7 samples, not 8). This is a known gap in ED5b; do not claim "8/8" for the exclusion plot — say "7/7 sample × technology cells tested".
- **Visual budget:** if the panels won't fit side-by-side at legible size, prioritize panel (a) at full width and shrink panel (b) to a sub-inset with only the diagonal-crossing arrows + the three labelled ρ numbers (0.66 → 0.79 HG002, 0.66 → 0.85 CHM13, 0.15 → 0.72 HG02148). The arrows are the message; the rest is supporting.

---

### slide_11_single_cell_3d

Source: `slides/v2/slide_11_single_cell_3d.md`

## Title
Single-cell 3D — and it works in haploid sperm

## Bullets
- **GM12878 Dip-C, 16 cells (Tan 2018, remapped to T2T-CHM13v2.0):** community arms 6.9% closer in 3D within community than between (W/B = 0.931, Wilcoxon p = 3.8 × 10⁻⁴, Mantel ρ = 0.296).
- **Sperm scHi-C, 20 cells (Xu et al. 2025):** **60% closer** within community (W/B = 0.401, Fisher p = 3.9 × 10⁻⁵¹) — in the haploid, hyper-condensed sperm nucleus.
- **Negative control — pseudo-community "S" of 7 zero-sharing arms** (chr2_p, chr3_p, chr5_p, chr8_q, chr11_q, chr14_q, chr18_q): 11% *farther* in GM12878, 40% *farther* in sperm. Sequence sharing is necessary; community label alone is not enough.
- Same pattern across diploid soma and haploid germline → 3D clustering of subtelomeres is **not** a Hi-C population artefact and **not** restricted to interphase chromatin.
- Sperm is the bridge: it puts the signal in the germline cell that actually transmits the recombination — setting up mouse meiotic Hi-C next.

## Primary figure
**Recommended:** `paper_prep/figures/fig3/figure_fig3.pdf`, **panel (c)** — per-cell C-community W/B vs S_all (negative control), already plotted for both GM12878 (16/16 cells with W/B < 1 inside C; 0/16 inside S_all) and sperm (20/20 vs 1/20). One panel, both datasets, both directions of effect.

**Alternative if a sperm-only zoom is preferred:** `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_mantel_scatter.pdf` paired with `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_16cells_phr_mantel_scatter.pdf` side-by-side. The two-panel composition below stitches them with matching axes and a shared S_all callout — pure ggplot2/cowplot, no SBATCH:

```r
# slide_11_dipc_sperm_pair.R — side-by-side Mantel panels for GM12878 + sperm
library(ggplot2); library(magick); library(grid); library(cowplot)

gm  <- ggdraw() + draw_image(magick::image_read_pdf(
  "/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_16cells_phr_mantel_scatter.pdf",
  density = 200))
spm <- ggdraw() + draw_image(magick::image_read_pdf(
  "/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/sperm_all20_mantel_scatter.pdf",
  density = 200))

label <- function(txt) ggdraw() + draw_label(txt, fontface = "bold", size = 12, hjust = 0.5)

panel <- plot_grid(
  plot_grid(label("GM12878 Dip-C (n = 16)\nW/B = 0.931 · ρ = 0.296 · p = 3.8e-04"), gm,
            ncol = 1, rel_heights = c(0.12, 1)),
  plot_grid(label("Sperm scHi-C (n = 20, Xu 2025)\nW/B = 0.401 · 60% closer · p = 3.9e-51"), spm,
            ncol = 1, rel_heights = c(0.12, 1)),
  ncol = 2
)

neg <- ggdraw() + draw_label(
  "Negative control: 7 zero-sharing arms (S_all)\nGM12878: 11% FARTHER  ·  Sperm: 40% FARTHER",
  size = 11, fontface = "italic"
) + theme(plot.background = element_rect(fill = "#FFF7E6", colour = "black", linewidth = 0.5))

plot_grid(panel, neg, ncol = 1, rel_heights = c(1, 0.12))
ggsave("slides/v2/slide_11_dipc_sperm_pair.pdf", width = 12, height = 6)
```

## Speaker notes
The bulk Hi-C signal you just saw could in principle be a population average artefact — many cells with weak or even absent contacts averaging into a coherent block. Single-cell 3D rules that out. Tan and colleagues' Dip-C reconstructs explicit 3D coordinates per allele in individual GM12878 cells; remapped to T2T-CHM13 with MAPQ-zero retention so we don't lose subtelomeric reads, 16 cells show community-member arms about 7% closer to each other than to non-members — Wilcoxon p of 3.8 × 10⁻⁴, Mantel ρ of 0.30. That's a small effect per cell, but it's there in essentially every cell. Now the punch line: we ran the same pipeline on twenty sperm cells from Xu and colleagues. Sperm is haploid, the chromatin is hyper-condensed, and the nuclear architecture is nothing like an interphase lymphoblast. The within-versus-between ratio is 0.40 — community arms are sixty percent closer in 3D — Fisher p of 3.9 × 10⁻⁵¹. Same direction, much stronger. The negative control closes the loop: a pseudo-community of the seven arms that share *no* subtelomeric sequence with anything else moves the *opposite* way — 11% farther in GM12878, 40% farther in sperm. Sequence sharing is necessary for clustering; a community label alone is not. So the 3D signal survives the bulk-to-single-cell test and survives the haploid-germline test. That last point is the bridge — sperm is the gamete that actually carries the recombination forward, and the next slide takes us into the meiotic cells where that recombination happens.

## Time budget
60 seconds.

## Notes for synthesizer
- **NEW slide vs v1:** v1 has no single-cell 3D content; this is a fresh contribution that lands two abstract claims at once — that the 3D signal is per-cell and that it generalizes to germline architecture. Both are needed for the abstract's "we hypothesize that these patterns are maintained by recombination facilitated by the physical proximity of subtelomeres" line.
- **Continuity inbound (slide 10):** previous slide is the bulk Hi-C / Pore-C 3D signal (community-level B/W < 1, Mantel rho on bulk maps). Open this slide with a single transition sentence: "Bulk Hi-C is a population average — does the signal survive at single-cell resolution?" Then dump straight into GM12878 → sperm.
- **Continuity outbound (slide 12):** next slide is mouse meiotic Hi-C (zygotene peak, lepto→zygo→pachy→diplo trajectory). The natural pivot is sperm → meiosis: "if the haploid product still shows it, what about the cells where the recombination happens?" Land the last bullet on that pivot — do **not** redo the meiosis story here.
- **Negative control is load-bearing.** The S_all pseudo-community result (11% / 40% *farther*) is the single best rebuttal to the "you're just measuring chromosome-territory crowding" objection. If the visual gets cropped, keep S_all in the bullet text at minimum.
- **Citations:** Tan et al. 2018 (Dip-C, GM12878); Xu et al. 2025 (sperm scHi-C). Both are in `paper_prep/synthesis/REFERENCES.bib`. Safe to cite verbally.
- **Provenance for numbers:** `end-to-end-report/report/06_dipc_validation.md` §"Community 3D enrichment (T2T)" and §"3D genome validation: sperm single-cell". Source TSVs live under `/moosefs/guarracino/HPRCv2/dipc_t2t/output_q0_XX/community_enrichment_k50/` (GM12878) and `/moosefs/guarracino/HPRCv2/dipc_t2t/sperm/enrichment_corrected/` (sperm). Per-cell PDFs are listed in **Primary figure**.
- **Scope reminder:** AUDIT_REPORT / CROSSWALK currently slot Dip-C and sperm into "SI" or "do not include" for the manuscript itself. For the *talk* they are central — they are what convinces the room that the 3D claim is robust. Do not let manuscript scope decisions shrink the talk version.
- **Don't oversell Mantel ρ on sperm.** Sperm Mantel ρ = 0.202 (p = 0.023, significant but modest); the headline number is W/B = 0.401. Keep emphasis on W/B + the 60% framing.
- **Time discipline:** 60 s is tight. Single figure, three numbers (GM12878 6.9%, sperm 60%, S_all neg-control), one bridge sentence to slide 12. If forced to cut, drop the Wilcoxon-vs-Mantel detail; never drop the negative control.

---

### slide_12_mouse_meiotic_zygotene_bouquet

Source: `slides/v2/slide_12_mouse_meiotic_zygotene_bouquet.md`

## Title

**Mouse meiotic Hi-C — the zygotene bouquet is where the 3D signal peaks**

Subtitle: bulk human Hi-C is somatic; mouse zygotene Hi-C (Zuo 2021) is the only meiotic 3D map available — and it is precisely where similarity-vs-contact lights up.

## Bullets

- **Bulk human Hi-C is mitotic** (LCLs, RPE-1, blood). The recombination we are explaining is **meiotic**. Mouse zygotene Hi-C from Zuo et al. 2021 is the only meiotic 3D dataset on a T2T-grade genome — *Patel-style 4-stage Hi-C: leptotene → zygotene → pachytene → diplotene*.
- **Mantel ρ peaks at zygotene** (50 kb, 1 Mb subtelomeric windows): **leptotene 0.687 / zygotene 0.718 / pachytene 0.683 / diplotene 0.577**. The per-PHR-pair Spearman matches: **ρ = 0.715, p = 4.4 × 10⁻⁵⁵, n = 344 inter-chromosomal pairs**.
- **Zygotene is the bouquet stage** — telomeres cluster at the LINC-anchored nuclear envelope to align homologs (Mefford 2002 / Linardopoulou 2005 framing in chapter 07). The 3D signal is strongest exactly when telomeres are physically clustered, then decays as the bouquet resolves through pachytene → diplotene.
- **Mouse T2T B6 + CAST (Francis 2025), 39 p-arm flanks have signal; 39 q-arm flanks have zero.** Mouse chromosomes are telocentric — q is centromere-proximal — confirming the signal lives in the subtelomeric end exactly as in human.
- **One-line takeaway: the human-LCL Hi-C signal is the somatic shadow of a meiotic phenomenon, and we can see the meiotic version directly in mouse — at the bouquet stage.**

## Primary figure

**Recommended:** `paper_prep/figures/fig4/figure_fig4.pdf`, **panel (d)** — mouse zygotene per-PHR-pair scatter (B6 + CAST T2T, Zuo 2021, 50 kb): Spearman ρ = 0.715, p = 4.4 × 10⁻⁵⁵, n = 344 inter-chromosomal pairs. Already publication-ready; no rebuild needed.

**Standalone alternative (zygotene scatter only, larger):** `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp/zuo2021_zygotene_phr_pair_scatter.pdf`.

**New companion panel — stage trajectory (the bouquet peak).** Panel 4d shows zygotene only; the *peak* claim needs all four stages on screen. Light ggplot2 inset, no SBATCH, runs in the agent worktree:

```r
# slide_12_meiotic_stage_trajectory.R — Mantel ρ across the 4 meiotic stages
# Output: slides/v2/slide_12_stage_trajectory.pdf  (~4 in × 2.5 in, intended as inset
# next to panel 4d so the audience sees zygotene as a peak, not a single number)
suppressPackageStartupMessages({ library(ggplot2) })

# Source: end-to-end-report/report/08_mouse.md, "Mouse meiotic Hi-C validation (1Mb window)"
# section, Mantel rho column at 50 kb resolution, 4 stages, B6+CAST per-haplotype.
df <- data.frame(
  stage = factor(c("leptotene","zygotene","pachytene","diplotene"),
                 levels = c("leptotene","zygotene","pachytene","diplotene")),
  rho   = c(0.687, 0.718, 0.683, 0.577),
  is_bouquet = c(FALSE, TRUE, FALSE, FALSE)
)

ggplot(df, aes(x = stage, y = rho, group = 1)) +
  geom_line(linewidth = 0.7, colour = "#444444") +
  geom_point(aes(colour = is_bouquet, size = is_bouquet)) +
  geom_text(aes(label = sprintf("%.3f", rho)), vjust = -1.0, size = 3.5) +
  scale_colour_manual(values = c(`FALSE` = "#1f77b4", `TRUE` = "#d62728"),
                      guide = "none") +
  scale_size_manual(values = c(`FALSE` = 2.4, `TRUE` = 4.0), guide = "none") +
  scale_y_continuous(limits = c(0.50, 0.78), breaks = seq(0.50, 0.75, 0.05)) +
  annotate("text", x = "zygotene", y = 0.76,
           label = "bouquet\n(telomeres clustered\nat nuclear envelope)",
           size = 3.0, hjust = 0.5, colour = "#d62728", fontface = "bold") +
  labs(x = "meiotic prophase stage (Zuo et al. 2021)",
       y = "Mantel ρ (similarity × Hi-C contact, 50 kb, 1 Mb window)",
       title = "Mouse meiotic Hi-C: zygotene peak") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"))

ggsave("slides/v2/slide_12_stage_trajectory.pdf", width = 4.0, height = 2.5)
ggsave("slides/v2/slide_12_stage_trajectory.png", width = 4.0, height = 2.5, dpi = 300)
```

Numbers in `df` are taken verbatim from `end-to-end-report/report/08_mouse.md` lines 100–103 (Mantel ρ column, 50 kb, 1 Mb window — the canonical row used in figure 4d's caption summary).

## Speaker notes

A skeptic could fairly ask: bulk Hi-C is mitotic — LCLs, RPE-1, PBMCs — but the recombination you are explaining happens in meiosis. Why should somatic 3D contact predict meiotic exchange? Mouse meiotic Hi-C is the answer. Zuo and colleagues in 2021 sorted four prophase stages — leptotene, zygotene, pachytene, diplotene — and produced stage-specific Hi-C on a genome we now have a T2T assembly for. Run the same Mantel test we used in human against B6 plus CAST mouse subtelomeres: ρ across the four stages is 0.687, **0.718**, 0.683, 0.577. Zygotene is the peak. Per-PHR-pair Spearman at zygotene is ρ = 0.715, p = 4.4 × 10⁻⁵⁵ across 344 inter-chromosomal pairs. Zygotene is the bouquet — telomeres are clustered at the LINC-anchored nuclear envelope while homologs align. The 3D signal we see in human LCLs is the somatic shadow of a meiotic event we can watch directly in mouse, at exactly the stage when telomeres are physically together. Cross-species generality is a bonus.

## Time budget

**60 seconds.** ~15 s on the somatic-vs-meiotic skeptic frame, ~25 s on the four-stage trajectory with the zygotene peak (point at the inset, say all four numbers — they tell the story), ~15 s on "this is the bouquet" with the LINC / Mefford / Linardopoulou anchor, ~5 s segue into whatever the synthesizer puts next (pedigree exchanges, the empirical proof that the contact does drive recombination).

## Notes for synthesizer

- **NEW slide vs v1.** v1 deck has no mouse / no meiotic content. This slide is genuinely new and lands the meiotic-3D bridge that the abstract's last sentence ("Hi-C-derived three-dimensional genome maps") implicitly relies on but does not explicitly defend. CROSSWALK §08_mouse flags this as **out-of-scope for the canonical Nature companion** but **in-scope for the talk** — the zygotene peak is the strongest existing evidence for the meiotic interpretation in C8. Do not promote it to a manuscript main figure.
- **Continuity inbound (slide 10/11):** slide 10 establishes the *bulk human* Hi-C / Pore-C signal (ρ ≈ 0.66, B/W ≈ 0.056, signal *strengthens* after acrocentric+sex exclusions) — i.e., the somatic 3D-mirrors-sequence result. Open this slide with the exact pivot: "the human Hi-C is somatic — the recombination is meiotic — what about a meiotic 3D map?" If slide 11 covers single-cell / Dip-C / sperm intermediates, this slide is the *meiotic capstone* of the 3D arc.
- **Continuity outbound (slide 13+):** the natural follow-on is the pedigree recombination evidence (Fig 4a — WashU 3-gen T2T, 92% within-Leiden patches). Land the segue: "we see the contact at the bouquet stage; here is the recombination it produces, one generation later." If the synthesizer keeps a single mouse slide and a single pedigree slide adjacent, the talk's central causal claim — sequence similarity → 3D contact at bouquet → exchange — is delivered in three slides (10 → 12 → 13).
- **The trajectory inset is load-bearing.** Panel 4d alone shows zygotene as a single number. The *peak* claim is `0.687 / 0.718 / 0.683 / 0.577` — leptotene rises into zygotene, then decays. Without the inset, the audience hears "zygotene = 0.715" and has to take the peak on faith. Keep the inset, even at small size, and orient panel 4d + inset side-by-side so the eye moves "lots of stages → which one peaks → here's the scatter for the peak."
- **Numbers to lock down:** Mantel ρ 0.687 / 0.718 / 0.683 / 0.577 (lepto/zygo/pachy/diplo, 50 kb, 1 Mb window). Per-PHR-pair Spearman ρ = 0.715, p = 4.4 × 10⁻⁵⁵, n = 344 inter-chromosomal pairs at zygotene. Source of truth: `end-to-end-report/report/08_mouse.md` §"Mouse meiotic Hi-C validation (1Mb window)" lines 87–113 (community-free) and 130–147 (community-based / multi-window). All numbers also appear in `paper_prep/synthesis/CROSSWALK.md` row 41.
- **Bouquet anchor.** "Bouquet" / "telomeres at the nuclear envelope" / "LINC complex" framing comes from Mefford 2002, Linardopoulou 2005, and Zuo 2021 — already synthesized in chapter 07 (Discussion C8 source). Do not re-derive the literature here; one phrase ("telomeres cluster at the nuclear envelope") is enough.
- **Telocentric caveat.** "39 p-arm signal / 39 q-arm zero" is genuinely informative — it tells the audience the pipeline is finding the right end. Do NOT spin it as "missing data"; mouse q-arms are centromere-proximal and biologically uninformative for subtelomeric homology. One sentence is enough.
- **Window-size / B/W robustness** (1 Mb → 2 Mb → 4 Mb, ρ 0.58–0.73) is in the substrate but not on the slide. If a reviewer in Q&A asks "is this just a 50 kb / 1 Mb artifact?", the answer is "Mantel ρ holds at 0.65–0.73 across 5 resolutions × 3 window sizes × 4 stages — see chapter 08 §window-size optimization." Do not put it on screen; a single line on the deck is too dense.
- **Do not call this "validation in mouse."** The framing in the pitch is *meiotic bridge*, not *cross-species validation*. The cross-species generality is a bonus, mentioned in one clause; the load-bearing claim is *meiotic 3D map shows the same effect at the bouquet*. CROSSWALK is explicit: mouse is not a manuscript main figure, but the bouquet-stage peak is the cleanest evidence for meiotic causation.
- **No SBATCH needed** for the inset; the R script reads four numbers hard-coded from the report and produces a 4 in × 2.5 in PDF in seconds inside the agent worktree.

---

### slide_13_pedigree_direct_evidence

Source: `slides/v2/slide_13_pedigree_direct_evidence.md`

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

---

### slide_14_gene_biology

Source: `slides/v2/slide_14_gene_biology.md`

## Title

Gene biology aside — DUX4, OR4F, TAR1 (the biology is interesting too)

## Bullets

- **DUX4 (FSHD locus).** Annotated on 18 q-arms across the pangenome, but only chr4q and chr10q (community **C1**) carry the full **D4Z4 macrosatellite** array (median **22** DUX4L copies; the FSHD-permissive 4qA haplotype lives here). All other arms: 0–2 copies — DUX4 has scattered, but the disease-relevant repeat unit has not.
- **OR4F (olfactory receptors).** 4 OR4F paralogs span 16 arms, **5,023** gene-copy entries total. Pseudogenisation runs as a clean per-arm gradient: **11.1% pseudogene at chr7p → 99.8% at chr15q** (population mean 62.1%). Same gene, same neighborhood — but the decay clock has run longer at one end.
- **TAR1 (telomere-associated repeat).** 21,544 entries across **94.6%** of all 15,668 PHR sequences and all 41 arms — universal except at PAR1 (chrX_p / chrY_p, **0.5%**). PAR1 has obligate meiotic crossover; satellite-mediated exchange anchors are evidently not required there.
- One sentence: **distinct biological histories — disease, decay, and exchange machinery — write themselves into the same subtelomeric architecture.**

## Primary figure

`slides/v2/slide_14_gene_biology.pdf` (this worktree) — produced by `slides/v2/slide_14_gene_biology.R`. Three panels, 16 × 5.2 in:

- **Panel a (DUX4):** boxplot of DUX4L copies per haplotype on chr4q and chr10q (the only arms with full D4Z4 arrays); red dashed line at C1 median = 22. Annotation reminds the audience that DUX4 is annotated on 18 q-arms but only C1 carries the macrosatellite.
- **Panel b (OR4F):** per-arm pseudogene fraction sorted ascending; chr7p (11.1%) and chr15q (99.8%) extremes labelled in red; population mean 62.1% as red dashed reference.
- **Panel c (TAR1):** per-arm TAR1 prevalence (% sequences carrying TAR1) sorted ascending; PAR1 arms (chrXp 0.3%, chrYp 1.1%) highlighted in blue against grey for autosomal arms (≥73%); all-PHR mean 94.6% as grey dashed reference.

Data sources (all already on disk; the R script reads them directly):

- `/moosefs/guarracino/HPRCv2/PHR_III/plots/d4z4_dux4l_by_community.tsv` (1,253 chr4q/chr10q haplotypes)
- `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv` (16 arms, 5,023 OR4F entries)
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_by_arm.tsv` (41 arms)

Cross-reference: `paper_prep/figures/ed4/figure_ed4.pdf` panels c (high-copy gene families) and d (OR4F pseudogenisation gradient) are the manuscript versions of the OR4F + DUX4 content; this slide's figure is the talk version (one row, three panels, larger labels) and adds the TAR1 panel that ED4 omits. Do **not** swap them — keep ED4 for the manuscript, slide_14 figure for the talk.

## Speaker notes

A quick aside before I move to the punchline. The talk has been about geometry — communities, exchange, scale — but the biology inside these regions is genuinely interesting in its own right, and worth flagging. Three vignettes, ten seconds each.

First — DUX4. We see DUX4 annotated across 18 q-arms in the pangenome. But the medically relevant biology — the D4Z4 macrosatellite array whose contraction causes facioscapulohumeral muscular dystrophy — lives only at chr4q and chr10q, the C1 community. Median 22 copies of DUX4L per haplotype on those two arms; everywhere else, just 0–2. So DUX4 has spread, but FSHD is geometrically constrained to one community out of fifteen.

Second — OR4F olfactory receptors. Four paralogs, sixteen arms, five thousand gene-copy entries. The pseudogenisation rate sweeps from eleven percent at chr7p to ninety-nine point eight percent at chr15q. Same gene family, distributed across the same subtelomeric exchange network — but the decay clock has been running for very different lengths of time at different ends.

Third — TAR1, the telomere-associated repeat. Ninety-four point six percent of all our subtelomeric sequences carry TAR1, across all forty-one arms. The one exception is PAR1 — chrX_p and chrY_p — at half a percent. PAR1 has obligate meiotic crossover. So the satellite seems to mark places that *need* a sequence-anchor for exchange; PAR1 doesn't, and TAR1 isn't there.

Three different biological readouts — disease, decay, exchange machinery — and all three write themselves into the architecture we just spent the talk mapping. **This is a digression from the core argument**, and if I'm running short I'll skip it. But if you remember one thing: this is the substrate the rest of human biology cares about — FSHD lives here, the olfactory repertoire ages here, and the exchange machinery itself leaves footprints here.

## Time budget

50 seconds. **Compressible to 0** — explicitly skippable if the previous slides have run long. If kept, run all three vignettes at ~15 s each plus 5 s framing.

## Notes for synthesizer

- **Framing — this is an ASIDE, not the core argument.** The slide opens "and the biology is interesting too" and the speaker note ends with "this is a digression from the core argument; if I'm running short I'll skip it." The synthesizer should preserve this skippable framing. The talk's main argument is geometry / exchange / population-scale (slides 1–13 and 15); slide 14 is a "look how rich this is" detour. Do not promote any of the three vignettes to a load-bearing claim of the talk.
- **Compressible to 0.** This is the cleanest slide to drop if the talk runs long. The 15-min slot is tight (15 slides, ~60 s each on average, with anchor slides — title, central PHR scale, PAR/PHR comparison, conclusions — needing more than 60 s). If timing is tight, slide 14 is the cut. The synthesizer should flag this explicitly to the speaker (Erik Garrison).
- **Numbers to lock down (single source of truth: ABSTRACT.md, CROSSWALK §C12, end-to-end-report 02_annotation.md and 03_gene_enrichment.md):**
  - DUX4: 18 q-arms (annotation), C1 = chr4q + chr10q (the D4Z4 community), median 22 DUX4L per haplotype on C1, 0–2 elsewhere. Andrea ch. 03 cites "DUX4L1–DUX4L44 (28 genes, 2 arms)" — the 28 distinct DUX4L paralog labels are the canonical CHM13 annotation; the 22-median is the per-haplotype copy count from `d4z4_dux4l_by_community.tsv`. These are not in conflict — different denominators.
  - OR4F gradient: 11.1% (chr7p) → 99.8% (chr15q), 16 arms, 5,023 entries, mean 62.1%. Verbatim from CROSSWALK §C12 / SURVEY_10_11_12 line 344.
  - TAR1: 94.6% of 15,668 sequences carry TAR1 across 41 arms; PAR1 = 0.3% (chrXp) and 1.1% (chrYp). 21,544 total TAR1 entries. From end-to-end-report 02_annotation.md ("TAR1 prevalence" subsection).
- **Figure provenance.** The slide figure is **slide-specific** (newly generated by `slide_14_gene_biology.R` in this worktree, using existing cluster TSVs — no SBATCH, no new data). It overlaps panels c+d of `paper_prep/figures/ed4/figure_ed4.pdf` (manuscript ED4) but adds the TAR1 panel and reformats for talk display. The synthesizer should not substitute ED4 verbatim — ED4 has 4 panels including a GO:BP plot and a copy-weighted comparison panel that are off-topic for this slide.
- **R script renders cleanly.** `Rscript slides/v2/slide_14_gene_biology.R` regenerates `.pdf` and `.png` from the cluster TSVs in ~3 s; no external dependencies beyond base R + readr/dplyr (both already in the project's R environment). Do not require patchwork / cowplot / gridExtra — they are NOT installed in the worktree's R environment.
- **Callbacks/setup.** The previous slide(s) (12, 13) cover within-community heterogeneity / cross-arm exchange status; slide 14 picks up "what's actually IN these regions?" framing without claiming any of it as a finding the speaker needs to defend. Slide 15 (closing / conclusions) should NOT cite slide 14 numbers as load-bearing — slide 14 is decorative.
- **Forward setup.** Slide 15 (concluding slide) returns to the geometry argument. Slide 14's three vignettes don't feed into it; slide 14 is a parallel "by the way" detour and the speaker explicitly returns to the main thread on slide 15.
- **PAR1 callback.** The TAR1 panel says "PAR1 is the only place that doesn't need a satellite-mediated exchange anchor." This implicitly connects to the PAR-vs-PHR framing established earlier in the deck (slide 05's "PAR2-class exchange landscape replicated 41 times"). The synthesizer can choose to make this callback more explicit if desired — but the bullet/note already gestures at it.
- **No new R packages required.** Verified the worktree has base R + ggplot2 + dplyr + readr but NOT patchwork / cowplot / gridExtra. Any future re-cut should stay within base R `layout()` / `par(mfrow=...)`.

---

### slide_15_concerted_evolution_thesis

Source: `slides/v2/slide_15_concerted_evolution_thesis.md`

## Title

Concerted evolution of human subtelomeres — what we saw, predicted, and recovered

## Bullets

- **Method (slide 3).** Implicit pangenome graph: wfmash all-vs-all over 18,827 telomere-anchored flanks, ~12% of all pairs evaluated — **230× above the Erdős-Rényi connectivity threshold** (`p* = log(n)/n ≈ 5.21×10⁻⁴`). No chromosome partitioning, no GFA: every haplotype is its own reference.
- **Empirical (slides 4–9).** **15,668 PHRs across 41/48 arms**, median 105 kb / mean 144 kb — PAR2-scale pseudohomology at nearly every chromosome end. Named clades: Xp/Yp & Xq/Yq via PARs, **acrocentric short arms** (C7, near-interchangeable), **10p–18p**, the big q-arm clade (**22q–21q–19q–1q–13q–17q**), and the **4q–10q DUX4** clade.
- **Mechanism (slides 10–12).** Hi-C/Pore-C/CiFi/Dip-C all independently recover community-structured 3D contacts (B/W 0.027–0.074; Mantel ρ=0.296, p=0.002; per-pair ρ=0.674 in CHM13 Hi-C). Median PHR (105 kb) sits at the base of a single meiotic loop — the **bouquet** is the predicted exchange venue.
- **Proof (slide 13).** **WashU T2T pedigree: 494/538 (92%) inter-chromosomal patches fall inside Leiden communities** — 133 gene-conversion-like, 16 crossover-like. The graph-derived community structure *predicts* where exchanges show up in real families.
- **Biology (slide 14).** D4Z4 / DUX4 (4q↔10q) is the disease-revealed instance of the same process — FSHD-modifying translocations are concerted evolution caught in the act.

> **Thesis: subtelomeres concertedly evolve through ongoing inter-chromosomal exchange — observable in pedigrees, predicted by 3D, recovered by an implicit pangenome graph across HPRC v2.**

## Primary figure

**Recommended:** `paper_prep/figures/ed8/figure_ed8.pdf` panel **(a)** — the four-link causal feedback loop (sequence sharing → 3D proximity → ectopic exchange → new shared segments → propagation), with edges color-coded by evidence type (solid blue = direct measurement, solid green = literature, dashed olive = inferred). This is the closer-friendly schematic; Andrea built it as the synthesis figure of the discussion (`end-to-end-report/report/07_integrated.md §1.6` → `paper_prep/figures/ed8/`).

**Why ed8(a) and not a new figure:** the closer should *recapitulate*, not introduce. Ed8(a) already encodes every pillar of the talk — sequence (slide 4–9), 3D (10–12), exchange (13–14) — as a single loop. No new compute, no SBATCH, no R needed.

If the synthesizer wants a thin annotation overlay that names each edge with its slide number (so the audience can map the loop back to what they saw), the following snippet wraps the existing PNG with four labels — runs locally in seconds, no data dependency:

```r
# slide_15_loop_callouts.R — overlay slide-number tags on the ed8(a) feedback loop
# Output: slide_15_loop_with_callouts.pdf  (drop-in for the slide background)
library(ggplot2); library(png); library(grid); library(cowplot)

panel_a <- readPNG("paper_prep/figures/ed8/figure_ed8.png")
bg      <- ggdraw() + draw_image(panel_a)

tag <- function(text, x, y) {
  ggdraw() + draw_label(text, x = 0.5, y = 0.5, hjust = 0.5, size = 9,
                        fontface = "bold", colour = "#1b3a6f") +
    theme(plot.background = element_rect(fill = "#FFF7E6",
                                         colour = "#1b3a6f", linewidth = 0.4))
}

# Four edge tags — coordinates are nominal; tune to the actual ed8(a) layout
ggdraw(bg) +
  draw_plot(tag("sequence sharing\n(slides 4–9)"),    x = 0.04, y = 0.78, width = 0.20, height = 0.07) +
  draw_plot(tag("3D proximity\n(slides 10–12)"),      x = 0.78, y = 0.78, width = 0.20, height = 0.07) +
  draw_plot(tag("ectopic exchange\n(slide 13 pedigree)"), x = 0.78, y = 0.10, width = 0.22, height = 0.07) +
  draw_plot(tag("propagation\n(slide 14 DUX4)"),      x = 0.04, y = 0.10, width = 0.20, height = 0.07)

ggsave("slides/v2/slide_15_loop_with_callouts.pdf", width = 12, height = 7)
```

(Synthesizer: ed8(a) is the recommended *primary* — the snippet above is optional polish, not a blocker. If the synthesizer skips the overlay, the bullets already do the slide-number tagging in words.)

## Speaker notes

This is the close. One slide, one breath, one thesis.

We started with a methodological commitment: the implicit pangenome graph. Wfmash over 18,827 telomere-anchored flanks, twelve percent of all pairs evaluated — and that twelve percent is two hundred and thirty times above the Erdős-Rényi threshold for graph connectivity. That is what licenses everything that followed: no chromosomal partitioning, every haplotype its own reference.

What that method showed us is that pseudohomology at PAR2 scale is replicated at nearly every chromosome end — fifteen thousand six hundred and sixty-eight PHRs across forty-one of forty-eight arms, median one hundred and five kilobases. We named the clades: PARs at Xp/Yp and Xq/Yq, the acrocentric short arms, ten-p with eighteen-p, the big q-arm clade — twenty-two q, twenty-one q, nineteen q, one q, thirteen q, seventeen q — and four q with ten q carrying DUX4.

We then asked: why are these in clades? Hi-C, Pore-C, CiFi, and single-cell Dip-C all independently recover community-structured three-dimensional contacts, with per-pair correlations as high as point-six-seven-four in CHM13. The median PHR fits at the base of a single meiotic loop — the bouquet is the predicted exchange venue.

And then the proof. The WashU T2T pedigree gives us four-hundred-ninety-four out of five-hundred-thirty-eight inter-chromosomal patches — ninety-two percent — falling *inside* the Leiden communities the graph predicted. That is the loop closing: sequence similarity, 3D proximity, and observed family-level exchange all agreeing.

D4Z4 and DUX4 are the disease-revealed instance — FSHD translocations are concerted evolution caught in the act on a clinically named locus.

The thesis: **subtelomeres concertedly evolve through ongoing inter-chromosomal exchange — observable in pedigrees, predicted by 3D, recovered by an implicit pangenome graph across HPRC v2.** Thank you.

## Time budget

**Target: 70 seconds.** Roughly 10 s on the method recap (Erdős-Rényi number stays in), 15 s on the empirical pillar (PHR scale + named clades — name 3–4 of them, don't list all), 15 s on the mechanism (3D convergence + bouquet, one Mantel number), 15 s on the pedigree proof (the 92 % number is the headline, do not skip it), 5 s on DUX4 as the disease-revealed instance, 10 s on the one-line thesis read aloud verbatim. The thesis sentence is the last thing the audience hears — do not paraphrase it on the fly.

## Notes for synthesizer

- **What slide 14 sets up.** Slide 14 lands DUX4/FSHD as the case study in concerted evolution caught in the act on a disease locus. Slide 15 must *not* re-explain DUX4; it cites it as the mechanistic exemplar in one bullet and moves on. The transition from 14 → 15 is "and that's not just one locus — here's the whole picture."
- **Callback discipline.** Every bullet on this slide is tagged to a prior slide number. **Do not drop the tags** in compression. The whole point of a closer is recapitulation; the audience needs the visual cue ("slide 13," "slides 10–12") to remember they already saw the evidence and trust the synthesis.
- **The thesis sentence is locked from the task spec.** "Subtelomeres concertedly evolve through ongoing inter-chromosomal exchange — observable in pedigrees, predicted by 3D, recovered by an implicit pangenome graph across HPRC v2." That exact wording is what the task asked for and what should appear on the slide. Synthesizer: this is a hard constraint, not a suggestion. Render as a blockquote / pull-quote, large type, bottom of slide.
- **Title-callback.** The manuscript title (slide 01) is "Concerted evolution and unorthodox recombination of human subtelomeres" (`paper_prep/synthesis/ABSTRACT.md`). Slide 15's headline echoes the "concerted evolution" half deliberately so the deck closes by closing the title — first slide and last slide bracket the same phrase.
- **Numbers to lock down (single source of truth: CROSSWALK + ABSTRACT).**
  - 18,827 telomere-anchored flanks; ~12% sampling; ER `p* ≈ 5.21×10⁻⁴`; ~230× threshold (slide 3)
  - 15,668 PHRs / 41 of 48 arms / median 105 kb / mean 144 kb (slide 5; matches `framing_synthesis.md`)
  - PAR2 ≈ 334 kb (slide 4 anchor)
  - Hi-C B/W 0.027–0.074, p 6.0e-18 to 9.4e-03; Mantel ρ=0.296, p=0.002; CHM13 per-pair ρ=0.674 (slides 10–12; from `07_integrated.md`)
  - WashU pedigree 494/538 = 92% within-community; 133 gene_conversion_like, 16 crossover_like (slide 13; from `14_pedigree_recombination.md`)
  - C1 = chr4_q/chr10_q, median 22 DUX4L, Mann-Whitney p=5.3e-6 (slide 14; from `07_integrated.md` D4Z4-CTCF section)
  - **All six numeric clusters above must agree with the slides they call back to.** If the synthesizer notices a drift between, e.g., slide 13's pedigree number and the number cited here, the slide-13 author is the source of truth — file a `wg msg` rather than silently editing this slide.
- **Figure choice rationale.** Ed8(a) is preferred because (i) it already exists, (ii) Andrea built it specifically as the discussion synthesis figure (`paper_prep/figures/ed8/caption.md`), (iii) it encodes the talk's argument as a single loop. The optional R overlay in *Primary figure* is a "nice to have" — if the synthesizer is short on time, ship ed8(a) bare.
- **Forward setup.** This is the last slide before Q&A. Do **not** add a "future work" / "thanks" slot inside this slide — the task spec is explicit (15 slides total, this is slide 15). Acknowledgements, if any, belong in a separate appendix slide that is not part of the 15-slide budget.
- **If compression is forced** (deck overruns and the synthesizer needs to cut): the order to drop is (1) the optional R overlay, (2) the named-clades enumeration in bullet 2 (drop 10p–18p and the big q-arm clade list, keep PAR + acro + 4q/10q), (3) the per-technology B/W numbers in bullet 3 (keep Mantel ρ=0.296 as the headline). **Do not drop:** the 230× ER number, the 92 % pedigree number, or the thesis pull-quote. Those three are the slide.
- **No SBATCH needed.** Ed8(a) PDF/PNG already exists at `paper_prep/figures/ed8/figure_ed8.{pdf,png}`. The optional overlay R script runs in seconds in the agent worktree and writes a single PDF.
- **No edits to other v2 slides.** Task spec: single output file, single commit. Cross-slide concerns are reported here for the synthesizer; this slide does not modify any neighbor.

---
