# coherence_check.md — BoG v2 deck

Cross-slide consistency / narrative-gap audit of `slides/v2/slide_01_*.md` … `slide_15_*.md`. Synthesised by `bog-v2-slides` (agent-809) on 2026-05-06.

Format: each entry is a flag (severity), the concrete inconsistency, the slides involved, and a recommended resolution. Severity is **Blocker** (talk loses meaning if not fixed), **Material** (audience would notice; should fix before render), or **Polish** (consistency hygiene; leave for post-talk).

---

## 1. Blocker — v1 deck PDF missing from this worktree

**Slides involved:** 02, 03, 09.

Three slide files reference `slides/20260204_Subtelomics_overview_EG.pdf` as a primary or recommended figure source (slide 02 page 2, slide 03 page 3, slide 09 page 10). The PDF is **not present** in this worktree (`Glob slides/20260204*` returns no matches). Confirmed 2026-05-06 by `Bash ls slides/`.

Each affected slide is downstream of the v1 deck for one specific figure:
- Slide 02: implicit-interval-tree diagram (page 2) — load-bearing visual; hard to render fresh without the schematic asset
- Slide 03: wave → dotplot → forest pipeline visual (page 3) — load-bearing for the methods anchor
- Slide 09: arm-level PCA / community scatter (page 10) — has on-disk alternate at `/moosefs/guarracino/HPRCv2/PHR_III/similarity/hprcv2.1Mb.subtelo.mds.arm-leiden-k15.communities.pdf` (verified to exist), so this one has the cleanest fallback

**Recommended resolution:** restore the v1 deck PDF to `slides/` from a backup or peer worktree before render. If that fails, slide 02 needs a freshly drawn interval-tree schematic (likely Keynote / TikZ); slide 03 can extract the wave-dotplot-forest from any v1 source the team has on Slack / a co-author tree; slide 09 can substitute the arm-leiden-k15 PDF directly. Decision needed from Erik before deck assembly.

---

## 2. Material — Slide 06 says "round out the six abstract-anchored clades" but enumerates only 5

**Slides involved:** 06 vs 07, 09, 13.

Slide 06 bullet 3 reads:
> "Acrocentric short arms (13p, 14p, 15p, 21p, 22p = **C7**) carry the heaviest, most uniform long-length distributions. PAR2 (Xq, Yq = **C14**) and PAR1 (Xp, Yp = **C15**) sit at the canonical pseudoautosomal scale. 4q and 10q (= **C1**, D4Z4 / DUX4) and the 10p–18p Linardopoulou pair (= **C2**) round out the **six** abstract-anchored clades."

That sentence lists **five** communities (C1, C2, C7, C14, C15) but says "six." The missing clade is **C6** — the tight q-arm clade {1q, 13q, 17q, 19q, 21q, 22q} — which slides 07, 09, and 13 all enumerate as the sixth abstract-anchored clade and which the abstract itself names verbatim ("a tightly linked clade involving 22q, 21q, 19q, 1q, 13q, and 17q").

C6 is **not visible as a distinct outlier on the slide-06 length-distribution facet grid** (its tail is moderate, not the clean "fat-right-tail" the slide narrates), which is presumably why slide 06's author left it out of the spoken walk. But the "round out the six" wording is internally inconsistent with that omission.

**Recommended resolution (lead-author call):** either (a) edit slide 06 bullet 3 to say "round out the **five** outlier clades visible on this panel; the sixth — C6, the q-arm clade — is moderate on length and only emerges in slide 09's PCA" (preserves the audience setup for slide 09), or (b) add a one-line note about C6 in the slide-06 callout legend. Option (a) is the lower-friction edit.

---

## 3. Material — PCA vs MDS labelling inconsistent across slides 08 / 09 / abstract

**Slides involved:** 08, 09; abstract C6.

The artifact is **classical MDS / PCoA** computed by `cmdscale(as.dist(jaccard_dist_df), eig=TRUE, k=5)` at `/moosefs/guarracino/HPRCv2/scripts/similarity/plot-similarity-subtelo.R` line 556. There is **no PCA artifact** on disk in the HPRCv2 tree (verified by Andrea's reports + cited in CROSSWALK §C6).

- **Slide 08 title:** *"All-vs-all in 2D — colored by chromosome (left) vs superpopulation (right)"* — body explicitly relabels the projection as MDS / PCoA and offers Erik two paths (a: relabel only; b: run a real PCA). Slide-08 author recommends path (a).
- **Slide 09 title:** *"All-vs-all PCA — 15 arm-level communities, named for the abstract's clades"* — still calls it "PCA" verbatim, with a body-note acknowledging the MDS/PCoA-vs-PCA gap and saying "for the talk this is fine — speaker should use 'PCA' verbally."
- **Abstract:** "**Principal component** and community detection analyses..." (the C6 wording).

The two adjacent slides (08 then 09) tell the same audience two different names for the same projection in 30 s. Audience confusion likely.

**Recommended resolution (lead-author call):** pick one terminology and propagate. Per slide-08's recommendation:
- **Path A (low friction):** change slide-09 title to "All-vs-all in 2D" or "All-vs-all MDS / PCoA — 15 arm-level communities…", drop "PCA" from the talk; rewrite the abstract's "Principal component and community detection analyses" to "MDS / PCoA and community detection analyses".
- **Path B (high friction):** run a real PCA on a haplotype × PHR-presence-binary feature matrix (REWRITE_PLAN TASK-16), generate fresh artifacts, keep "PCA" everywhere. Conclusion does not change; cost is engineering time.

For the BoG talk, **path A is the safe call**. The abstract is not locked (per CROSSWALK §C6) — relabel-only is the cheapest fix that makes 08 and 09 consistent.

---

## 4. Material — "Missing introvert arms" count differs between slide 11 and slides 05–09

**Slides involved:** 11 vs 05, 06, 07, 08, 09.

- Slides 05, 06, 07, 08, 09: **6** missing introvert arms — `2p, 3p, 5p, 8q, 11q, 14q` (no inter-chromosomal PHR detected, absent from the 41×41 matrix).
- Slide 11 negative-control "S_all" pseudo-community: **7** zero-sharing arms — adds `chr18_q` to the list above.

The two lists are not in conflict — chr18_q has n=1 in C15 (PAR1) per slide 09's clade table, so it is *not* an introvert arm in the matrix sense, but it *is* zero-sharing for the Dip-C / sperm pseudo-community on slide 11. Different definitions, different roles, but same audience hearing both lists 5 minutes apart.

**Recommended resolution:** add a half-second clarification in slide 11's speaker notes — "S_all extends the missing-introvert list with chr18_q because chr18_q has only n=1 sharing event, which is below the Dip-C window threshold." Do not edit the slide files (they are read-only per task spec); flag for the speaker.

---

## 5. Material — Hi-C numbers in slide 15 closer don't exactly match slide 10

**Slides involved:** 15 vs 10.

Slide 15's "Mechanism" bullet recaps:
> "Hi-C/Pore-C/CiFi/Dip-C all independently recover community-structured 3D contacts (B/W 0.027–0.074; Mantel ρ=0.296, p=0.002; per-pair ρ=0.674 in CHM13 Hi-C)."

Slide 10 (which slide 15 calls back to) cites:
- HG002 Pore-C bulk B/W = **0.056, p = 3.9 × 10⁻⁸⁵**
- Mantel ρ across 8 datasets ranging 0.15 → 0.66; specific samples HG002 0.66 → 0.79, CHM13 0.66 → 0.85 after exclusions
- Slide 11: GM12878 Dip-C Mantel ρ = **0.296, p = 3.8 × 10⁻⁴**

Slide 15's "B/W 0.027–0.074" comes from CROSSWALK §C7 (different per-haplotype range across samples; the 0.056 single number on slide 10 falls inside that range). Slide 15's "Mantel ρ=0.296, p=0.002" appears to be Dip-C / GM12878 from slide 11 — but slide 11 cites p = 3.8 × 10⁻⁴, not p = 0.002. The 0.002 value isn't directly traceable to a slide-10 or slide-11 number from the fanout files.

CROSSWALK §C7 chapter-05 numbers: HG002 B/W = **0.027 at 50 kb**, p = 4.0 × 10⁻⁶⁶; CHM13 B/W = **0.071** at 50 kb, p = 6.0 × 10⁻¹⁸. So slide 15's "0.027–0.074" is the per-haplotype range, not the slide-10 headline.

**Recommended resolution (lead-author call):** either (a) bring slide 15's headline numbers in line with slide 10's stated headlines (0.056 / 3.9e-85 / Mantel ρ 0.66 → 0.79), so the closer recap matches the slide the audience just saw; or (b) leave slide 15 as a recap of the *full* CROSSWALK / chapter-05 range and add a one-line gloss in slide 15's notes that "the 0.027–0.074 range covers per-haplotype values; the headline number on slide 10 was the cross-haplotype HG002 Pore-C single value." Option (a) is less surprising for the audience.

---

## 6. Polish — Slide 13 says "slides 5 and 9" for the Leiden-community context

**Slides involved:** 13.

Slide 13 speaker notes: "communities that show up in the family are exactly the communities we named in **slides 5 and 9**…"

The clade vocabulary is **introduced in slide 06** (named outliers on the length-distribution panel), **defined in slide 07** (NJ tree + heatmap), and **enumerated in slide 09** (PCA keystone with the abstract-anchored legend). Slide 05 is the PHR-scale slide; it does **not** name communities.

**Recommended resolution:** change to "slides 7 and 9" (NJ tree + PCA keystone) or "slides 6, 7, and 9" (the full naming arc). Polish-level — the audience won't trip on it during the talk, but the speaker notes drift from the actual deck order. Speaker should mentally substitute when delivering.

---

## 7. Polish — Slide 06 says "as you'll see on slide 9" — slide 9 is the keystone

**Slides involved:** 06 → 09.

Slide 06's reframe bullet: "the shape of the per-arm length distribution **recapitulates the community partition you are about to see on slide 09**."

Slide 09 is correctly identified as the PCA keystone — but the deck has slide 07 (heatmap + NJ tree) and slide 08 (PCA-by-chrom + superpop) **between** slide 06 and slide 09. The speaker should pre-cue this with something like "you'll see this confirmed three slides from now on the keystone PCA" rather than "the next slide."

**Recommended resolution:** speaker delivers, no edit. Speaker should say "you'll see this confirmed in a few slides on the keystone PCA" or similar to preserve the forward-pointer without claiming it's immediate.

---

## 8. Polish — DUX4 number propagation across slides 13/14/15

**Slides involved:** 13, 14, 15.

Three different DUX4 numbers, all internally consistent but emphasising different things:
- **Slide 13**: chr4q→chr10q gene-conversion in PAN028 maternal at **score 0.957** (a single observed event in the family)
- **Slide 14**: median **22 DUX4L copies** per haplotype on chr4q+chr10q in C1 (population-level copy-number distribution)
- **Slide 15**: cites the slide-14 number (median 22, Mann-Whitney p=5.3e-6); does **not** cite the slide-13 0.957

The pedigree event (0.957) and the population copy-number median (22) are different facts about the same locus, so the differing emphases are correct. Slide 15's omission of the 0.957 is a small narrative-callback miss — slide 13's punchline is exactly that locus, and slide 15 could close the loop with one phrase.

**Recommended resolution (polish):** speaker can say in slide 15 delivery: "DUX4 — the disease-revealed instance, where we caught a 4q→10q gene conversion in the family at score 0.957." Speaker delivery only; do not edit slides.

---

## 9. Polish — Time budget overrun (980 s vs 900 s target — see SLIDES_v2_PLAN.md §1)

**Slides involved:** all.

Sum of per-slide time budgets is 980 s (16:20). Target is ≤ 900 s (15:00). Overrun = 80 s.

Per-slide budgets, in order: 30, 50, 80, 70, 60, 50, 80, 70, 80, 80, 60, 60, 90, 50, 70.

This is not a coherence inconsistency per se — the budgets are internally consistent and each slide author justified their own time — but the deck does not fit a 15-minute slot as currently specified. Triage candidates ranked by slide-author-flagged compressibility:
1. **Slide 14** (50 s) — compressible to **0** per slide-14 author. Drops 50 s, leaves 30 s overrun.
2. **Slide 13** (90 s) — could be cut to 70 s by dropping the CEPH1463 paragraph per slide-13 author. Drops 20 s.
3. **Slide 03** (80 s) — could be cut to 65 s by compressing the ER-callout discussion. Drops 15 s.

Cuts (1)+(2) alone bring the deck to 910 s — still 10 s over but inside Q&A round-off. (1)+(2)+(3) brings it to 895 s, comfortably inside. See `SLIDES_v2_PLAN.md` §1 for the detailed decision matrix.

**Recommended resolution (lead-author call):** decide which of the three cuts to commit to. Default recommendation: cut slide 14 entirely if the talk runs long live; pre-rehearse with slide 13's CEPH1463 paragraph as a soft-cut.

---

## 10. Polish — slide 08 axes labels: "Dimension 1 (16.05%)" vs slide 09 "PC1 (16.05%)"

**Slides involved:** 08, 09.

Slide 08's relabel decision (MDS/PCoA per §3 above) means axes should read "Dimension 1" / "Dimension 2" with variance-explained percentages. Slide 09 inherits the v1 deck PCA scatter axes ("PC1 / PC2"). If §3 path A is chosen, slide 09 should be relabeled too — otherwise the same scatter shows different axis labels on consecutive slides.

**Recommended resolution:** tied to §3 decision. If path A, relabel both slides "Dim 1 / Dim 2 (variance-explained)"; if path B, relabel both "PC1 / PC2 (variance-explained)" after re-rendering. Either way, slides 08 and 09 must agree on axis labels.

---

## 11. Polish — Mouse meiotic figure-vs-inset coherence (slide 12)

**Slides involved:** 12.

Slide 12 instructs the synthesizer to place the four-stage trajectory inset (lepto/zygo/pachy/diplo Mantel ρ trajectory) **next to** Fig 4d (zygotene scatter). Without the inset, the audience hears "zygotene = 0.715" and has to trust the peak claim. The inset is **not yet rendered** in `slides/v2/` — see `figure_manifest.md`.

**Recommended resolution:** flagged to slide-12 author and synthesizer; render the inset (`Rscript slide_12_meiotic_stage_trajectory.R` from the slide content) before deck assembly. This is a load-bearing render, not optional polish.

---

## Summary — narrative gaps and callbacks

The narrative arc (see `SLIDES_v2_PLAN.md` §2) **does land**: every callback in the closer (slide 15) is supported by an earlier slide; every clade enumerated on slide 09 is substantively defined on slide 06 or slide 07; the 3D arc (slides 10 → 11 → 12) builds cleanly from bulk → single-cell → meiotic → pedigree (slide 13 proof). Items §1–§5 are concrete consistency issues to resolve before the talk; §6–§11 are polish-level.

**Decisions blocking finalization (in order of urgency):**
1. v1 deck PDF restoration (§1)
2. Slide 06 "six clades" vs five-listed (§2)
3. PCA-vs-MDS terminology (§3) — and the cascade to axis labels (§10)
4. Slide-15 Hi-C numbers vs slide-10 headlines (§5)
5. Time budget triage (§9)
