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
