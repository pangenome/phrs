// BoG 2026 v2 deck — built by build-bog-v2-2 (agent-813) from slides/v2/slide_NN_*.md
// Compile: typst compile --root /moosefs/erikg/phrs/.wg-worktrees/agent-813 \
//          slides/v2/_typst/deck.typ slides/v2/BoG_2026.pdf
// Speaker: Erik Garrison · Biology of Genomes 2026 (CSHL).

#import "@preview/polylux:0.4.0": *

// ---------- Page + typography ------------------------------------------------
// Standard 16:9 talk format (10in × 5.625in ≈ 25.4cm × 14.29cm).
#set page(
  width:  25.4cm,
  height: 14.29cm,
  margin: (x: 1.0cm, y: 0.7cm),
)
#set text(font: ("Liberation Sans", "DejaVu Sans"), size: 14pt)
#set par(leading: 0.55em)
#show heading.where(level: 1): it => block(
  below: 0.4em,
  text(weight: "bold", size: 22pt, fill: rgb("#1b3a6f"), it.body),
)
#show heading.where(level: 2): it => block(
  below: 0.3em,
  text(weight: "bold", size: 16pt, fill: rgb("#1b3a6f"), it.body),
)

// ---------- Helpers ----------------------------------------------------------
#let abs_root = "/moosefs/erikg/phrs/.wg-worktrees/agent-813"

#let bog-slide(title: none, time: none, body) = slide[
  // time: budget comment, e.g. "30s"
  #if title != none {
    block(
      stroke: (bottom: 0.5pt + rgb("#1b3a6f")),
      inset: (bottom: 4pt),
      below: 8pt,
      text(weight: "bold", size: 18pt, fill: rgb("#1b3a6f"), title),
    )
  }
  #body
  #if time != none {
    place(
      bottom + right,
      dx: -0.2cm,
      dy: -0.1cm,
      text(size: 8pt, fill: gray, [time budget: #time]),
    )
  }
]

#let bullets(..items) = list(
  marker: ([▪], [•]),
  spacing: 0.55em,
  ..items.pos().map(it => it),
)

// Speaker notes live inside `/* SPEAKER … */` block comments below each slide.
// They are not rendered but are preserved in this .typ source for Erik.

// =============================================================================
// SLIDE 01 — Title
// time: 30s
// =============================================================================
#slide[
  #set align(center + horizon)
  #text(size: 28pt, weight: "bold", fill: rgb("#1b3a6f"))[
    Concerted evolution and unorthodox recombination
  ]
  #v(0.3em)
  #text(size: 28pt, weight: "bold", fill: rgb("#1b3a6f"))[
    of human subtelomeres
  ]
  #v(1.2em)
  #text(size: 18pt)[Andrea Guarracino · Erik Garrison]
  #v(0.6em)
  #text(size: 14pt, fill: gray)[Companion to HPRC v2 (Nature, in submission)]
  #v(0.4em)
  #text(size: 14pt)[Inter-chromosomal subtelomeric relationships at HPRC v2 scale —]
  #linebreak()
  #text(size: 14pt)[466 near-complete haplotypes]
  #v(0.8em)
  #text(size: 13pt, fill: rgb("#444444"))[Biology of Genomes · Cold Spring Harbor · May 2026]
]
/* SPEAKER
Title slide. Land the full title once, then translate it for a BoG audience that
  includes folks who don't work on subtelomeres: what we actually did is survey
  inter-chromosomal relationships between subtelomeres at population scale, and
  what we found is that this looks like concerted evolution — gene-conversion-like
  and crossover-like exchange between non-homologous chromosomes — happening more
  broadly than previously appreciated. Frame the talk as the companion to HPRC v2:
  same data, new question. Co-authored with Andrea Guarracino. Manuscript in
  submission to Nature. Fifteen-minute talk; thirty-second slide; move quickly.
*/

// =============================================================================
// SLIDE 02 — Implicit interval tree
// time: 50s
// =============================================================================
#bog-slide(
  title: [Implicit interval tree — the data structure under the implicit pangenome graph],
  time: [50s],
)[
  #grid(
    columns: (1.4fr, 2.2fr),
    column-gutter: 12pt,
    [
      #set text(size: 11pt)
      #bullets(
        [Each pairwise alignment becomes an interval `[start, end)` on its target sequence (CIGAR + target range)],
        [Li \& Rong (2020) *implicit interval tree*: sorted array, tree topology implied by index — O(log n) overlap queries, near-zero memory overhead],
        [One interval tree per target sequence; the union over all sequences is an *interval forest*],
        [That forest, queried by transitive closure, *is* our implicit pangenome graph — no graph is ever materialized],
        [Foundation for the rest of the talk: every plot you see comes from queries against this structure],
      )
    ],
    align(center + horizon, image("v1_page_02-02.png", height: 10cm)),
  )
]
/* SPEAKER
One slide of plumbing before the biology — the rest of the talk is just queries
  against this object. Li & Rong 2020: implicit interval tree, O(log n) overlap,
  zero pointer memory. One tree per sequence → interval forest. IMPG (github.com/
  pangenome/impg) does the transitive closure. Keep this picture in your head —
  it's the substrate for every plot that follows.
*/

// =============================================================================
// SLIDE 03 — IMPG workflow (the implicit pangenome graph)
// time: 80s
// =============================================================================
#bog-slide(
  title: [The implicit pangenome graph — wfmash all-vs-all *is* the graph],
  time: [80s],
)[
  #grid(
    columns: (1.4fr, 2.0fr),
    column-gutter: 10pt,
    [
      #set text(size: 10.5pt)
      #bullets(
        [All-vs-all wfmash (`-p 95`) over n = 18,827 telomere-anchored 500 kb flanks → 18,827 PAFs. The union of edges *is* the implicit pangenome graph — no GFA, no construction step.],
        [Index = interval forest, one implicit interval tree per sequence (Li \& Rong 2020). O(n log n) build, O(log n + k) interval queries.],
        [Query = `impg query -x` (transitive closure) — chase chains of overlapping mappings outward from any seed; the graph is never materialized.],
        [Sampling is dense, not sparse. wfmash evaluates ~12% of C(n,2) ≈ 1.77×10⁸ pairs at full alignment cost. ER threshold p\* = log(n)/n ≈ 5.21×10⁻⁴; *12% is ~230× above p\** → densely connected w.h.p., closure reaches genome-wide.],
        [Cite IMPG: https://github.com/pangenome/impg (Garrison et al.)],
      )
    ],
    [
      #image("v1_page_03-03.png", height: 6.5cm)
      #v(0.2em)
      #align(right, image("../slide_03_er_callout.png", height: 4cm))
    ],
  )
]
/* SPEAKER
Methods anchor — slow down. We do NOT build a GFA. wfmash all-vs-all on 18,827
  telomere-anchored 500 kb flanks; the resulting PAF set IS the graph. Index it as
  an interval forest. Query via IMPG transitive closure. Erdős-Rényi threshold for
  n=18,827 is p*=log(n)/n ≈ 5.2e-4; 12% sampling is 230× above that → graph is
  densely connected. That's what licenses "no chromosomal partitioning."
*/

// =============================================================================
// SLIDE 04 — Genome-wide identity heatmap
// time: 70s
// =============================================================================
#bog-slide(
  title: [Genome-wide identity heatmap — interchromosomal homology at PAR2 scale],
  time: [70s],
)[
  #grid(
    columns: (1.2fr, 2.2fr),
    column-gutter: 12pt,
    [
      #set text(size: 11pt)
      #bullets(
        [466 HPRCv2 haplotypes vs CHM13, per-position max identity to any matching chromosome (100 kb windows)],
        [Most of each chromosome is silent; dense red bands appear where assemblies reach the telomeres],
        [chr18 q-arm inset: tight subtelomeric block of >98% identity matches from many chromosomes],
        [These blocks span *10s–100s of kb — comparable in scale to PAR2 (~334 kb on Xq/Yq)*],
        [PAR2-scale pseudo-homology at nearly every subtelomere → motivates the all-vs-all view that follows],
      )
    ],
    align(center + horizon, image("/paper_prep/figures/fig1/figure_fig1.png", height: 10cm)),
  )
]
/* SPEAKER
Empirical foundation. For every 100 kb window across 466 haplotypes we plot max
  alignment identity to any other chromosome. Most chromosomes silent; telomeres
  light up. chr18q inset: tight band of >98% identity hits from many chromosomes,
  spanning 10s–100s kb. Reframe — PAR2 is ~334 kb. We are seeing PAR2-scale
  pseudo-homology at nearly every subtelomere.
*/

// =============================================================================
// SLIDE 05 — Interchromosomal similarities (n-chromosomes per region)
// time: 60s
// =============================================================================
#bog-slide(
  title: [Interchromosomal similarities — n-chromosomes per region (HPRCv2)],
  time: [60s],
)[
  #grid(
    columns: (1.1fr, 2.4fr),
    column-gutter: 12pt,
    [
      #set text(size: 10.5pt)
      #bullets(
        [Orange traces: number of unique chromosomes matching each 100 kb window across CHM13; CEN / PAR / PHR / XTR painted as background reference],
        [Spikes are not confined to PARs — every subtelomere lights up, plus centromeres and the acrocentric short arms],
        [*PHR scale: median 105 kb, mean 144 kb, range 5–500 kb* (15,668 PHR sequences across 41/48 chromosome arms)],
        [Same length scale as PAR2 (~334 kb) — but present at nearly every chromosome end, not just X/Y],
        [Read the plot as: a *PAR2-class exchange landscape replicated 41 times* across the human genome],
      )
    ],
    align(center + horizon, image("/p_num_chromosomes_wide.png", height: 10cm)),
  )
]
/* SPEAKER
Same view as last time — number of unique chromosomes per 100 kb window across
  CHM13. 466 haplotypes → 15,668 PHRs spanning 41 of 48 arms. Median length 105
  kb, mean 144, range 5–500. PAR2 ≈ 334 kb. So a typical subtelomeric PHR is on
  the same order as PAR2 — and unlike PAR2, this is happening at nearly every
  chromosome end. Extended pseudohomology at nearly all subtelomeres, comparable
  in scale to canonical PARs, but replicated dozens of times across the genome.
*/

// =============================================================================
// SLIDE 06 — Length distributions per arm + named clades callout
// time: 50s
// =============================================================================
#bog-slide(
  title: [Length distributions per arm — the outliers are named clades],
  time: [50s],
)[
  #grid(
    columns: (2.0fr, 1.4fr),
    column-gutter: 10pt,
    [
      #image("img/slide_06_length_dist-1.png", height: 11.5cm)
    ],
    [
      #set text(size: 9pt)
      #bullets(
        [Faceted histograms (v1 slide 6, verbatim): one panel per arm; *p-arm = blue, q-arm = orange*; pink = introvert arms (2p, 3p, 5p, 8q, 11q, 14q)],
        [Bulk arms cluster around the population PHR scale — *median 105 kb, mean 144 kb*],
        [*The fat-right-tail outliers are the abstract's clades* (C1, C2, C7, C14, C15)],
        [Pink panels are biological signal, not missing data — arms with no cross-chrom hits],
        [*The shape of these distributions already recapitulates the community partition on slide 9*],
      )
      #v(0.4em)
      #image("../slide_06_clade_callouts.png", width: 100%)
    ],
  )
]
/* SPEAKER
Same panel as last time — one histogram per chromosome arm, p blue / q orange,
  pink for introvert arms with zero cross-chromosomal hits. Outlier facets =
  abstract's clades. C7 acrocentric (13p, 14p, 15p, 21p, 22p — fully homogenized),
  C14 PAR2 (Xq/Yq), C15 PAR1 (Xp/Yp), C1 4q–10q DUX4, C2 10p–18p Linardopoulou.
  The histograms already know the clades before we cluster.
*/

// =============================================================================
// SLIDE 07 — All-vs-all heatmap + NJ tree (the cladistic structure)
// time: 80s
// =============================================================================
#bog-slide(
  title: [All-vs-all at the arm level — heatmap + NJ tree, with named clades],
  time: [80s],
)[
  #grid(
    columns: (1fr, 1fr),
    column-gutter: 6pt,
    align(center + horizon, image("/paper_prep/figures/fig1/figure_fig1.png", height: 10cm)),
    align(center + horizon, image("/paper_prep/figures/nj_tree_arms/nj_tree_annotated.png", height: 10cm)),
  )
  #v(-0.3em)
  #align(center, text(size: 8pt, fill: gray, style: "italic")[
    Same 41×41 arm-level Jaccard distance matrix; left = Fig 1 (panel c) clustered heatmap (Leiden k=15 + UPGMA), right = NJ tree (rooted at acrocentric MRCA, 1000-rep bootstrap). Six abstract-named clades are monophyletic with 100% bootstrap.
  ])
]
/* SPEAKER
All-vs-all at the arm level. 41×41 matrix; cyan boxes = 15 Leiden communities;
  UPGMA dendrogram on top recovers 14/15. Right panel: neighbor-joining tree on
  the SAME matrix, rooted at acrocentric MRCA. Every clade name in the abstract —
  PAR1, PAR2, acrocentric short arms, 10p–18p, 4q–10q DUX4, the tight q-arm clade
  (22q/21q/19q/1q/13q/17q) — is monophyletic with 100% bootstrap. Three
  algorithms (Leiden, UPGMA, NJ) → same six clades. Cladistic structure is real.
*/

// =============================================================================
// SLIDE 08 — MDS / PCoA, by chromosome vs by superpopulation
// time: 70s
// =============================================================================
#bog-slide(
  title: [All-vs-all in 2D — by chromosome (left) vs by superpopulation (right)],
  time: [70s],
)[
  #grid(
    columns: (1fr, 1fr),
    column-gutter: 6pt,
    align(center + horizon, image("img/slide_08_mds_chrom.png", height: 7.5cm)),
    align(center + horizon, image("img/slide_08_mds_superpop.png", height: 7.5cm)),
  )
  #v(0.1em)
  #set text(size: 8.5pt)
  #bullets(
    [Identical points (n = 15,668 flanks); left colored by source chromosome, right by 1KGP superpopulation (AFR/AMR/EAS/EUR/SAS)],
    [Left clusters are *arm-community shaped* (D4Z4 4q+10q, acrocentric p, PAR1, q-arm clade); right shows superpopulations *mixed across all clusters*],
    [Population structure is real but secondary: Hudson Fst mean 0.044, AFR vs non-AFR 0.10–0.15 — *arm-community first, population structure within*],
  )
]
/* SPEAKER
Same all-vs-all object as the heatmap, projected into 2D — classical MDS / PCoA
  on the Jaccard distance matrix. One scatter, two colorings. Left: clusters are
  arm-community shaped — D4Z4 (4q/10q), acrocentric p-arms, PAR1, the q-arm
  clade. Right: same points, recolored by superpop — colors mix across every
  cluster. Population structure is real but secondary. Hudson Fst mean 0.044,
  AFR-vs-non-AFR 0.10–0.15. Arm-community first, population within.
*/

// =============================================================================
// SLIDE 09 — PCA by community (the keystone slide)
// time: 80s
// =============================================================================
#bog-slide(
  title: [PCA — 15 arm-level communities, named for the abstract's clades],
  time: [80s],
)[
  #grid(
    columns: (1.1fr, 1.4fr),
    column-gutter: 10pt,
    align(center + horizon, image("v1_page_10-10.png", height: 11cm)),
    [
      #set text(size: 8pt)
      #bullets(
        [Same PCA layout as v1: each point is an arm; clusters = Leiden k=15 communities],
        [*Every named clade in the abstract is a community on this plot* — PAR1=C15, PAR2=C14, acrocentric=C7, 22q-clade=C6, 10p–18p=C2, DUX4=C1],
        [Three zones: *PARs* (C14+C15), *concerted-exchange core* (C6+C7+C3), *4q–10q DUX4* (C1)],
        [PC1 16.05% / PC2 11.2%; UPGMA k=14 recovers 12/15 — partition is not a Leiden artefact],
        [*Everything that follows is evidence about clusters you can already see here*],
      )
      #v(0.2em)
      #image("../slide_09_clade_legend.png", height: 6cm)
    ],
  )
]
/* SPEAKER
Keystone of the talk. Every clade word in the abstract is a colored cluster.
  PCA on the 41×41 Jaccard arm distance matrix. Walk the audience around in
  abstract order: PAR2 = C14 (Xq/Yq), PAR1 = C15 (Xp/Yp), acrocentric = C7 (all
  five p-arms), the novel q-arm clade = C6 (1q/13q/17q/19q/21q/22q exact match),
  10p–18p = C2 Linardopoulou, 4q–10q DUX4 = C1. The abstract's vocabulary maps
  onto a single empirical structure; everything that follows is evidence for it.
*/

// =============================================================================
// SLIDE 10 — Hi-C / Pore-C bulk + Mantel exclusions
// time: 80s
// =============================================================================
#bog-slide(
  title: [Hi-C / Pore-C confirm communities are 3D — signal *strengthens* with confound exclusions],
  time: [80s],
)[
  #grid(
    columns: (1fr, 1fr),
    column-gutter: 6pt,
    align(center + horizon, image("/paper_prep/figures/fig3/figure_fig3.png", height: 7.5cm)),
    align(center + horizon, image("/paper_prep/figures/ed5/figure_ed5.png", height: 7.5cm)),
  )
  #v(0.1em)
  #set text(size: 8.5pt)
  #bullets(
    [*HG002 Pore-C contact matrix*, 50 kb, 77 arm-haplotypes ordered by sequence community: B/W = 0.056, p = 3.9×10⁻⁸⁵],
    [*Bulk Mantel ρ positive in 7/8 datasets* — CHM13 ρ=0.66, HG002 Hi-C ρ=0.66, Pore-C ρ=0.49],
    [*Exclude acrocentric+sex+strong communities → ρ goes UP*: HG002 0.66→0.79, CHM13 0.66→0.85, *HG02148 0.15→0.72* (marginal becomes strong)],
    [At 10 kb, community-free per-pair correlations reach ρ = 0.83 (NA19036), 0.81 (HG02148)],
    [Sequence similarity predicts 3D contact — and the signal *strengthens* after confound exclusions, so it's *generic subtelomeric homology*, not nucleolar / PAR clustering],
  )
]
/* SPEAKER
From sequence to nucleus. Left: HG002 Pore-C inter-arm contacts at 50 kb,
  ordered by sequence community. Diagonal blocks light up — within-community
  contacts 18× higher than between, p ≈ 1e-85. Right: bulk Mantel similarity ×
  contact, full matrix vs (no acrocentric + no sex). 7/7 above identity. HG002
  0.66→0.79, CHM13 0.66→0.85, HG02148 0.15→0.72. Generic subtelomeric homology
  drives 3D contact, not nucleolar artifact.
*/

// =============================================================================
// SLIDE 11 — Single-cell 3D (Dip-C + sperm)
// time: 60s
// =============================================================================
#bog-slide(
  title: [Single-cell 3D — and it works in haploid sperm],
  time: [60s],
)[
  #grid(
    columns: (1.5fr, 1.5fr),
    column-gutter: 10pt,
    align(center + horizon, image("/paper_prep/figures/fig3/figure_fig3.png", height: 11cm)),
    [
      #set text(size: 9pt)
      #bullets(
        [*GM12878 Dip-C, 16 cells (Tan 2018, T2T-CHM13v2.0):* community arms 6.9% closer in 3D within community than between (W/B = 0.931, Wilcoxon p = 3.8×10⁻⁴; Mantel ρ = 0.296)],
        [*Sperm scHi-C, 20 cells (Xu 2025): 60% closer* within community (W/B = 0.401, Fisher p = 3.9×10⁻⁵¹) — in the haploid, hyper-condensed sperm nucleus],
        [*Negative control* — pseudo-community S of 7 zero-sharing arms (2p, 3p, 5p, 8q, 11q, 14q, 18q): *11% farther* in GM12878, *40% farther* in sperm],
        [Same pattern in diploid soma + haploid germline → not a Hi-C artefact, not restricted to interphase],
        [*Sperm is the bridge* — puts the signal in the gamete that carries the recombination, sets up mouse meiotic Hi-C next],
      )
    ],
  )
]
/* SPEAKER
Bulk Hi-C could be a population artefact. Single-cell rules that out. Tan 2018
  Dip-C, GM12878 (16 cells, T2T-remapped): community arms ~7% closer in 3D within
  community. Mantel ρ = 0.30. Now the punchline: 20 sperm cells from Xu 2025.
  Haploid, hyper-condensed. W/B = 0.40 — 60% closer; Fisher p ≈ 4e-51. Negative
  control: 7 zero-sharing arms (S_all) move the *opposite* way, +11% / +40%
  farther. Sequence sharing necessary; community label alone insufficient.
*/

// =============================================================================
// SLIDE 12 — Mouse meiotic Hi-C — zygotene bouquet peak
// time: 60s
// =============================================================================
#bog-slide(
  title: [Mouse meiotic Hi-C — the zygotene bouquet is where the 3D signal peaks],
  time: [60s],
)[
  #grid(
    columns: (1.4fr, 1.4fr),
    column-gutter: 10pt,
    align(center + horizon, image("/paper_prep/figures/fig4/figure_fig4.png", height: 11cm)),
    [
      #image("../slide_12_stage_trajectory.png", height: 5.5cm)
      #v(0.1em)
      #set text(size: 7.5pt)
      #bullets(
        [Bulk human Hi-C is *mitotic*; the recombination is *meiotic*. Mouse zygotene Hi-C (Zuo 2021) is the only meiotic 3D map on a T2T-grade genome — 4 stages: lepto/zygo/pachy/diplo.],
        [*Mantel ρ peaks at zygotene*: 0.687 / *0.718* / 0.683 / 0.577. Per-PHR-pair Spearman ρ = 0.715, p = 4.4×10⁻⁵⁵, n = 344],
        [*Zygotene = the bouquet stage* — telomeres clustered at the LINC-anchored nuclear envelope to align homologs (Mefford 2002 / Linardopoulou 2005)],
        [Mouse is telocentric: 39 p-flanks have signal, 39 q-flanks have zero — the pipeline finds the right end],
        [*Human LCL Hi-C is the somatic shadow of a meiotic phenomenon — visible directly in mouse, at the bouquet*],
      )
    ],
  )
]
/* SPEAKER
The skeptic: bulk Hi-C is mitotic; recombination is meiotic. Zuo 2021 sorted 4
  prophase stages: lepto / zygo / pachy / diplo on a T2T-grade mouse. Run the
  same Mantel: ρ = 0.687 / 0.718 / 0.683 / 0.577. Zygotene peaks. Per-PHR Spearman
  at zygotene ρ = 0.715, p = 4.4e-55, n = 344. Zygotene is the bouquet —
  telomeres clustered at LINC-anchored nuclear envelope while homologs align.
  Cross-species generality is a bonus.
*/

// =============================================================================
// SLIDE 13 — WashU pedigree — 92% within community
// time: 90s
// =============================================================================
#bog-slide(
  title: [Caught in the act — three generations of a T2T pedigree show ongoing exchange],
  time: [90s],
)[
  #grid(
    columns: (1.6fr, 1.4fr),
    column-gutter: 10pt,
    align(center + horizon, image("img/slide_13_untangle-1.png", height: 10.5cm)),
    [
      #set text(size: 8pt)
      #bullets(
        [*WashU T2T pedigree* (Cechova 2025), 3 generations, 4 individuals: PAN010 → PAN027 → PAN028 + paternal grandfather PAN011. Every haplotype T2T; *odgi untangle* compares each child flank against its parent.],
        [*538 high-quality inter-chromosomal patches; 494 (92%) fall inside HPRC-v2 Leiden communities* (built from 233 unrelated samples — slide 9). The graph predicts where exchange shows up; the family delivers the events.],
        [*133 ectopic gene-conversion-like sandwich tracts* at score ≥ 0.81; 96 at the perfect 1.000/1.000 ceiling. *16 crossover-like* events (left + right flanks resolve to *different* haplotypes — meiotic crossover signature).],
        [*C7 acrocentric traffic dominates*; non-acrocentrics land in the named clades: chr18p↔chr10p (C2), chr3q↔chr9q (C3), chrXp↔chrYp (C15 PAR1), and *DUX4 chr4q→chr10q at 0.957 in PAN028 maternal* (C1)],
        [*CEPH1463 replication* (Porubsky 2025): 11 parent features in same Leiden communities, confirmed by *both* hifiasm and verkko — different family, two assemblers, same communities],
      )
    ],
  )
]
/* SPEAKER
Direct empirical proof. WashU pedigree, 4 T2T individuals, 3 generations. Push
  each child haplotype through odgi untangle vs its parent. The figure: PAN027
  maternal hap1 painted onto PAN010. Diagonal stripes = correctly inherited.
  Off-color stripes = inter-chromosomal patches — gene-conversion-like events.
  538 HQ patches; 494 (92%) inside Leiden communities built from 233 unrelated
  samples. 133 gene-conversion-like, 96 at perfect score, 16 crossover-like.
  DUX4 chr4q→chr10q at 0.957 in PAN028 maternal — the disease locus, in a normal
  family. CEPH1463: 11 features confirmed by both assemblers, same communities.
  Concerted evolution caught in the act, three generations of one family.
*/

// =============================================================================
// SLIDE 14 — Gene biology aside (DUX4, OR4F, TAR1)
// time: 50s (compressible to 0)
// =============================================================================
#bog-slide(
  title: [Gene biology aside — DUX4, OR4F, TAR1 (the biology is interesting too)],
  time: [50s],
)[
  #grid(
    columns: (1.6fr, 1.0fr),
    column-gutter: 10pt,
    align(center + horizon, image("../slide_14_gene_biology.png", height: 8.5cm)),
    [
      #set text(size: 9pt)
      #bullets(
        [*DUX4 (FSHD locus).* Annotated on 18 q-arms, but only chr4q and chr10q (community *C1*) carry the full *D4Z4 macrosatellite* (median *22* DUX4L copies; FSHD-permissive 4qA lives here). Elsewhere: 0–2 copies.],
        [*OR4F (olfactory receptors).* 4 paralogs span 16 arms, 5,023 entries. Pseudogenisation runs as a clean per-arm gradient: *11.1% (chr7p) → 99.8% (chr15q)*; population mean 62.1%.],
        [*TAR1 (telomere-associated repeat).* 21,544 entries across *94.6%* of all 15,668 PHR sequences and all 41 arms — universal *except at PAR1* (Xp/Yp, 0.5%). PAR1 has obligate meiotic crossover; no satellite anchor needed.],
        [*Distinct biological histories — disease, decay, exchange machinery — write themselves into the same subtelomeric architecture*],
      )
    ],
  )
]
/* SPEAKER
Quick aside before the close. Three vignettes, ten seconds each. DUX4 — 18
  q-arms have it, but only chr4q + chr10q carry the FSHD-relevant D4Z4 array
  (median 22 copies). OR4F — pseudogenisation runs 11% → 99.8% across the
  subtelomeric exchange network; same gene, different decay clock. TAR1 — present
  at 94.6% of PHRs across 41 arms, except PAR1 (0.5%) which has obligate meiotic
  crossover and doesn't need a satellite anchor. Disease, decay, exchange
  machinery — all three write themselves into the same subtelomeric architecture.
  Skippable if running long.
*/

// =============================================================================
// SLIDE 15 — Concerted evolution thesis (closer)
// time: 70s
// =============================================================================
#bog-slide(
  title: [Concerted evolution of human subtelomeres — what we saw, predicted, and recovered],
  time: [70s],
)[
  #grid(
    columns: (1.4fr, 1.4fr),
    column-gutter: 10pt,
    [
      #set text(size: 9.5pt)
      #bullets(
        [*Method (slide 3).* Implicit pangenome graph: wfmash all-vs-all over 18,827 telomere-anchored flanks; ~12% of pairs evaluated — *230× above the Erdős-Rényi threshold* (p\* ≈ 5.21×10⁻⁴). No partitioning, no GFA.],
        [*Empirical (slides 4–9).* *15,668 PHRs across 41/48 arms*, median 105 kb / mean 144 kb — PAR2-scale pseudohomology at nearly every chromosome end. Named clades: PARs (Xp/Yp, Xq/Yq), acrocentrics (C7), 10p–18p, the q-arm clade *22q–21q–19q–1q–13q–17q*, and *4q–10q DUX4*.],
        [*Mechanism (slides 10–12).* Hi-C/Pore-C/CiFi/Dip-C all recover community-structured 3D contacts (B/W 0.027–0.074; Mantel ρ=0.296, p=0.002; CHM13 per-pair ρ=0.674). Median PHR sits at the base of a single meiotic loop — the *bouquet*.],
        [*Proof (slide 13).* *WashU pedigree: 494/538 (92%) inter-chromosomal patches inside Leiden communities* — 133 gene-conversion-like, 16 crossover-like.],
        [*Biology (slide 14).* D4Z4 / DUX4 (4q↔10q) is the disease-revealed instance — FSHD-modifying translocations are *concerted evolution caught in the act*.],
      )
    ],
    align(center + horizon, image("/paper_prep/figures/ed8/figure_ed8.png", height: 8.5cm)),
  )
  #v(0.2em)
  #block(
    fill: rgb("#FFF7E6"),
    stroke: 0.5pt + rgb("#1b3a6f"),
    inset: 6pt,
    radius: 3pt,
    width: 100%,
    text(size: 10pt, weight: "bold", fill: rgb("#1b3a6f"))[
      Thesis: subtelomeres concertedly evolve through ongoing inter-chromosomal exchange — observable in pedigrees, predicted by 3D, recovered by an implicit pangenome graph across HPRC v2.
    ],
  )
]
/* SPEAKER
Close. One slide, one breath, one thesis. Method: implicit pangenome graph, 230×
  above ER threshold. Empirical: 15,668 PHRs across 41/48 arms, named clades.
  Mechanism: 3D community structure across Hi-C / Pore-C / Dip-C / sperm; bouquet
  is the predicted exchange venue. Proof: 92% of WashU pedigree inter-chrom
  patches inside Leiden communities. Biology: DUX4 is the disease-revealed
  instance. Thesis verbatim. Thank you.
*/
