// BoG 2026 — Annotated Review Deck
// Layout: 16:9 widescreen, figure LEFT, explanation RIGHT
// One logical slide per page (composite slides split into two pages)
// Build: typst compile review_deck.typ ../BoG_2026_review.pdf

#set page(
  width: 13.33in,
  height: 7.5in,
  margin: (x: 0.42in, y: 0.32in),
)
#set text(size: 11pt, lang: "en")
#set par(justify: true, leading: 0.6em)

// --- color palette ---
#let col-title   = rgb("#1a3a6b")
#let col-why-bg  = rgb("#fff8e8")
#let col-why-bar = rgb("#d4820a")
#let col-note-bg = rgb("#f7f7f7")
#let col-hdr-bg  = rgb("#dce8f7")
#let col-cap     = rgb("#888888")

// --- helper: compact header bar ---
#let page-header(num, label) = block(
  fill: col-hdr-bg,
  width: 100%,
  inset: (x: 0.08in, y: 0.04in),
  radius: 2pt,
)[#text(size: 7.5pt, fill: col-title.lighten(20%))[*Slide #num* #h(0.15in) #label]]

// --- helper: "why this slide matters" callout ---
#let why-box(body) = block(
  fill: col-why-bg,
  width: 100%,
  inset: (x: 0.1in, y: 0.08in),
  radius: 3pt,
  stroke: (left: 3pt + col-why-bar),
)[#text(size: 9pt)[#strong[Why this slide matters: ]#body]]

// --- normal two-column slide layout ---
#let slide(
  num, label, fig-path, title, notes, why,
  caption: "",
  fig-fraction: 55%,
) = {
  page-header(num, label)
  v(0.04in)
  grid(
    columns: (fig-fraction, 100% - fig-fraction),
    gutter: 0.22in,
    // LEFT: figure
    align(horizon + center)[
      #image(fig-path, fit: "contain", width: 100%, height: 5.85in)
      #if caption != "" {
        v(0.03in)
        text(size: 7pt, fill: col-cap)[#caption]
      }
    ],
    // RIGHT: explanation
    [
      #text(size: 23pt, weight: "bold", fill: col-title, hyphenate: false)[#title]
      #v(0.09in)
      #block(
        fill: col-note-bg,
        width: 100%,
        inset: (x: 0.1in, y: 0.09in),
        radius: 3pt,
      )[#set par(leading: 0.65em); #text(size: 10.5pt)[#notes]]
      #v(0.08in)
      #why-box(why)
    ],
  )
}

// --- text-only slide layout (slides 01 and 15) ---
#let text-slide(num, label, title, notes, why) = {
  page-header(num, label)
  v(0.12in)
  align(center)[
    #text(size: 30pt, weight: "bold", fill: col-title)[#title]
  ]
  v(0.22in)
  {[#set par(leading: 0.7em); #text(size: 12pt)[#notes]]}
  v(0.25in)
  why-box(why)
}


// ══════════════════════════════════════════════════════════════════════════════
// PAGE 1 — Slide 01: Title (text-only)
// ══════════════════════════════════════════════════════════════════════════════

#text-slide(
  "01",
  "Title slide — text only",
  "Concerted evolution and unorthodox recombination of human subtelomeres",
  [
    Title slide for a 15-minute BoG 2026 talk (companion to the HPRC v2 Nature submission). Authors: Andrea Guarracino & Erik Garrison.

    The central reframing: what we actually did is survey *inter-chromosomal* relationships between subtelomeres at population scale — 466 near-complete haplotypes — and what we found is that this looks like *concerted evolution*: gene-conversion-like and crossover-like exchange between non-homologous chromosomes, happening more broadly than previously appreciated.

    Frame the talk as the companion to HPRC v2: same data, new question. The HPRC v2 main paper hands us 466 near-complete haplotypes; we use them to re-examine a question genetics parked since the 1990s — how related are the ends of different chromosomes? — and answer it at scale.

    Slide 01 sets up: (a) HPRC v2 framing, (b) the "concerted evolution / ongoing recombination" thesis the closing slide 15 must land. Move quickly to motivation — 30 seconds budget.
  ],
  [This slide sets up the dual thesis (concerted evolution + unorthodox recombination) that slide 15 will close, and anchors the talk as a companion to the HPRC v2 Nature submission — framing every subsequent result as "same data, new question."],
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 2 — Slide 02: Implicit interval tree
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "02",
  "Implicit interval tree — data structure foundation",
  "assets/s02_interval_tree.png",
  "Implicit interval tree — the data structure under the implicit pangenome graph",
  [
    One slide of plumbing before the biology — the rest of the talk is just queries against this object.

    An alignment is an interval on a target sequence: start, end, CIGAR. Li & Rong (2020) showed you can index a set of intervals as an *implicit* tree — a sorted array where tree topology is recovered from the index, not stored as pointers. O(log n) overlap queries, essentially zero memory beyond the intervals themselves.

    We build one tree per target sequence. The union is an interval forest. Walked by transitive closure, that forest *is* the implicit pangenome graph — the same object a pggb graph represents, never materialized. The implementation is IMPG (github.com/pangenome/impg); `impg query -x` does the transitive lookup that powers everything downstream.

    So when a later slide says "we queried the pangenome at chr18 q-arm," what's literally happening is `impg query` walking these trees. Keep this picture in your head — it's the substrate for every plot that follows.
  ],
  [This data structure, one per sequence, is the query interface for the entire pangenome — no GFA is ever built; the forest of implicit trees *is* the graph.],
  caption: "Source: v1 deck slide 2 | Li & Rong 2020 cgranges / IMPG",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 3 — Slide 03: IMPG workflow + Erdős–Rényi callout
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "03",
  "Methods anchor — IMPG workflow + Erdős–Rényi connectivity",
  "assets/s03_impg_workflow.png",
  [The implicit pangenome graph — wfmash all-vs-all *is* the graph (no GFA, no construction step)],
  [
    This is the methods anchor — slow down here. We run wfmash all-vs-all on *n = 18,827* telomere-anchored 500 kb flanks; the resulting PAF set *is* the graph. Each PAF edge is one pairwise alignment. We index it as an interval forest — one implicit interval tree per sequence, in the Li-and-Rong 2020 sense from slide 02. To query we use IMPG, `impg query -x`, which walks transitive closure: start anywhere, chase overlapping mappings outward, collect what's reachable. No graph constructed, nothing to break.

    Is 12% sampling enough? The Erdős-Rényi connectivity threshold for n = 18,827 is p\* = log(n)/n ≈ 5.2×10⁻⁴. Twelve percent is *230× above that* — the random subgraph is densely connected with high probability — closure from any subtelomere reaches every other one. *That* is what licenses "no chromosomal partitioning." Everything downstream rides on this.

    #v(0.05in)
    #box(stroke: 0.5pt + col-cap, radius: 3pt, inset: (x:0.07in, y:0.05in))[
      #image("assets/s03_er_callout.png", width: 100%, height: 1.3in, fit: "contain")
      #text(size: 7pt, fill: col-cap)[ER callout: wfmash 12% sampling is ~230× above the p\* connectivity threshold]
    ]
  ],
  [The Erdős-Rényi argument is the methodological keystone: 230× above threshold means any `impg query -x` walk reaches every subtelomere in the genome — no chromosomal partitioning needed.],
  caption: "Source: v1 deck slide 3 | IMPG: github.com/pangenome/impg",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 4 — Slide 04: Genome-wide identity heatmap
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "04",
  "First result — PAR2-scale pseudohomology at every chromosome end",
  "assets/s04_fig1.png",
  "Genome-wide identity heatmap — interchromosomal homology at PAR2 scale",
  [
    This is the empirical foundation of the talk. For every 100 kb window across 466 HPRCv2 haplotypes we plot the maximum alignment identity to any *other* chromosome. Most of each chromosome is silent — alignments only hit the same chromosome, as expected. But at the telomeres, dense red bands appear wherever assemblies reach the chromosome end.

    These are the inter-chromosomal exchange blocks. The chr18 q-arm inset zooms in on one: a tight band at the very end where many other chromosomes match at over 98% identity, extending tens to hundreds of kilobases inward from the telomere.

    Here is the reframe — this scale is what PAR2 looks like on the sex chromosomes. PAR2 is about *334 kb*. We are seeing *PAR2-scale pseudohomology at nearly every chromosome end*. That dramatically expands the known scope of pseudohomologous regions in the human genome, which is the thesis the next slides will quantify.

    The natural next question: *how many* chromosomes are mixing here? — which slide 05 quantifies.
  ],
  [Fig 1a (manuscript-quality): the chr18q inset shows PAR2-scale pseudohomology replicated genome-wide — this single figure reframes the entire subtelomeric landscape from "X/Y curiosity" to a general human genome phenomenon.],
  caption: "paper_prep/figures/fig1/figure_fig1.png — panel (a) genome-wide + chr18q inset",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 5 — Slide 05: Interchromosomal similarities (n-chromosomes per region)
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "05",
  "Quantifying the landscape — 15,668 PHRs across 41 arms",
  "assets/s05_interchrom.png",
  "Interchromosomal similarities — n-chromosomes per region (HPRCv2)",
  [
    Same view as v1 slide 5 — number of unique chromosomes per 100 kb window across CHM13, with CEN/PAR/PHR/XTR painted as background. Last time: notice the orange spikes at chromosome ends. This time: put numbers on it.

    Across 466 HPRCv2 haplotypes we recovered *15,668 pseudohomologous regions (PHRs)* spanning *41 of 48* chromosome arms. Median length *105 kb*, mean *144 kb*, range 5–500 kb. To anchor that scale: PAR2 is ~334 kb. So a typical subtelomeric PHR is on the same order as PAR2 — and unlike PAR2, this is happening at nearly *every* chromosome end.

    The takeaway is one sentence: *extended pseudohomology at nearly all subtelomeres, comparable in scale to canonical pseudoautosomal regions, but replicated dozens of times across the genome.* That is the central observation the rest of the talk explains — what these communities look like, who is in them, and why they exist.

    These numbers (15,668 / 41 of 48 / median 105 kb / mean 144 kb) anchor every downstream slide.
  ],
  [The PHR count and scale (15,668 PHRs, median 105 kb ≈ PAR2 size, 41/48 arms) is the central quantification the rest of the talk builds on — named clades, 3D contacts, and pedigree exchanges all refer back to this object.],
  caption: "p_num_chromosomes_wide.pdf — v1 slide 5 figure, preserved verbatim",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 6 — Slide 06: Length distributions per arm — outliers are named clades
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "06",
  "Per-arm length distributions — shape encodes the community structure",
  "assets/s06_length_dist.png",
  "Length distributions of inter-chromosomal matches per arm — the outliers are named clades",
  [
    Faceted histograms — one panel per chromosome arm; p-arm = blue, q-arm = orange; pink fill marks the six introvert arms with no inter-chromosomal hits (2p, 3p, 5p, 8q, 11q, 14q). Most arms cluster around the population scale: median ~105 kb, mean 144 kb, range 5–500 kb.

    What to notice: the *shape* of the outlier facets — the ones with the heaviest right tails. Those are not noise, those are the abstract's named clades. *Acrocentric short arms* (13p, 14p, 15p, 21p, 22p — C7): fully homogenized, rDNA-adjacent. *Xq/Yq (C14)*: PAR2, ~334 kb scale. *Xp/Yp (C15)*: PAR1. *4q/10q (C1)*: DUX4/D4Z4, long-tail and copy-number diverse. *10p/18p (C2)*: Linardopoulou 2005 pair.

    Pink panels (2p, 3p, 5p, 8q, 11q, 14q) are a biological signal, not missing data — those arms simply do not participate in the inter-chromosomal exchange landscape.

    *One-sentence reframe: the shape of the per-arm length distribution already recapitulates the community partition you will see on slide 09.*
  ],
  [The fat-right-tail outlier facets are exactly the abstract's six named clades (C1 DUX4, C2 10p–18p, C7 acrocentric, C14 PAR2, C15 PAR1) — the histograms "know" the clades before any clustering is applied.],
  caption: [/moosefs/guarracino/HPRCv2/PHR_III/plots/all-vs-all.1Mb.p95.id95.len_length_dist_by_chr_arm.pdf — v1 verbatim
#h(0.1in)#image("assets/s06_clade_callouts.png", width: 60%, height: 0.6in, fit: "contain")],
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 7 — Slide 07a: All-vs-all heatmap (composite slide, page 1 of 2)
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "07a",
  "All-vs-all heatmap (composite 1/2) — Leiden k=15 + UPGMA dendrogram",
  "assets/s04_fig1.png",
  "All-vs-all heatmap — 41×41 arm-level Jaccard distance with 15 Leiden community boxes",
  [
    This is the all-vs-all picture at the arm level. The matrix is 41×41 — every chromosome arm we have a signal on, against every other arm. Cell color is Jaccard distance on the pangenome graph; the *cyan boxes* are the 15 Leiden communities; the dendrogram on top is the original UPGMA which pulled out 14 of the 15 blocks.

    Fig 1 panel (c) is shown here (full Fig 1 displayed; panel c occupies the right-lower quadrant). The 41×41 arm-level matrix is computed from the full 15,668-flank × 15,668-flank Jaccard similarity matrix by taking median similarity per arm pair.

    *Same figure, page 2 shows the NJ tree on the same 41×41 matrix* — rooted at the acrocentric MRCA, six abstract clades labeled in bold color with 100% perturbation-bootstrap support. Three algorithms (Leiden k=15, UPGMA k=14, NJ) all recover the same clade partition.

    Missing introvert arms (2p, 3p, 5p, 8q, 11q, 14q) are absent from the 41×41 matrix by construction — their absence is a biological signal.
  ],
  [The heatmap shows community structure at arm resolution; the UPGMA dendrogram independently recovers 14/15 Leiden clusters — three-algorithm agreement is the robustness story that answers "is the clustering a Leiden artefact?"],
  caption: "paper_prep/figures/fig1/figure_fig1.png — panel (c) arm-level heatmap (see 07b for NJ tree companion)",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 8 — Slide 07b: NJ tree with named clades (composite slide, page 2 of 2)
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "07b",
  "NJ tree with named clades (composite 2/2) — six abstract clades, 100% bootstrap",
  "assets/s07b_nj_tree.png",
  "NJ tree — same 41×41 Jaccard matrix, rooted at acrocentric MRCA, six abstract clades labeled",
  [
    Neighbor-joining tree on the *same* 41×41 arm-level Jaccard distance matrix. Rooted at the acrocentric short-arm clade (monophyletic, stable orientation). Every clade name in the abstract — PAR1, PAR2, acrocentric short arms, 10p–18p, the tight q-arm clade of 22q/21q/19q/1q/13q/17q, and 4q–10q DUX4 — is a *monophyletic block* on this tree, in bold color, with *100% bootstrap support* at the named-clade root.

    Bootstrap: 1,000 perturbation reps, Gaussian noise σ = 25% of off-diagonal IQR. Named-clade MRCAs never broke. Deeper backbone edges show 32–90% support — the *named clades* are robust; relative ordering is not.

    *NJ ↔ Leiden 1:1 mapping:* PAR1 = C15, PAR2 = C14, acrocentric_p = C7, 10p–18p = C2, tight q-arm = C6, DUX4 = C1. Three algorithms (Leiden/UPGMA/NJ), same answer — the clades are real.

    *This is the highest-value v1→v2 swap*: v1 had only the heatmap with an unlabeled dendrogram; v2 adds the NJ tree the abstract explicitly names.
  ],
  [The NJ tree closes the critical gap from v1: every clade the abstract names is now monophyletically defined with 100% bootstrap support, providing the visual vocabulary for slides 08–15.],
  caption: "paper_prep/figures/nj_tree_arms/nj_tree_annotated.png — upstream nj-tree-from task, commit 602a9d3",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 9 — Slide 08a: PCA/MDS colored by chromosome (composite 1/2)
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "08a",
  "MDS by chromosome (composite 1/2) — arm-community structure is visible",
  "assets/s08a_mds_chrom.png",
  [All-vs-all in 2D — colored by chromosome: arm-community clusters are visible],
  [
    Classical MDS (cmdscale on 1 − Jaccard, k = 5) of the 15,668 × 15,668 Jaccard distance matrix from pggb + odgi similarity. Each point is one subtelomeric flank; color encodes source chromosome (chr1–22, X, Y); shape encodes p/q arm.

    *Left panel result:* points group by Leiden community (15 communities over 41 arms). Visible clusters match arm-community structure — D4Z4 chr4_q + chr10_q, acrocentric p-arms (13/14/15/21/22p), PAR1 (Xp+Yp), PAR2 (Xq+Yq), the 6-arm tight q-arm clade (1q/13q/17q/19q/21q/22q), etc.

    This is the 15,668-flank unfolding of the 41-arm matrix from slide 07. PC1 = 16.05%, PC2 = 11.2% — same axes as v1 slide 10. Missing introvert arms (2p, 3p, 5p, 8q, 11q, 14q) absent throughout.

    *See page 08b for the same projection colored by superpopulation* — the key result: population structure is real but secondary to arm-community structure.
  ],
  [The arm-community clustering in MDS space shows that the 15 Leiden communities are real geometric clusters in the 15,668-flank similarity space, not artefacts of the 41×41 aggregation.],
  caption: "hprcv2.1Mb.subtelo.mds.color-by-chromosome.png — Andrea's pipeline, cmdscale on Jaccard",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 10 — Slide 08b: PCA/MDS colored by superpopulation (composite 2/2)
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "08b",
  "MDS by superpopulation (composite 2/2) — population is secondary to arm-community",
  "assets/s08b_mds_superpop.png",
  [All-vs-all in 2D — colored by superpopulation: clusters do not split by ancestry],
  [
    *Same MDS projection, same points* — only the coloring changes. Recolored by 1KGP superpopulation (AFR/AMR/EAS/EUR/SAS).

    *Right panel result:* the arm-community clusters do NOT split by AFR/AMR/EAS/EUR/SAS — population structure is *real but secondary*. Hudson Fst on cross-arm affinity (ch. 04): mean Fst = 0.044; AFR vs non-AFR pairs 0.10–0.15, non-AFR pairs −0.05 to 0.02. The AFR-deepest split mirrors the human out-of-Africa tree (novel contribution #19) — a population signal exists, but it lives at finer scale than the dominant arm-community signal driving the global 2D layout.

    *The key message:* what arm your subtelomere comes from matters far more than what population you're from. The inter-chromosomal exchange network is an ancient, pre-AFR-divergence structure that has been maintained across human history.

    Combining slides 08a + 08b: the sequence community partition is not a population structure artefact.
  ],
  [The superpopulation-colored view proves that the arm-community clustering is chromosome-identity-driven, not population-driven — ruling out demographic confounding as an explanation for the community structure.],
  caption: "hprcv2.1Mb.subtelo.mds.color-by-superpopulation.png — same cmdscale projection as 08a",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 11 — Slide 09: PCA communities — keystone slide
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "09",
  "Keystone slide — 15 arm-level communities, abstract clades named",
  "assets/s09_pca_communities.png",
  "All-vs-all PCA — 15 arm-level communities, named for the abstract's clades",
  [
    The keystone slide: each point is one chromosome arm in PHR-similarity space; clusters = Leiden k=15 communities on the 41×41 arm-level Jaccard distance matrix. *Every named clade in the abstract is a community on this plot.*

    - *PAR1 = C15* (Xp/Yp + 18q outlier)
    - *PAR2 = C14* (Xq/Yq)
    - *Acrocentric short arms = C7* (13p,14p,15p,21p,22p)
    - *Tight q-arm clade = C6* — exact 6-arm match: 22q/21q/19q/1q/13q/17q
    - *10p–18p Linardopoulou pair = C2*
    - *4q–10q DUX4 = C1*

    Three interpretive zones: *PAR-driven* (lower-left, C14+C15), *concerted-exchange PHR core* (center, C6+C7+C3), *DUX4/D4Z4* (right, C1). PC1=16.05%, PC2=11.2%. Missing introvert arms (2p,3p,5p,8q,11q,14q): same six arms as slides 07–08.

    After this slide the rest of the talk is *evidence for* what is grouped here.
  ],
  [This is the single slide that ties methods to abstract: every clade word the audience read in the abstract now has a colored cluster, an arm list, and a recombinational interpretation — making slides 10–15 comprehensible as evidence for this map.],
  caption: [v1 deck slide 10 (v1_page_10-10.png) — Leiden k=15 arm-level communities
    #h(0.05in)#image("assets/s09_clade_legend.png", width: 85%, height: 0.8in, fit: "contain")],
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 12 — Slide 10a: Hi-C / Pore-C bulk contact matrix (composite 1/2)
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "10a",
  "Sequence → nucleus: bulk contact matrix confirms 3D community structure (composite 1/2)",
  "assets/s10_fig3.png",
  [Hi-C / Pore-C — sequence communities are physically co-localized in 3D nuclei],
  [
    Moving from sequence to nucleus. *Left panel (Fig 3a):* the HG002 Pore-C inter-arm contact matrix at 50 kb, 77 arm-haplotypes, ordered by Leiden sequence community. The diagonal blocks are precisely the communities we built from pangenome graph similarity — within-community contacts are *18-fold higher* than between, p ≈ 10⁻⁸⁵. Sequence-defined communities are physical.

    *Right panel (bulk Mantel):* across eight HPRC datasets, arms with more similar subtelomeres also contact each other more often. HG002 Hi-C and CHM13 each ρ = 0.66; Pore-C ρ = 0.49; HG02148 ρ = 0.15 (marginal).

    *Robustness check:* strip out chr13–22 p/q-arms, chrX, chrY, and the four strongest communities (D4Z4, acrocentric p, PAR1, PAR2) — every sample's ρ *goes up*: HG002 to 0.79, CHM13 to 0.85, HG02148 to 0.72. The signal is not a nucleolar artifact — it is a generic property of subtelomeric homology.

    *See 10b for the Mantel exclusion figure.*
  ],
  [The Pore-C contact matrix shows that sequence-community membership predicts 3D proximity at p ≈ 10⁻⁸⁵ — the sequence clustering is not arbitrary but reflects genuine nuclear organization.],
  caption: "paper_prep/figures/fig3/figure_fig3.png — Fig 3 (panels a+c shown; see 10b for ED5b Mantel)",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 13 — Slide 10b: Mantel exclusion robustness (composite 2/2)
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "10b",
  "Mantel robustness — signal strengthens when known confounds are removed (composite 2/2)",
  "assets/s10b_ed5.png",
  [Mantel test: ρ before vs after acrocentric + sex-chr + strong-community exclusion],
  [
    *Extended Data Fig 5 panel (b):* the bulk Mantel test comparing subtelomeric Jaccard similarity versus inter-chromosomal contact, for each of 8 HPRC Hi-C/Pore-C datasets, shown *before* and *after* excluding the strongest communities (acrocentric short arms, PAR1, PAR2) plus sex chromosomes.

    The skeptic's worry: the Mantel signal is just acrocentric nucleolar clustering (chr13/14/15/21/22 all co-localize for rDNA replication) or pseudoautosomal contact (X/Y obligate crossover). So strip those out — every sample's ρ *goes up*, not down.

    Key numbers after exclusion: HG002 Hi-C 0.79 (up from 0.66), CHM13 Hi-C 0.85 (up from 0.66), HG02148 0.72 (up from 0.15). Drop to 10 kb resolution and treat individual sequence pairs without community labels → 0.83 in NA19036, 0.81 in HG02148.

    *Conclusion:* the Mantel signal is a *generic* property of subtelomeric sequence homology, not a confound-driven artefact. The more carefully you look, the stronger it gets.
  ],
  [The Mantel exclusion test is the decisive rebuttal to "you're just measuring nucleolar clustering" — removing the obvious drivers (acrocentrics, sex chrs) makes the signal stronger, proving it is a generic subtelomeric property.],
  caption: "paper_prep/figures/ed5/figure_ed5.png — Extended Data Fig 5b Mantel ρ before/after exclusion",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 14 — Slide 11: Single-cell 3D (Dip-C + haploid sperm)
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "11",
  "Single-cell 3D — signal survives per-cell and in haploid sperm",
  "assets/s10_fig3.png",
  "Single-cell 3D — and it works in haploid sperm",
  [
    The bulk Hi-C signal could be a population average artefact — many cells with weak contacts averaging into a coherent block. Single-cell Dip-C rules that out.

    *GM12878, 16 cells (Tan 2018, remapped T2T-CHM13):* community-member arms ~7% closer to each other than to non-members — Wilcoxon p = 3.8×10⁻⁴, Mantel ρ = 0.30. Small per-cell effect, present in essentially every cell.

    *Punch line:* same pipeline on *20 sperm cells (Xu 2025)*. Sperm is haploid, chromatin hyper-condensed, nuclear architecture nothing like interphase. Within/between ratio = 0.40 — community arms *60% closer* in 3D — Fisher p = 3.9×10⁻⁵¹. Same direction, much stronger.

    *Negative control:* a pseudo-community of the 7 arms sharing no subtelomeric sequence moves the *opposite* way — 11% farther in GM12878, 40% farther in sperm. Sequence sharing is *necessary* for clustering; a community label alone is not.

    The 3D signal survives the bulk→single-cell test and the haploid germline test. Sperm is the bridge to meiosis on slide 12.
  ],
  [Sperm single-cell 3D (W/B = 0.40, p = 3.9×10⁻⁵¹) proves the community co-localization is present in the gametic genome — the haploid cell that carries recombination forward still shows the same proximity pattern.],
  caption: "paper_prep/figures/fig3/figure_fig3.png — Fig 3 panel (c) single-cell enrichment (W/B per cell)",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 15 — Slide 12: Mouse meiotic Hi-C — zygotene bouquet peak
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "12",
  "Meiotic Hi-C in mouse — zygotene bouquet is where the 3D signal peaks",
  "assets/s12_fig4.png",
  [Mouse meiotic Hi-C — the zygotene bouquet is where the 3D signal peaks],
  [
    Bulk Hi-C is mitotic — LCLs, RPE-1, PBMCs — but the recombination we are explaining happens in meiosis. Mouse meiotic Hi-C is the answer. Zuo et al. (2021) sorted four prophase stages (leptotene, zygotene, pachytene, diplotene) and produced stage-specific Hi-C on a genome we now have a T2T assembly for.

    Run the same Mantel test against B6 + CAST mouse subtelomeres. ρ across stages: leptotene 0.687, *zygotene 0.718*, pachytene 0.683, diplotene 0.577. *Zygotene is the peak.*

    Per-PHR-pair Spearman at zygotene: ρ = 0.715, p = 4.4×10⁻⁵⁵ across 344 inter-chromosomal pairs.

    Zygotene is the *bouquet* — telomeres clustered at the LINC-anchored nuclear envelope while homologs align. The 3D signal we see in human LCLs is the somatic shadow of a meiotic event we can watch in mouse, at exactly the stage when telomeres are physically together.

    The stage trajectory inset (right) shows the rise → peak → fall pattern across all four stages.

    #image("assets/s12_trajectory.png", width: 100%, height: 0.9in, fit: "contain")
  ],
  [The zygotene peak (ρ = 0.718) directly links the sequence community structure to the meiotic bouquet — the stage when subtelomeres are physically clustered and the recombination that produces ongoing exchange must happen.],
  caption: "paper_prep/figures/fig4/figure_fig4.png — Fig 4 panel (d) mouse meiotic; inset = stage_trajectory.png",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 16 — Slide 13: Pedigree direct evidence
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "13",
  "Direct proof — ongoing exchange in three generations of a T2T pedigree",
  "assets/s13_pedigree.png",
  "Caught in the act — three generations of a T2T pedigree show ongoing subtelomeric exchange",
  [
    Proof that exchange is *ongoing*, not frozen. Substrate: WashU pedigree (Cechova 2025) — PAN010 grandmother, PAN011 grandfather, PAN027 mother, PAN028 daughter, all telomere-to-telomere. Push each child haplotype through `odgi untangle` against its parent in the implicit pangenome graph.

    The figure: PAN027's maternal haplotype painted onto mother PAN010. Diagonal stripes = correctly inherited self-flank. Off-color stacks across the acrocentric short arms = gene-conversion-like patches from other chromosome ends.

    *Numbers:* 538 high-quality inter-chromosomal patches total. *494/538 (92%)* fall inside Leiden communities built from 233 independent HPRC v2 individuals. The graph from the population *predicts* where exchange happens in the family.

    *133* events have textbook gene-conversion-like sandwich pattern. *16* are crossover-like (left+right haplotype flip across an inter-chromosomal patch — real meiotic crossover). Communities: C7 (acrocentric), C2 (10p↔18p), C15 (PAR1), and C1 — *DUX4, chr4q onto chr10q* at score 0.957 in PAN028. Disease-named exchange in a normal family.

    *Replication:* CEPH1463 pedigree, 4 generations, both hifiasm+verkko assemblers — 11 events survive the double-assembler filter.
  ],
  [494/538 (92%) of observed inter-chromosomal exchange patches fall inside sequence-community boundaries predicted by independent population data — the community structure is causally predictive of where meiotic exchange occurs in real families.],
  caption: "PAN027.maternal_hap1_from_PAN010_mother.untangle.pdf (odgi untangle ribbon, WashU pedigree)",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 17 — Slide 14: Gene biology — DUX4, OR4F, TAR1
// ══════════════════════════════════════════════════════════════════════════════

#slide(
  "14",
  "Biology of the exchange network — DUX4, OR4F, TAR1",
  "assets/s14_gene_biology.png",
  "Gene biology aside — DUX4, OR4F, TAR1 (the biology is interesting too)",
  [
    A quick aside before the closer. Three vignettes showing the biology inside these regions — each reads out a different aspect of the exchange network.

    *DUX4:* annotated across 18 q-arms in the pangenome, but the medically relevant biology (D4Z4 macrosatellite contraction → FSHD muscular dystrophy) lives only at chr4q and chr10q — the C1 community. Median *22 DUX4L copies* per haplotype on those two arms; everywhere else 0–2. DUX4 has spread, but FSHD is geometrically constrained to one community of fifteen.

    *OR4F olfactory receptors:* four paralogs, 16 arms, ~5,000 gene-copy entries. Pseudogenization rate sweeps from 11% at chr7p to 99.8% at chr15q. Same gene family across the same exchange network — the decay clock has been running for very different lengths of time at different ends.

    *TAR1 (telomere-associated repeat):* 94.6% of all subtelomeric sequences carry TAR1, across all 41 arms. The one exception: PAR1 (chrX_p + chrY_p) at 0.5%. PAR1 has obligate meiotic crossover — no sequence anchor needed for exchange; TAR1 is absent. This marks places that *need* a sequence anchor for exchange.

    This is compressible to zero if time is short, but DUX4/FSHD is the disease-revealed instance of the talk's central process.
  ],
  [The three gene systems (DUX4 = disease readout, OR4F = decay clock, TAR1 = exchange machinery marker) each write themselves into the community architecture — showing the biological stakes of the sequence-sharing network beyond genomics.],
  caption: "slides/v2/slide_14_gene_biology.pdf (rendered from slide_14_gene_biology.R — three-panel: DUX4/OR4F/TAR1)",
)

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 18 — Slide 15: Closer (part A) — pillar summary table
// ══════════════════════════════════════════════════════════════════════════════

#{
  page-header("15", "Thesis closer — pillar summary")
  v(0.12in)
  align(center)[
    #text(size: 22pt, weight: "bold", fill: col-title)[Concerted evolution of human subtelomeres — what we saw, predicted, and recovered]
  ]
  v(0.22in)
  align(center)[
    #block(width: 7.5in)[
      #set text(size: 10pt)
      #set par(leading: 0.65em)
      #table(
        columns: (1.75in, 5.75in),
        stroke: none,
        inset: (x: 0.09in, y: 0.12in),
        [*Method (slide 3):*],
        [Implicit pangenome graph: wfmash all-vs-all over 18,827 telomere-anchored flanks, ~12% of all pairs evaluated — *230× above the Erdős-Rényi connectivity threshold* (p\* = log(n)/n ≈ 5.21×10⁻⁴). No chromosome partitioning, no GFA.],
        [*Empirical (slides 4–9):*],
        [15,668 PHRs across 41/48 arms, median 105 kb / mean 144 kb — PAR2-scale pseudohomology at nearly every chromosome end. Named clades: Xp/Yp & Xq/Yq via PARs, *acrocentric short arms* (C7), *10p–18p*, the big q-arm clade (*22q–21q–19q–1q–13q–17q*), and *4q–10q DUX4*.],
        [*Mechanism (slides 10–12):*],
        [Hi-C/Pore-C/Dip-C/sperm all recover community-structured 3D contacts (Mantel ρ = 0.296, p = 0.002; per-pair ρ = 0.674 in CHM13 Hi-C). Median PHR (105 kb) sits at the base of a single meiotic loop — the *bouquet* is the predicted exchange venue.],
        [*Proof (slide 13):*],
        [WashU T2T pedigree: *494/538 (92%)* inter-chromosomal patches fall inside Leiden communities — 133 gene-conversion-like, 16 crossover-like.],
        [*Biology (slide 14):*],
        [D4Z4/DUX4 (4q↔10q) is the disease-revealed instance of the same process.],
      )
    ]
  ]
}

#pagebreak()

// ══════════════════════════════════════════════════════════════════════════════
// PAGE 18b — Slide 15: Closer (part B) — thesis statement
// ══════════════════════════════════════════════════════════════════════════════

#{
  page-header("15", "Thesis closer — synthesis")
  v(1.4in)
  align(center)[
    #block(width: 9.5in)[
      #block(
        fill: col-why-bg,
        inset: (x: 0.35in, y: 0.32in),
        radius: 6pt,
        stroke: (left: 6pt + col-why-bar),
        width: 100%,
      )[
        #text(size: 17pt, weight: "bold")[Thesis: ]#text(size: 17pt)[_Subtelomeres concertedly evolve through ongoing inter-chromosomal exchange — observable in pedigrees, predicted by 3D, recovered by an implicit pangenome graph across HPRC v2._]
      ]
    ]
  ]
  v(0.5in)
  why-box([The closing slide recapitulates all five pillars (method → empirical → mechanism → proof → biology) in 70 seconds and ends on the verbatim thesis sentence — the title's "concerted evolution" echoed from slide 01 to slide 15.])
}
