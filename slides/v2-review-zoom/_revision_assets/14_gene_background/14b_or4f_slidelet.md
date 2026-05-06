# Slidelet 14b spec: OR4F

Status: candidate slidelet spec, not a rendered deck asset.

## Review Goal

Keep the useful OR4F pseudogenization gradient, but connect it explicitly to the
gene-enrichment outputs and the broader subtelomeric duplicon story.

## One-Line Takeaway

OR4F olfactory receptor genes are widespread subtelomeric duplicon markers; the
PHR analysis turns that known subtelomeric gene family into a population-scale
per-arm pseudogenization gradient.

## On-Slide Copy

Title:

OR4F: the olfactory-receptor decay gradient

Main visual label:

OR4F pseudogene fraction by chromosome arm

Callouts:

- 5,023 OR4F annotations across 16 arms.
- 11.1% pseudogene at chr7p.
- 99.8% pseudogene at chr15q.
- Population mean: 62.1%.
- OR4F genes are detected across 7 communities; OR4F5 and OR4F8P each appear on 14 arms.

Footer:

Why this matters: OR4F is not the main PHR mechanism; it is a readable marker of
how shared subtelomeric gene families decay differently on different chromosome
ends.

## Candidate Layout

Preserve the current bar-gradient concept:

1. Full-width sorted bar plot from low to high pseudogene fraction.
2. Use one neutral color for most arms, red only for chr7p and chr15q.
3. Put the mean as a thin dashed horizontal line with a small label.
4. Add a narrow right-side evidence box:
   "10 OR4F genes, 7 communities; OR4F5/OR4F8P on 14 arms".
5. Keep the text beside the plot, not on top of it.

Do not overstate Fisher enrichment. The local gene-enrichment survey says the
OR4F pattern is a qualitative presence pattern after multiple testing, while
copy-weighted ORA is a separate synthesis choice. The slide can mention
"copy-rich olfactory signal" in notes, but the on-slide claim should stay with
observed distribution and pseudogenization.

## Speaker Notes

OR4F is useful because it is easy to understand and visually strong. These are
olfactory receptor family members that live in subtelomeric duplicon sequence.
In our annotations, OR4F appears broadly: 10 family members across 7
communities, with OR4F5 and OR4F8P each present on 14 arms. The striking thing
for the slide is not only presence, but decay. Across 5,023 annotations, the
pseudogene fraction spans from 11.1% at chr7p to 99.8% at chr15q. Same gene
family, similar subtelomeric neighborhood class, very different decay history.

Suggested spoken version:

"OR4F gives us a cleaner gene-family readout. It is an olfactory receptor
family spread across subtelomeric duplicons. The slide keeps the gradient
because it works: from 11 percent pseudogene on chr7p to essentially all
pseudogene on chr15q. That is PHR biology too: shared gene content is not
static; it ages differently on different chromosome ends."

## Source Anchors

- `slides/v2/slide_14_gene_biology.md:8` gives the existing OR4F slide claim:
  four paralogs, 16 arms, 5,023 entries, 11.1% to 99.8% gradient, 62.1% mean.
- `slides/v2/slide_14_gene_biology.R:69` through
  `slides/v2/slide_14_gene_biology.R:98` implements the current OR4F gradient.
- `paper_prep/figures/ed4/caption.md:7` says high-copy OR4F variants are part
  of the top high-copy coding/pseudogene families.
- `paper_prep/figures/ed4/caption.md:9` gives the manuscript OR4F gradient:
  5,023 entries across 16 arms, 11.1% to 99.8%, mean 62.1%.
- `paper_prep/figures/ed4/sources.tsv:6` gives the off-tree source table:
  `/moosefs/guarracino/HPRCv2/PHR_III/plots/or4f_pseudogene_fraction.csv`.
- `subtelomeric_analysis_report.md:511` through
  `subtelomeric_analysis_report.md:515` connects OR4F to gene enrichment:
  10 OR4F family genes across 7 communities; OR4F5 and OR4F8P on 14 arms.
- `subtelomeric_analysis_report.md:526` maps Ambrosini block 2 at chr15q to
  OR4F17/OR4F4 in C8.
- `paper_prep/surveys/SURVEY_03_gene_enrichment.md:21` through
  `paper_prep/surveys/SURVEY_03_gene_enrichment.md:29` gives the enrichment
  caveat: qualitative presence patterns do not survive BH correction.
- `paper_prep/surveys/SURVEY_03_gene_enrichment.md:141` through
  `paper_prep/surveys/SURVEY_03_gene_enrichment.md:144` records the
  copy-weighted ORA exploration, including the 598-fold olfactory signal.
- `paper_prep/surveys/SURVEY_10_11_12_limits_summary_lit.md:343` through
  `paper_prep/surveys/SURVEY_10_11_12_limits_summary_lit.md:344` locks the
  D4Z4 and OR4F quantitative sub-confirmations.
- `subtelomeric_analysis_report.md:2386` lists Linardopoulou 2005 as the local
  literature anchor for OR4F21/block 5 subtelomeric sharing.

## Provenance Notes

The current visible `s14_or4f.png` is an exact-copy crop from agent-878. The
underlying plot concept is still useful and should be preserved, but any new
visual should be regenerated from `or4f_pseudogene_fraction.csv` or documented
as a deliberate crop. See `README.md` in this directory for the full git
timeline and blob hashes.
