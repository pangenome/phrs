# SURVEY 02 — Annotation section

Source: `end-to-end-report/report/02_annotation.md`
Scope: gene + repeat (TAR1) annotation pipeline, internal (TTAGGG)n islands, terminal telomere length per community.

---

## 1. Key findings with metrics

**Annotation coverage**
- Inputs: 464 haplotype-specific Liftoff GFF3 + CHM13 (465 total). HG002 annotated separately with JHU Liftoff v0.6 on HG002v1.1.
- 18,827 PHR sequences → **15,668 retained** (3,158 removed for no inter-chromosomal signal; chimeric chr18_q already excluded upstream).
- **173,881 gene annotations** (374 unique genes) across 232 individuals and 39 arms; chr7_q and chr12_q yield zero genes (PHR confined to 5–25 kb tip, below Liftoff threshold).
- **21,544 TAR1 entries** across 14,816 sequences (94.6%) and all 41 arms.

**TAR1 prevalence**
- Present in all 41 arms; >99% on most autosomal arms.
- Near-absent on PAR1: chrX_p 1/327 (0.3%), chrY_p 1/92 (1.1%) — consistent with obligate meiotic crossover, not repeat-mediated exchange (Rouyer et al. 1986).
- Acrocentric p-arms intermediate (73.1–78.9%); highest density chr18_p (4.00 copies/seq) and chr9_q (2.82).

**TAR1 positional distribution**
- 66.9% of 21,544 entries within 10 kb of telomere; 70.3% within 25 kb (median 0.3 kb, mean 44.1 kb).
- Most arms median <1 kb from telomere; deeper TAR1 in acrocentric p-arms (chr13_p–chr22_p median 179–196 kb), chr9_q (223 kb), chr8_p (140 kb), chr18_p (87 kb), chr19_p (69 kb), chr6_p (57 kb).

**TAR1 per community (table in source)**
- C15 (PAR1) effectively TAR1-free (0.5%).
- C7 (acrocentric p) lowest non-PAR (76.9%).
- C2 (chr10_p/chr18_p) highest density (mean 2.51 copies/seq).

**Internal (TTAGGG)n islands**
- **18,352 islands across 8,321 sequences (53.1%)**; all 41 arms. Median 79 bp, mean 102 bp; max 22 islands per sequence.
- Top arms by count: chr20_q (1,765), chr12_q (1,149), chr16_p (898), chr18_p (851).
- 16/21 p-arms have median island position <500 bp from telomere; 5 p-arms with deep PHR show deeper islands (chr11_p 153 kb, chr6_p 99 kb, chr20_p 65 kb, chr18_p 54 kb, chr16_p 1.1 kb).

**TTAGGG island boundary enrichment test**
- KS stat=0.37, p<1e-300; mean fractional position 0.54; 42.5% within 5 kb of PHR outer boundary (7,806/18,352).
- Caveat: tested PHR outer boundary, NOT internal duplicon-to-duplicon boundaries (Ambrosini et al. 2007 claim untested).

**TTAGGG island length distribution**
- Mode 50–74 bp (46.0%), monotonically decreasing; only 7.9% in the 150–199 bp range Ambrosini et al. (2007) reported as their mode. Difference is methodological (50 bp floor + degenerate motif search vs single-reference detection).

**TTAGGG island motif composition**
- 296,406 hexamer instances. Canonical TTAGGG 52.3%; variants TGAGGG 19.0%, TTGGGG 16.0%, TCAGGG 12.7%.
- Only 32.2% of islands "pure canonical" (≥80% TTAGGG+CCCTAA); 47.2% variant-dominant (<50% canonical).
- Same three dominant variants as Ambrosini et al. (2007). Linardopoulou et al. (2005): degenerate telomeric repeats enriched at 4% of subtelomeric DSB sites vs 0.5% background.

**Island count by exchange status**
- 8,321 sequences with islands: 1,569 cross-arm (18.9%), 6,752 self-arm (81.1%). Mean 2.08 vs 2.24. Mann-Whitney U z=−1.89, **p=0.045** — marginal, small effect; conclusion: no meaningful effect of cross-arm status on island count.

**Terminal telomere tract length by community**
- Kruskal-Wallis H=100.89, **p=3.2e-15** across communities. Medians 7,638 bp (C10/chr17_p) → 9,418 bp (C13/chr4_p). Range 470–33,826 bp.
- Longest mean: C13 chr4_p (9,640 bp), C15 PAR1 (9,266 bp). Shortest: C9 chr16_p (8,360), C10 chr17_p (8,178).
- Telomere length × island count (n=8,321 with island): Spearman ρ=−0.056, p=2.7e-7 (significant but trivial effect size).

---

## 2. Existing figures/plots referenced (paths)

The `02_annotation.md` source **references no figure files**. The only annotation-adjacent plots already on disk live in the heterogeneity directory and re-use TAR1 data:

- `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/within_arm_heterogeneity_11_tar1_prevalence.{pdf,png}` — TAR1 prevalence (used by section 04).
- `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/within_arm_heterogeneity_12_tar1_cross_vs_self.{pdf,png}` — TAR1 cross-arm vs self-arm comparison (used by section 04).

No figures exist for: TAR1 positional distribution, internal (TTAGGG)n islands (length, motif composition, position, count by exchange status), or terminal telomere length per community.

---

## 3. Existing CSVs / data files (paths)

**Inputs**
- `/moosefs/guarracino/HPRCv2/PHR_III/hprc_annotations/*.gff3.gz` — 464 haplotype-specific Liftoff GFF3 (HPRC index covers 462; HG002 added separately).
- `/moosefs/guarracino/HPRCv2/PHR_III/hprc_repeatmasker/*.RepeatMasker.bed.gz` — 465 haplotype-specific repeat annotations.
- `/moosefs/guarracino/HPRCv2/PHR_III/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz` (+ `.fai`, `.gzi`) — sequences scanned for internal islands and (presumably) terminal telomere TSV.

**Annotation intersection product**
- `/moosefs/guarracino/HPRCv2/PHR_III/annotations/subtelomeric_annotations.1Mb.rds` — gene + TAR1 annotations intersected with PHR regions.

**TAR1 outputs**
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_summary.tsv` — TAR1 prevalence per community (drives §1.3 table).
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_by_arm.tsv` — per-arm/community TAR1 counts.
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/tar1_positional/tar1_positional_overall.tsv` — distance-from-telomere distribution.
- `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/tar1_positional/tar1_positional_per_arm.tsv` — per-arm TAR1 positional summary.
- `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/cross_arm_tar1_comparison.tsv` — TAR1 in cross-arm vs self-arm sequences.

**TTAGGG island outputs**
- `/moosefs/guarracino/HPRCv2/PHR_III/ttaggg_analysis/boundary_enrichment.tsv` — KS / binomial boundary test inputs.
- `/moosefs/guarracino/HPRCv2/PHR_III/ttaggg_analysis/length_distribution.tsv` — island length histogram bins.
- `/moosefs/guarracino/HPRCv2/PHR_III/ttaggg_analysis/motif_composition.tsv` — per-island canonical/variant fractions.
- `/moosefs/guarracino/HPRCv2/PHR_III/ttaggg_analysis/per_island_hexamers.tsv` — raw hexamer-instance table.
- `/moosefs/guarracino/HPRCv2/PHR_III/ttaggg_analysis/cross_arm_island_count.tsv` — island count by cross-arm status.
- `/moosefs/guarracino/HPRCv2/PHR_III/heterogeneity/island_exchange_status.tsv` — per-sequence island count + exchange status (joined).

**Terminal telomere data** — referenced as `.telo.tsv` in source (per-sequence terminal telomere lengths). Likely upstream of `hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz` trim step; source path not stated in this section.

---

## 4. Methods used

- **Annotation pipeline**: HPRC remote Liftoff GFF3 index for 462 haplotypes + JHU Liftoff v0.6 on HG002v1.1 for the remaining 2 → CHM13 added → 465 GFF3. RepeatMasker BED extracted to a TAR1-only set; CHM13 converted to PanSN naming; HG002 converted bigBed → 10-column BED. Intersection of GFF3/TAR1 with PHR regions written to a single `.rds`.
- **TAR1 positional**: per-entry distance to the telomeric end of the host PHR region; per-arm summary statistics (median, mean).
- **Internal (TTAGGG)n island detection**: `seqkit locate` for canonical TTAGGG plus three variants (TGAGGG, TCAGGG, TTGGGG) on all 15,668 PHR sequences; overlapping hits merged with 12 bp gap tolerance; filtered to 50–1000 bp.
- **Boundary enrichment**: KS test for non-uniformity of fractional position within PHR; binomial test for boundary-proximal counts (5 kb window).
- **Motif composition**: hexamer counting per island; canonical fraction = (TTAGGG + CCCTAA) / total hexamers.
- **Cross-arm vs self-arm**: Mann-Whitney U on island counts.
- **Telomere length × community**: Kruskal-Wallis across community medians; Spearman correlation between terminal telomere length and per-sequence island count (restricted to n=8,321 with ≥1 island).
- **Statistical tools and significance thresholds** are not explicitly stated in the section (assume R/Python defaults).

**Scripts (cited in source)**
- `/moosefs/guarracino/HPRCv2/scripts/preprocessing/preprocess-subtelomeric-annotations.R` — annotation intersection.
- `/moosefs/guarracino/HPRCv2/scripts/community/analyze-tar1-positional.R` — TAR1 positional.
- `/moosefs/guarracino/HPRCv2/scripts/community/analyze-ttaggg-islands.py` — island boundary, length, cross-arm count.
- `/moosefs/guarracino/HPRCv2/scripts/community/analyze-ttaggg-motifs.py` — hexamer composition.
- `/moosefs/guarracino/HPRCv2/scripts/community/analyze-island-exchange-status.py` — exchange-status linkage.
- `/moosefs/guarracino/HPRCv2/scripts/community/ttaggg_boundary_enrichment.py` — KS / binomial test.
- `/moosefs/guarracino/HPRCv2/scripts/community/telomere_length_by_community.py` — Kruskal-Wallis on telomere length.

---

## 5. Open gaps for figure (re)generation

1. **No telomere TSV path documented.** The Kruskal-Wallis result depends on a `.telo.tsv` file the section does not locate. To regenerate any telomere-length figure we need the canonical path (e.g. RUKKI/seqtk-telo output upstream of `hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz`).
2. **Liftoff version inconsistency.** Source mentions JHU Liftoff v0.6 for HG002 only; the upstream HPRC annotation index Liftoff version is not stated. Methods will need to disclose both.
3. **TTAGGG search regex/code** is partially described (motif list + 50 bp floor + 12 bp merge gap) but the exact `seqkit locate` invocation, strand handling, and filtering of overlaps with the **terminal** telomere array vs internal islands are not in the section. A reproducible figure caption needs that.
4. **Boundary test scope mismatch.** The KS/binomial test addresses PHR outer boundary; Ambrosini et al.'s claim concerns *internal* duplicon boundaries. A direct test requires per-PHR duplicon-block coordinates, which are not produced in this section. Either generate them or restate the claim narrowly.
5. **Acrocentric p-arm "intermediate" prevalence (~73–79%)** is reported as a single statistic per arm without per-haplotype variability or rDNA-array adjacency context. A density-vs-distance-to-rDNA panel may be warranted.
6. **Telomere ↔ island Spearman ρ=−0.056** is statistically significant but trivially small; conclusion language ("biological significance unclear") should be paired with a scatter or bin-summary figure to communicate the effect-size verdict honestly.
7. **TTAGGG length distribution comparison with Ambrosini et al.** is asserted methodologically; a side-by-side density plot (or stratified histogram by minimum-length filter) would let the reader verify the 50 bp threshold accounts for the mode shift.
8. **Variant motif fractions** (TGAGGG/TTGGGG/TCAGGG) are reported population-aggregated; per-arm or per-community breakdown is not provided and may be informative.
9. **Per-sequence island count distribution** (currently summarised as median 1–4, max 22) lacks a histogram artifact.

---

## 6. Suggested figures

**Convention.** [PROD] = already produced on disk; [GEN] = to generate; main vs Extended/SI judged for a Nature manuscript focused on subtelomeric exchange, where annotation supports but is not the central story.

### Main-text candidates
- **Fig M1. TAR1 + internal-island landscape across 41 chromosome arms** [GEN]
  Multi-panel arm-ordered plot. (a) TAR1 prevalence per arm (highlight PAR1 absence, acrocentric drop, autosomal saturation); (b) TAR1 positional density (telomere-proximal mode + acrocentric/chr9_q deep tail); (c) Internal (TTAGGG)n island count per arm; (d) Median island distance from telomere per arm.
  *Caption.* "Subtelomeric repeat content is structured along chromosome arms. (a) TAR1 prevalence is near-saturating on autosomes (>99%), absent in PAR1, and intermediate (73–79%) on acrocentric p-arms. (b) TAR1 sits within 10 kb of the telomere on most arms; acrocentric p-arms and chr9_q carry deeper TAR1 reflecting their large PHR regions. (c) Internal (TTAGGG)n islands are most abundant on chr20_q, chr12_q, chr16_p, chr18_p. (d) Most p-arms place islands near the telomere; arms with large PHRs (chr11_p, chr6_p, chr20_p, chr18_p, chr16_p) extend islands deeper. n=15,668 PHR sequences across 465 near-complete assemblies."

- **Fig M2. TTAGGG island composition reveals degenerate telomeric remnants** [GEN]
  (a) Length distribution histogram with mode at 50–74 bp; overlay Ambrosini et al. (2007) reported 150–200 bp range; (b) Hexamer-composition stacked bar (TTAGGG / TGAGGG / TTGGGG / TCAGGG); (c) Distribution of canonical-fraction per island (showing only 32.2% are ≥80% canonical).
  *Caption.* "Internal (TTAGGG)n islands are short and motif-degenerate. (a) Population-scale island lengths (n=18,352) are dominated by 50–74 bp tracts (46.0%); only 7.9% lie in the 150–200 bp range that Ambrosini et al. (2007) called the mode. (b) Hexamer composition is 52.3% canonical and 47.7% variant, with the same three variant motifs (TGAGGG, TTGGGG, TCAGGG) Ambrosini et al. identified. (c) 47.2% of islands are variant-dominant (<50% canonical), consistent with degenerate remnants of duplicated telomeric sequence."

### Extended Data candidates
- **ED1. TAR1 prevalence by community** [GEN] — bar chart of the §1.3 table (C15 PAR1 0.5% → C2 chr10_p/chr18_p 99.8%, mean 2.51 copies); reuse `/moosefs/guarracino/HPRCv2/PHR_III/enrichment/community_tar1_summary.tsv`. (Note: section 04 has a related plot at `…/heterogeneity/within_arm_heterogeneity_11_tar1_prevalence.{pdf,png}` [PROD]; verify whether to re-cut for annotation framing or cross-reference.)
- **ED2. PHR-boundary enrichment of internal islands** [GEN] — KS-style ECDF of fractional island position within PHR, with 5 kb boundary window highlighted (KS=0.37, p<1e-300; 42.5% within 5 kb). Show why the conclusion is "non-uniform but not strongly polarized."
- **ED3. Terminal telomere length by community** [GEN] — boxplot of `.telo.tsv` lengths grouped by community, ordered by median (C10 7,638 → C13 9,418 bp); annotate Kruskal-Wallis H=100.89, p=3.2e-15.

### Supplementary candidates
- **SI1. Island count by cross-arm vs self-arm exchange status** [GEN] — violin/box of per-sequence island counts (cross 2.08 vs self 2.24, U test z=−1.89, p=0.045); explicitly call out small effect size.
- **SI2. Telomere-length × island-count scatter** [GEN] — restricted to n=8,321 with ≥1 island; Spearman ρ=−0.056. Underline near-null effect.
- **SI3. Liftoff annotation coverage matrix** [GEN] — heatmap (samples × arms) of gene annotation count, flagging chr7_q and chr12_q empties.
- **SI4. TAR1 cross-arm vs self-arm** [PROD partial] — already at `…/heterogeneity/within_arm_heterogeneity_12_tar1_cross_vs_self.{pdf,png}`; cross-reference rather than regenerate.
- **SI5. Per-arm TTAGGG island position density (small-multiples)** [GEN] — 41-panel facet using `tar1_positional_per_arm.tsv`-style inputs derived from `per_island_hexamers.tsv`, useful to support the "deep on chr11_p / chr6_p" claim.
- **SI6. Variant hexamer fraction per arm/community** [GEN] — heatmap of TGAGGG / TTGGGG / TCAGGG / TTAGGG fractions; supports the methodological-vs-biological framing for the Ambrosini comparison.

---

## 7. Suggested talk slide takeaways (15-min talk)

The annotation section is supporting material rather than the main result; budget **2 slides** (≈90 s).

- **Slide A — "Every subtelomere carries the same building blocks."**
  - 173,881 gene annotations (374 unique genes) and 21,544 TAR1 entries across 15,668 PHR sequences from 465 near-complete assemblies.
  - TAR1 in **all 41 arms**, near-saturating outside PAR1; PAR1 is essentially TAR1-free (chrX_p 0.3%, chrY_p 1.1%) — consistent with PAR1 being a crossover region, not a repeat-mediated exchange region.
  - Visual: Fig M1 (TAR1 prevalence + position landscape).

- **Slide B — "Internal (TTAGGG)n islands are short, degenerate, ubiquitous."**
  - 18,352 islands in 53.1% of sequences; median 79 bp; only 32% are pure canonical TTAGGG; same three variant motifs Ambrosini et al. (2007) flagged.
  - Best read as ancient telomeric remnants embedded by past subtelomeric duplication, not a current exchange-driving feature (no significant cross-arm vs self-arm difference, p=0.045 with effect ≈ 0).
  - Visual: Fig M2 (length histogram + canonical-fraction distribution).

**Optional cut-throughs if time permits**
- Acrocentric p-arms have intermediate TAR1 prevalence (~73–79%) and the deepest TAR1 (medians 179–196 kb), reflecting their unusually long PHR regions.
- Terminal telomere length differs across communities (Kruskal-Wallis p=3.2e-15) but the correlation with internal island count is trivial (ρ=−0.056) — flag as "interesting, mechanism unclear" to defer in a 15-min talk.

**What NOT to claim on stage**
- Do **not** assert "internal islands mark duplicon boundaries" — the test in this section is against the PHR *outer* boundary; Ambrosini's internal-duplicon claim remains untested here.
- Do **not** present cross-arm vs self-arm island difference as a real finding; the p=0.045 result has trivial effect size and is excluded from the conclusion.
