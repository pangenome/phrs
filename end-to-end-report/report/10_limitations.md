## Sample composition and limitations

### Sample composition

**Key metrics.** 232 HPRCv2 individuals: AFR = 67 (28.9%), EAS = 52 (22.4%), AMR = 44 (19.0%), SAS = 37 (15.9%), EUR = 32 (13.8%). Population-level findings (the cross-arm affinity section) should be interpreted with this imbalance in mind. Within-superpopulation heterogeneity is not modeled.

### Methodological limitations

1. **Identity threshold**: The 95% minimum identity threshold for both all-vs-all alignment and inter-chromosomal region detection means that more divergent inter-chromosomal homology (older exchanges) is invisible. The 15.9% cross-arm affinity rate and the set of 7 arms with no inter-chromosomal signal represent lower bounds, not complete inventories. Ambrosini et al. (2007) identified a bimodal distribution of duplicon identity with peaks at 91% and 98%; the 95% threshold captures primarily the recent high-identity peak and misses most of the older 91% peak, which includes many olfactory receptor and immunoglobulin-related duplicon families. Note that Ambrosini et al. (2007) distinguished "subtelomere-only" duplicon blocks (Table 1, found exclusively at subtelomeric regions) from "subterminal" duplicon families (Table 2, positioned adjacent to terminal TTAGGG tracts but sometimes having non-subtelomeric copies as well); this distinction is not modeled in the present analysis, which treats all inter-chromosomal similarity uniformly.

2. **Flank size**: The 500 kb maximum flank extraction may truncate longer subtelomeric similarity regions. No sensitivity analysis on this parameter was performed.

3. **Region length threshold**: The 3 kb minimum output region length and 5 kb window/step size set a floor on detectable inter-chromosomal regions. Shorter shared segments would be missed.

4. **Assembly quality**: Subtelomeric regions are among the most difficult to assemble. Assembly gaps, collapses, or errors near telomeres could affect sequence content and inter-chromosomal signal. The telomere-presence filter mitigates but does not eliminate this risk. The chr18_q chimeric contig in NA18982#1 (the inter-chromosomal detection section) illustrates this concern: JBKABS010000018.1 fuses chr18 with 966 kb of chrX PAR1 across a 100 bp NNN scaffold join, and no separate chrX contig exists in this haplotype.

5. **Community detection resolution**: The 50-community Leiden solution (k-NN = 75, resolution = 0.8) is one of a family of possible solutions. Different resolution parameters would yield different community numbers and compositions. The solution was selected by modularity optimization within a 5–50 community range, but modularity does not guarantee biological correctness.

6. **Small sample sizes**: Some population structure findings rest on small sample sizes: chrY_p has only 10 cross-arm sequences out of 92 total in the population enrichment test (p_adj=0.028), and chr15_p discordance is based on 6/22 individuals (27.3%). These should be considered preliminary pending validation with larger or independent datasets.

7. **Exchange timing**: Cross-arm affinity demonstrates that exchange has occurred but cannot date individual events. High discordance rates are consistent with recurrent exchange but could also reflect a single ancient event still segregating. Distinguishing these scenarios requires trio/family data or population genetic modeling.

8. **Somatic exchange in cell lines**: Mefford & Trask (2002, citing van Overveld et al. 2000) noted that "some people are mosaic for 4q/10q subtelomeric translocations, which indicates that subtelomeric sequences can interchange in somatic cells." Since most HPRCv2 assemblies derive from lymphoblastoid cell lines (LCLs), some cross-arm affinity — particularly at chr4_q/chr10_q (C1) — could reflect somatic exchange during cell culture rather than germline polymorphism. This caveat applies to all LCL-derived assemblies and is difficult to quantify without matched blood-derived controls.

### 3D validation limitations

9. **Somatic vs meiotic context**: All 3D data (Hi-C, Pore-C, Dip-C) captures somatic interphase organization. The meiotic bouquet stage — when ectopic recombination between subtelomeric arms would occur — has never been captured by Hi-C in humans. The observed somatic 3D signal is interpreted as a residual of meiotic chromosome organization (Rabl configuration — the retained centromere-telomere polarity from the preceding cell division, where centromeres cluster at one nuclear pole and telomeres at the other), but this interpretation is indirect. Human meiotic Hi-C remains the single most informative missing experiment.

10. **Sample size (N=6 Hi-C)**: Five diploid HPRC samples plus one haploid (CHM13) are sufficient to demonstrate that the 3D signal exists and is reproducible, but insufficient for population-level claims about 3D variation. CHM13 shows no significant communities, primarily due to reduced power (37 arms vs 75 in diploid samples).

11. **GM12878 cell line**: The Dip-C data uses GM12878 (an EBV-transformed B-lymphoblastoid cell line), which has an abnormal karyotype and may not represent normal nuclear organization. PBMC results provide a primary-cell control but with fewer cells.

12. **hg19/T2T coordinate incompatibility**: Dip-C data (Tan et al. 2018) uses hg19 coordinates; PHR regions are defined on CHM13/T2T. Coordinate projection via impg partially mitigates this, but coarse 50 kb resolution introduces noise, particularly near assembly gaps and subtelomeric regions where hg19 and T2T differ most.

13. **Multi-mapping at PHR intervals**: Duplicated PHR sequences cause Hi-C/Pore-C reads to multi-map between community partner arms, inflating apparent inter-chromosomal contacts. All four technologies (Hi-C, Pore-C, CiFi, Dip-C/sperm) disable MAPQ filters and retain multimappers, but each multimapped read keeps exactly **one randomly-chosen alignment** — aligners run in default mode (no `-k`/`--all`), so reads are not duplicated across all valid positions. The noise is symmetric (equal probability of misplacement in either direction), so aggregate enrichments hold but individual pair-level contacts in repetitive regions are unreliable. The flanking-region control (unique sequence) addresses this by demonstrating that the signal persists — and is stronger — in regions without multi-mapping. However, the PHR-specific enrichment values (B/W 0.027–0.074) cannot be separated from multi-mapping contribution.

14. **Confound controls**: The comprehensive exclusion analysis (the acrocentric exclusion control) tests five exclusion sets (acrocentric p-arms, sex chromosomes, both, all acrocentric p+q + sex, and strongest communities C1+C7+C14/C15) and demonstrates that the Mantel rho INCREASES when removing these arms (HG002: 0.657→0.790; HG02148: 0.152→0.720), ruling out nucleolar association, PAR sharing, and D4Z4 as drivers. Rabl configuration (centromere-telomere polarity causing generic telomere clustering) is addressed by the flanking analysis (the flanking analysis section): if Rabl drove the signal, flanking regions farther from the telomere tip should show weaker signal, but they show stronger signal. Chromosome size effects (small chromosomes intermingle more) are addressed by the Mantel test (the Mantel test section), which tests a continuous correlation between sequence similarity and Hi-C contact across all arm pairs — a size confound would add noise but not create a correlation between Jaccard distance and contact frequency. The significant Mantel results (HG002 Hi-C rho = 0.66, Pore-C rho = 0.49; the Mantel test section) demonstrate that the 3D signal is specifically tied to subtelomeric sequence similarity.

15. **Parameter sensitivity**: The 50 kb resolution and 500 kb flanking window were selected based on optimization but not subjected to formal sensitivity analysis. Different resolutions or window sizes could yield different enrichment values. Multi-resolution analysis at all 5 mcool resolutions (5kb, 10kb, 20kb, 50kb, 100kb) across human, RPE-1, and mouse systems (the resolution sensitivity section, the acrocentric exclusion control, the RPE-1 validation section.1, the mouse flanking Hi-C section) demonstrates that the core signal is robust to resolution choice.

16. **Fragmented assemblies produce NaN flanking values**: HG02148 and NA19036 assemblies are sufficiently fragmented at some subtelomeric regions that flanking-region coordinates fall outside contig boundaries, producing NaN values in the community-free flanking correlation. These samples are excluded from flanking analyses but included in PHR-based analyses.

17. **Dip-C cell 12 duplicate**: Cell 12 produces identical 3dg output to cell 10 (shared SRR7226706 long-insert library). Cell 12 is excluded from all analyses (16 cells used in the Dip-C section).

18. **Mouse 1 Mb PHRs mostly fill the extraction window**: At the 1 Mb extraction scale, mouse PHR regions have a median length of ~980 kb, meaning the PHR nearly saturates the window. Even larger windows (1.5–2 Mb) may reveal additional inter-chromosomal similarity beyond the current detection boundary, particularly for the p-arm community (C1, 16 chromosomal arms) where sharing extends deep into the chromosome.

---

The preceding sections established subtelomeric community structure from sequence similarity. The following sections test whether this structure has a physical counterpart in 3D nuclear organization.

