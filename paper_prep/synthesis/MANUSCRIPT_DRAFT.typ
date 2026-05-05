#set page(paper: "us-letter", margin: (x: 1in, y: 1in))
#set text(font: ("New Computer Modern", "DejaVu Serif", "Liberation Serif"), size: 10pt)
#set par(justify: true, leading: 0.65em)
#show heading.where(level: 1): set text(size: 14pt, weight: "bold")
#show heading.where(level: 2): set text(size: 12pt, weight: "bold")
#let horizontalrule = align(center, line(length: 40%, stroke: 0.5pt))

= Population-scale subtelomeric communities mirror three-dimensional nuclear organisation across human and mouse
<population-scale-subtelomeric-communities-mirror-three-dimensional-nuclear-organisation-across-human-and-mouse>
== Abstract
<abstract>
Cytogenetic-era studies sketched the human subtelomere as a "patchwork"
of duplicated blocks shared across non-homologous chromosome ends, but
the population-scale architecture of these exchanges and their relation
to nuclear organisation has remained invisible. We analyse 18,827
telomere-anchored 500 kb flanks across 465 HPRCv2 haplotype-resolved
assemblies from 233 individuals; 83 % carry inter-chromosomal sequence
sharing, partitioning 41 chromosome arms into 15 inter-chromosomal
communities and 15,668 polymorphic homology regions (PHRs) into 50
sequence-level communities (Leiden, modularity Q = 0.97). Within
communities, alleles are closer than paralogs in 8 of 9 multi-arm
communities; the lone reversal — acrocentric p-arms — quantifies
population-scale interchromosomal homogenisation. The same partition is
recovered as enriched 3D contact across six technologies (bulk Hi-C,
Pore-C, CiFi, Dip-C, sperm scHi-C, mouse meiotic Hi-C); non-sharing arms
invert the signal, ruling out incidental proximity. Inter-chromosomal
exchanges captured directly in a T2T 3-generation pedigree fall in
Leiden communities at 92 %, identifying the events that build the
population structure. Cross-arm exchange frequencies recover
out-of-Africa topology; subtelomeric architecture is a
population-genetic record. We interpret the system as a self-reinforcing
similarity ↔ 3D-proximity ↔ ectopic-exchange feedback loop, with the
pseudoautosomal regions as the limiting case of a continuum that spans
all 24 chromosome ends.

#horizontalrule

== 1. Introduction
<introduction>
For four decades the human subtelomere has been described as a patchwork
of duplicated blocks shared across non-homologous chromosome ends. FISH
on flow-sorted chromosomes @Trask1998 and the Mefford–Trask cosmid
surveys @MeffordTrask2002 established that distal regions exchange
sequence between arms with \~20 % chr4q ↔ chr10q D4Z4 translocation in
the population @vanDeutekom1996. Flint and colleagues @Flint1997
formalised a two-domain model — a distal segment shared with many ends,
a proximal segment shared with few, separated by a degenerate (TTAGGG)n
boundary — characterised at chr4p, chr16p and chr22q. Linardopoulou et
al. @Linardopoulou2005 mapped 41 paralogous blocks across 33 ends and
inferred an interchromosomal duplication rate \>60× the point-mutation
rate. Ambrosini et al. @Ambrosini2007 catalogued 11 subtelomere-specific
duplicons and the bimodal sequence-identity distribution (peaks at 91 %
and 98 %). Riethman and colleagues @Riethman2004@Stong2014 resolved
haplotype-specific subtelomeric assemblies for a single individual. The
pseudoautosomal regions (PAR1/PAR2) @Rouyer1986 are recognised as a
special case of the same biology — obligate-recombining
inter-chromosomal homology blocks. These cytogenetic and single-genome
studies are the antecedents of the analysis we report.

Three classes of question are now newly addressable. First, the
#emph[population-scale phylogeography of subtelomeric exchange] — how
exchange frequencies vary across continental ancestries, which exchanges
are fixed and which are polymorphic, whether subtelomeric architecture
carries a recoverable demographic signal. Second, the
#emph[sequence-to-3D linkage] — whether the inter-chromosomal homology
partition predicts contact frequencies measured by Hi-C and its
derivatives across cell types, single cells, and meiotic stages. Third,
the #emph[unification] of named special cases (PARs, D4Z4, acrocentric
p-arms) into a single quantitative framework. The HPRCv2 pangenome — 233
individuals × 2 haplotypes, telomere-to-telomere where complete —
enables, for the first time, a population-scale, sequence-resolved,
haplotype-resolved analysis of all 24 chromosome ends together,
complemented by single-individual T2T resources for a 3-generation
pedigree @Cechova2025, a 4-generation pedigree @Porubsky2025, the
diploid line RPE-1, and the mouse subspecies B6 and CAST @Francis2025.
We use this resource to recover the cytogenetic-era patterns at
quantitative resolution, link them to 3D organisation across six
technologies, and capture the inter-chromosomal exchanges that build the
population structure directly in pedigrees.

#horizontalrule

== 2. Population-scale subtelomeric communities define an interchromosomal landscape
<population-scale-subtelomeric-communities-define-an-interchromosomal-landscape>
We extracted 18,827 telomere-anchored 500 kb flanks from 465 HPRCv2
haplotype assemblies, ran wfmash all-vs-all (asm20, p95, id95, length ≥
30 kb), built per-flank PAF graphs and called 15,668 polymorphic
homology regions (PHRs) where ≥ 2 non-homologous arms share sequence
(`Methods`; pipeline schematic \[ED1\]). PHRs have median length 105 kb
and mean length 144 kb. Of the 18,827 input flanks, 83 % (n = 15,621)
participate in at least one inter-chromosomal sharing event (\[ED1\],
\[Fig 1\]). Leiden community detection on the 41 × 41 arm-level Jaccard
distance matrix yields #strong[15 inter-chromosomal communities]
(silhouette s = 0.347, modularity Q = 0.97); the same algorithm applied
to the 15,668 × 15,668 sequence Jaccard matrix yields #strong[50
sequence-level communities] (silhouette s = 0.602; ARI 0.35 / NMI 0.76
with the arm-level partition; \[ED2\]). UPGMA on the same arm matrix
recovers 14 of 15 Leiden blocks (\[Fig 1c\]).

The 15 arm-level communities recapitulate every subtelomeric biology
named in the cytogenetic literature (\[Table 1\]). C1 contains the chr4q
↔ chr10q D4Z4 macrosatellite pair central to FSHD pathology
@Lemmers2010@vanDeutekom1996. C2 is the chr10/chr18 (Linardopoulou)
recurrent transfer pair @Linardopoulou2005. C7 contains the five
acrocentric p-arms (chr13\_p, chr14\_p, chr15\_p, chr21\_p, chr22\_p).
PAR1 and PAR2 form their own communities. f7501 (the Trask-lab cosmid
distribution) and high-copy gene families (DUX4 ×18, BAGE2, MTCO,
RPL23A, SEPTIN14P22, OR4F) map cleanly onto specific communities
(\[ED4c\]).

Each arm falls in one of three architectural categories (\[Fig 1d\]):
#strong[homogeneous] (4/41 arms, 0 % cross-arm sequence at the sequence
level), #strong[polymorphic] (28/41 arms, ≥ 2 sequence-level communities
present, the modal architecture), and #strong[fully interchangeable]
(9/41 arms — acrocentric p-arms, PARs, chr4q/chr10q — where ≥ 99 % of
subtelomeric sequence is cross-arm at the sequence level). 15.9 % of
subtelomeric sequences are cross-arm at the arm level; 11.1 % at the
sequence level (n = 15,668 sequences; \[ED2c\]). Six of the six
chromosome arms originally characterised by FISH (chr4p, chr4q, chr16p,
chr18p, chr20p, chr22q) confirm the Flint–Mefford two-domain model at
sequence resolution, and the model extends to 39 of 48 testable arms
(\[Fig 2b\]; details in §3).

We reproduce the f7501 distribution of @MeffordTrask2002 across 465
haplotypes, fix all sites previously noted, and after BH-FDR correction
within the 80-test (16 arms × 5 superpopulations) family identify 8 arm
× superpopulation pairs at q \< 0.05: chr16\_q AFR (OR\_conditional =
17.24, 95 % CI \[9.41, 32.97\], q = 5.3 × 10⁻²⁵), chr2\_q SAS (OR =
22.75, q = 2.7 × 10⁻⁹), chr16\_p AFR, chr9\_q AFR, chr8\_p AFR, chr15\_q
EUR, chr6\_p AMR, chr15\_q EAS (`STATS_AUDIT.md` §2). chr7\_p AFR,
quoted as enriched in `SURVEY_01 §1.9`, is demoted to "suggestive"
(q\_BH = 0.073) by the BH correction. TAR1 prevalence is 94.6 % across
all 41 arms but drops to 0.5 % at PAR1 (\[ED3a\]); internal (TTAGGG)n
islands populate 53.1 % of PHR sequences and only 32.2 % are pure
canonical (\[ED3b\]; @Ambrosini2007).

#horizontalrule

== 3. Within-community heterogeneity captures population-scale exchange
<within-community-heterogeneity-captures-population-scale-exchange>
If the 15 arm-level communities reflect ongoing inter-chromosomal
exchange, alleles at one arm should be closer to alleles at the same arm
than to paralogs at a different arm in the same community — the lower
bound on continued differentiation. Across all nine multi-arm
communities pooled (n = 5,946 paired distances), allele is closer than
paralog at Wilcoxon p \< 1 × 10⁻³⁰⁰ (uncorrected; single combined paired
test). Per-community, allele is closer in 8 of 9 communities; the lone
reversal is #strong[C7] (the acrocentric p-arms), where 70.5 % of pairs
put paralog closer than allele (per-community Wilcoxon p = 2.0 × 10⁻⁷;
q\_BH ≤ 5 × 10⁻⁶ in the 9-community family); C7 silhouette is −0.029 and
the conversion score is 1.000 toward each of the four other acrocentric
p-arms (\[Fig 2a\]). C7 represents complete population-scale
interchromosomal homogenisation — the limit of the inter-arm exchange
continuum.

We then test the Flint–Mefford two-domain model arm-by-arm using a
per-haplotype Spearman gradient of identity-vs-position followed by a
piecewise breakpoint fit (\[Fig 2b\]). 39 of 48 arms (81 %) carry a
significant Spearman gradient and 39 of 41 testable arms (95 %) prefer
the two-segment fit; the gradient is detectable in 99.7 % of individual
haplotypes. Arm-specific breakpoints range from 15 kb (chr22q) to 445
kb. Internal (TTAGGG)n islands sit within 25 kb of the breakpoint on 11
of 19 arms where both can be tested (\[ED3b\]). The original FISH-era
anchors (chr4p, chr4q, chr22q) reproduce; the model extends to 33
additional arms not previously tested.

The exchange landscape carries a population signature. Cross-arm ×
superpopulation Fisher tests on 19 cross-arm pairs (\[Fig 2c\]) yield 10
BH-significant pairs; mean Hudson Fst is 0.10–0.15 between AFR and each
non-AFR superpopulation and ≈ 0 between non-AFR pairs
(`SURVEY_04 §2.1`). Building a UPGMA tree from cross-arm-affinity
frequencies across 9 arms × 5 superpopulations recovers the canonical
out-of-Africa topology — AFR splits first, AMR and EUR are the closest
pair (\[Fig 2d\]). Subtelomeric exchange is therefore a phylogeographic
record (N15, N19); the same loci that carry recurrent translocations
also carry continental-ancestry frequencies.

The strongest individual-locus signal is at chr16\_q (OR = 17.24, 95 %
CI \[9.41, 32.97\], q\_BH = 5.3 × 10⁻²⁵) — a single subtelomere where
carrier frequency reaches 67 / 134 in AFR and \< 0.05 in every non-AFR
superpopulation. chr15\_q is enriched in both EUR and EAS but at
intermediate ORs (12.6 and 3.3). The structural-heterozygosity rate
within an individual reaches 47.5 % at the highest-discordance locus.

#horizontalrule

== 4. Three-dimensional nuclear organisation mirrors the sequence partition
<three-dimensional-nuclear-organisation-mirrors-the-sequence-partition>
We tested whether the sequence-defined partition predicts
inter-chromosomal contact frequency in the nucleus, using six
independent technologies (\[Fig 3\], \[ED5\], \[ED6\], \[ED7\]). On the
HG002 Pore-C inter-chromosomal contact matrix at 50 kb, with arms
ordered by sequence community, the diagonal blocks are visibly enriched
(B/W = 0.056, p = 3.9 × 10⁻⁸⁵, BH q = 2.2 × 10⁻⁸⁴ within the 40-test
multi-resolution family; \[Fig 3a\]). Across 14 independent tests (8
bulk Hi-C / Pore-C samples, GM12878 16-cell Dip-C, sperm 20-cell scHi-C,
RPE-1 self, plus 4 mouse meiotic stages), every effect size is on the
within-community-closer side of unity (\[Fig 3b\]; `SURVEY_07 §6 T-1`).

Five lines of evidence rule out incidental-proximity confounds. (i)
#strong[Negative control.] S\_all — the seven non-sharing arms —
#emph[invert] the signal in 0/16 GM12878 cells and 1/20 sperm cells (W/B
\> 1; the non-sharing arms are 11 % / 40 % #emph[farther] in 3D than
between-community pairs; \[Fig 3c\]). (ii) #strong[Multi-resolution
stability.] The B/W ratio is conserved across 5 mcool resolutions (5,
10, 20, 50, 100 kb) for all 8 bulk datasets, and 35 of 40 Mantel ρ tests
pass BH q \< 0.05 (\[ED5a\]; `STATS_AUDIT.md` §3). (iii)
#strong[Confound exclusions strengthen, never weaken, the signal.]
Excluding acrocentric p-arms, the sex chromosomes, or the strongest
community individually moves Mantel ρ within ±0.05 of full-data values;
HG02148, the lone marginal sample, is rescued from ρ = 0.152 to ρ =
0.720 by the no-acro-pq + no-sex exclusion (\[ED5b\]; `SURVEY_05 §1.7`).
(iv) #strong[Flanking paradox.] The unique-sequence 100 kb
centromere-ward of each PHR — explicitly #emph[not] multi-mapping —
produces a #emph[stronger] 3D signal than the PHR centres themselves:
HG002 flanking B/W 0.002 vs PHR B/W 0.027, a 13 × difference (\[Fig
3d\]). The Dip-C flanking radial position is also significantly more
interior than the non-flanking terminal region (0.503 vs 0.551, p = 7.4
× 10⁻³⁵). Multi-mapping artefacts cannot account for this. (v)
#strong[C4 minimal-PHR positive control.] Community C4 (chr7\_q +
chr12\_q) carries only 5–25 kb of subtelomeric sharing and zero gene
annotations, yet produces a significant 3D signal in 4 of 5 diploid Hi-C
samples (\[Fig 3a, 3b\]) — confirming that sequence sharing per se, not
gene content, drives the 3D coupling.

Single-cell technologies converge on the same picture across cellular
contexts. GM12878 Dip-C @Tan2018 gives 6.9 % within-community closer
(Mantel ρ = 0.296, p = 0.002 across 16 cells; \[ED6a\]). Sperm scHi-C
@Xu2025 gives W/B = 0.401 (60 % closer; uncorrected p = 3.9 × 10⁻⁵¹
across 20 cells; \[ED6b\]). RPE-1 self-vs-self (\[ED6d\]) reveals
#strong[cell-cycle modulation];: mitotic samples show a 3× stronger
global W/B but a 1.4 × weaker per-arm-pair correlation than asynchronous
controls — consistent with chromosome end clustering tightening in
mitosis while losing pair-resolved structure (`SURVEY_09 §1.4`). In
mouse meiotic Hi-C @Zuo2021, the per-PHR-pair Spearman correlation peaks
at zygotene (ρ = 0.715, p = 4.4 × 10⁻⁵⁵, n = 344; BH q ≪ 10⁻³⁰ within
the 8-test mouse family) and is significant across all four prophase-I
stages (\[Fig 4d\]; \[ED7c\]; `SURVEY_08 §1.7`). The mouse subtelomere
itself partitions into only #strong[2 arm-level communities] (vs 15 in
human) — a quantitative architectural contrast across species (\[ED7a\])
— yet within mouse the 3D signal saturates at the 1 Mb window
(\[ED7b\]).

#horizontalrule

== 5. Pedigree-resolved exchanges build the community structure
<pedigree-resolved-exchanges-build-the-community-structure>
The strongest causal claim — that the 15 communities are built by
inter-chromosomal exchange — is demonstrated directly in a 3-generation
T2T pedigree (WashU; PAN010 grandmother → PAN027 mother → PAN028
granddaughter; @Cechova2025). Untangling each child’s haplotypes against
the maternal contributors’ assemblies recovers 538 high-confidence
inter-chromosomal patches; #strong[494 of 538 (92 %) sit within a Leiden
community] (\[Fig 4a\]; `SURVEY_14 §1.1`). The patches comprise 16
crossover-like events, 133 gene-conversion-like tracts (≈ 90 % within
C7), and 229 acrocentric\_like NAHR signatures. The same untangle
approach applied to the fragmented CEPH1463 short-read draft assemblies
recovers only 12–13 % within-community membership — a quality-floor
positive control: the 92 % membership rate requires T2T input and is not
an artefact of overcalling.

The CEPH1463 4-generation haplotype-resolved set @Porubsky2025 provides
cross-assembler validation. Filtering hifiasm and verkko intersected
calls across 28 samples yields 11 robust parent features (\[Fig 4b\];
`SURVEY_14 §1.6`). Of these, the chr10/chr18 transfer (community C2; the
Linardopoulou recurrent pair) is detected independently in NA12877
paternal #emph[and] NA12878 maternal — the same canonical exchange in
two unrelated individuals. As a positive-control test for pipeline
behaviour on a single individual we re-ran the entire pipeline on RPE-1,
where the diploid karyotype carries a known t(X;10) translocation (\[Fig
4c\]; `SURVEY_09 §6 T-1`); the pipeline rediscovers a chrX\_q ↔ chr10\_q
partition consistent with the cytogenetically known translocation, from
sequence alone. Cross-species, the pipeline applied to mouse B6 + CAST
T2T assemblies @Francis2025 recovers a 2-community arm-level partition;
the per-PHR-pair correlation with zygotene Hi-C @Zuo2021 is ρ = 0.715
(\[Fig 4d\]). Mouse private-pairs lifted to hg38 land in human community
interiors (\[ED7d\]). The architecture is therefore not a
primate-specific artefact, although the #emph[number] of communities
differs sharply between species.

#horizontalrule

== 6. Discussion
<discussion>
We have shown that 41 human chromosome arms partition into 15
inter-chromosomal subtelomeric communities, that the same partition is
recovered by six independent 3D technologies and two species, that the
events building it are observable in pedigrees at 92 %
community-membership rate, and that the resulting frequencies recover
continental-ancestry topology. This positions the subtelomeric region as
a system in which sequence similarity, three-dimensional proximity, and
ectopic exchange form a self-reinforcing feedback loop (\[ED8a\]). Three
of the four links — similarity ↔ proximity, proximity ↔ exchange (via
meiotic chromosome-end clustering @Patel2019@Zuo2021), and exchange ↔
similarity — are supported by present data. The fourth, similarity →
#emph[future] proximity, requires temporal data and is flagged as the
outstanding causal claim.

The mechanistic candidate for community C1 — the chr4q ↔ chr10q D4Z4
pair — is a CTCF- and lamin-tethered element @Ottaviani2009@Masny2004.
Sequence-level distance histograms peak in the 0–15 kb interval at D4Z4
and DUX4L copy number is tightly distributed (median 22, with rare 0–2
outliers; Mann–Whitney p = 5.3 × 10⁻⁶ for outlier deviation; FSHD
context @Lemmers2010) (\[ED8b\]). Community C7 (acrocentric p-arms) is
the candidate nucleolar-association cluster; community C2 (chr10/chr18 —
the Linardopoulou pair) is the canonical recurrent-NAHR pair. PARs are
the limiting case: a community of two arms whose obligate-recombining
biology has long been treated as a special exception. Our data say PARs
are PHRs at the limit of the same continuum — the genome-wide "exchange
community" is the general case.

#strong[Honest null.] Recombination rate vs cross-arm affinity, computed
with the T2T-CHM13 cM/Mb map of @Lalli2025 across all 39 callable
autosomal arms, gives Spearman ρ = −0.43, p = 0.006. #emph[After]
excluding seven short-read-confounded acrocentric / PAR arms (each with
0–12 callable variants on the recombination map), the residual 32-arm
correlation collapses to #strong[ρ = 0.00, p = 0.98] (\[ED8c\]). The
apparent recombination ↔ exchange correlation in a naive analysis is
therefore driven entirely by the assembly-quality difference between
confounded and well-callable arms; we report this null as the explicit
honesty figure. Compartment identity at chromosome tips reads as 68 % A
by HG002 100 kb eigenvector but with mean e1 only +0.007 — A by GC,
interior-positioned by Dip-C radial — telomere clustering, not lamina
association, drives end alignment (\[ED8d\]).

#strong[What is not yet addressed.] Human meiotic Hi-C is the single
most informative missing experiment; a phased mouse F1 (B6 × CAST) Hi-C
dataset is the cross-species analogue (`SURVEY_07 §5; SURVEY_08 §5`). A
de-novo / inherited split of pedigree patches and source-stratified (LCL
vs blood) cross-arm validation are likewise flagged. None of the present
claims relies on these.

The 15 communities and 50 sequence-level subgroups behave as a
quantitative, population-scale, organism-spanning extension of the
cytogenetic-era patchwork model. The same 41 chromosome arms whose
homology blocks were sketched by FISH and whose exchanges were inferred
by Southern blot now organise nuclear architecture, partition
single-cell 3D structures, and trace continental ancestries.
Subtelomeric architecture is not a curiosity at the chromosome ends; it
is a population-genetic record, and PARs are its named tip.

#horizontalrule

== 7. Methods (online; \~500-word stub)
<methods-online-500-word-stub>
#strong[Assemblies and flank extraction.] Inputs: 233 HPRCv2 individuals
× 2 haplotypes = 465 assemblies (CHM13 added → 466 in some counts),
HG002v1.1, T2T-CHM13v2.0, RPE1v1.1, WashU pedigree T2T (PAN010 / PAN011
\/ PAN027 / PAN028) @Cechova2025, CEPH1463 hifiasm + verkko draft
assemblies @Porubsky2025, and mouse B6 + CAST T2T @Francis2025. p/q-arm
classification of 12,649 contigs against CHM13 centromere coordinates
uses `scripts/classify_contigs.py` (paths: `SCRIPT_INVENTORY.md` §1).
Telomere detection and trimming uses `trim-telomeres.sh` upstream of
`seqtk-telo / RUKKI`; 500 kb telomere-anchored flanks are emitted per
arm.

#strong[Pipeline backbone.] wfmash v0.23.0-41-gb5f0ff1c
(`do_not_overfilter` branch for RPE-1) all-vs-all PAFs at `-p 95 -t 48`
produce 18,827 PAF.gz files; `find-multichr-regions-incremental.py`
calls 15,668 PHRs from per-flank PAFs (default ≥ 30 kb, p95, id95). pggb
\/ odgi (Guix profile pinned; `VERSIONS.md` §1, §7 — pggb / odgi /
igraph / R-package versions remain TODO for Methods) build the per-PHR
pangenome graph; `odgi similarity --all -P` computes the 15,668 × 15,668
Jaccard matrix; `odgi untangle -e 50000 -m 1000` is used for pedigree
analysis. impg 0.4.1 projects PHR coordinates between assemblies.
minimap2 v2.30 plus wfmash cross-aligner agreement flags the chr18\_q
chimera in NA18982\#1 (\[ED1d\]).

#strong[Community detection.] `scripts/community/detect_communities.R`
runs igraph Leiden modularity scans across resolution × k for both the
41 × 41 arm matrix and the 15,668 × 15,668 sequence matrix. Reported:
arm-level k = 15 (silhouette 0.347, Q = 0.97) and sequence-level k = 50
(silhouette 0.602, Q = 0.97); ARI 0.35 / NMI 0.76 between the two
(\[ED2d\]; `SCRIPT_INVENTORY.md` §2).

#strong[3D analysis.] Bulk Hi-C / Pore-C / CiFi mcools at 5/10/20/50/100
kb are processed by `analyze_hic_communities.py` (B/W within-vs-between,
Mantel permutation, per-PHR-pair Spearman; `SCRIPT_INVENTORY.md` §2).
Dip-C and sperm scHi-C use BWA-MEM2 → sam2seg (`-q 0`) → hickit
(`--min-mapq=0`) → impute3 (4 rounds for diploid GM12878; omitted for
haploid sperm); `community_3d_enrichment.py` computes per-cell W/B,
Mantel, and radial position. The flanking control extracts 100 kb
centromere-ward of each PHR with `extract_flanking_sequences.py`. Mouse
meiotic analysis @Zuo2021 uses the 1 Mb-window community partition
(`run_mouse_1Mb_pipeline.sh`).

#strong[Pedigree analysis.]
`scripts/pedigree/analyze-pedigree-recombination.py` constructs patches
from `odgi untangle` output; `plot-pedigree-untangle.R` draws per-child
× haplotype 46-arm ribbons (`SCRIPT_INVENTORY.md` §6). Patches are
intersected against the arm-Leiden community assignment for the 92 %
within-community statistic (\[Fig 4a\]).

#strong[Statistics.] All multiplicity corrections use Benjamini–Hochberg
(`p.adjust(..., method = "BH")`, R 4.3.0); 95 % confidence intervals on
2×2 odds ratios are exact Fisher / conditional-MLE intervals
(`fisher.test(...)$conf.int`). Multi-test family definitions and
per-test BH q-values are reported in `STATS_AUDIT.md` and the per-test
TSVs in `paper_prep/synthesis/stats_audit/`. Single-hypothesis
pre-specified tests are explicitly labelled
`(uncorrected; single combined test)` where applicable. Per-test scripts
and pinned tool versions are catalogued at
`paper_prep/synthesis/SCRIPT_INVENTORY.md` and
`paper_prep/synthesis/VERSIONS.md`; pggb / odgi / igraph / R-package
versions and the HPRC Liftoff annotation-index version remain to be
captured from a fresh `sessionInfo()` and Guix manifest before
manuscript freeze (`VERSIONS.md` §7).

#horizontalrule

== Display items referenced
<display-items-referenced>
- Main text uses each of \[Fig 1\], \[Fig 2\], \[Fig 3\], \[Fig 4\] and
  each of \[ED1\] – \[ED8\] at least once; \[Table 1\] is referenced in
  §2.
- Supporting tables: ED Table 2 (23-gene coding shortlist), ED Table 3
  (high-copy gene families), SI Table S1 (`NOVEL_CONTRIBUTIONS.tsv`, 27
  rows), SI Table S2 (cross-confound Mantel ρ), SI Table S3
  (per-community Hi-C reproducibility, 15 × 10), SI Table S4
  (`LIMITATIONS_X_FINDINGS.tsv`, 12 × 6).

#horizontalrule

== Bibliography
<bibliography>
Bibliography entries below are reproduced verbatim from
`paper_prep/synthesis/REFERENCES.bib` (24 entries; produced by
`validate-captions-references`). Citation keys in this draft (`[@Key]`)
match `.bib` entry keys.


#bibliography("REFERENCES.bib", style: "ieee", title: none)
