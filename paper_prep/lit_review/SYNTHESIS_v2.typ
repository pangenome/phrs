= SYNTHESIS\_v2 — Readable history of human subtelomere research, anchored to the Guarracino & Garrison Nature manuscript
<synthesis_v2-readable-history-of-human-subtelomere-research-anchored-to-the-guarracino-garrison-nature-manuscript>
#strong[Author:] lit-review-v2 (agent-937), 2026-05-06. Supersedes
`SYNTHESIS.md` (lit-review-synthesis, agent-877). #strong[Anchors:]
`paper_prep/synthesis/ABSTRACT_nature.md` (primary) and
`ABSTRACT_BoG.md` (talk variant). #strong[Bibliography:]
`paper_prep/synthesis/REFERENCES_v3.bib` (295 entries; merged from
`REFERENCES.bib`, the prior `REFERENCES_v2.bib`, and the 14 topical
`topic_NN_*.bib` files; deduplicated by DOI and (first-author surname,
year, normalized-title prefix) with canonical (prefix-free) keys
preferred so citations resolve). #strong[Companion files:]
`CHRONOLOGY_v2.md` (timeline table), `GAPS_v2.md` (open questions and
missing-substrate notes), `paper_prep/synthesis/CROSSWALK.md` (claim ↔
Andrea-report substrate map). #strong[Versioning note:] v1 deliverables
(`SYNTHESIS.md`, `CHRONOLOGY.md`, `GAPS.md`, `REFERENCES_v2.bib`) are
retained as historical record. v2 is the canonical reference moving
forward; the substantive change from v1 is that topics 05 (acrocentric
short arms / rDNA / Robertsonian, refreshed v3-opus pass) and 07
(concerted evolution / NAHR / gene conversion, refreshed v3-opus pass)
supply richer per-paper substrate to Parts II and III, and the
bibliography now incorporates the final versions of all 14 topical bibs.

#line(length: 100%)

== Foreword — how to read this document
<foreword-how-to-read-this-document>
This document is a single, end-to-end narrative history of the field
that the Guarracino & Garrison subtelomere manuscript is entering. It is
written for Erik to read in order, before a Board-of-Governors talk or a
manuscript-revision session. It is organised into five parts. Parts I–IV
each consolidate three to five of the originally-planned topical reviews
(01–14) into a chronological field history; transitions between
sub-topics are marked in the prose so each part still reads like a
unified essay. Part V then steps back to ask how the present claims
(C1–C8 of the abstract) sit in this history: which are old discoveries
reproduced at population scale, which are new, and which depend on
framing decisions that the abstract has chosen.

Citations use BibTeX keys from `REFERENCES_v3.bib`
(e.g.~`[@MeffordTrask2002]`); claim labels follow the C1–C8 convention
defined in `CROSSWALK.md`. All 14 topical reviews are committed
(`topic_01_*` cytogenetic foundations, `topic_02_*` subtelomere
structure, `topic_03_*` pseudohomology concept, `topic_04_*`
sex-chromosome PARs, `topic_05_*` acrocentric / rDNA / Robertsonian,
`topic_06_*` D4Z4 / DUX4 / FSHD, `topic_07_*` concerted evolution /
NAHR, `topic_08_*` meiotic bouquet / nuclear-envelope tethering,
`topic_09_*` Hi-C / Pore-C / single-cell / meiotic 3D, `topic_10_*`
pangenome graphs / IMPG, `topic_11_*` pedigree-based recombination
detection, `topic_12_*` HPRC v1 / v2 population pangenomes, `topic_13_*`
subtelomere population genetics, `topic_14_*` OR / OR4F gradient) under
`paper_prep/lit_review/`; their bibliographies feed `REFERENCES_v3.bib`
(295 entries) directly. Topics 05 and 07 each received a v3-opus refresh
after the v1 synthesis; that refresh is fully integrated below in Parts
II and III. The narrative names the canonical primary papers for each
part; the topical reviews provide the corresponding per-paper depth and
full citation context. `GAPS_v2.md` records bibliographic and
substantive open questions that the synthesis surfaced for future
revision passes.

#line(length: 100%)

== Part I — Cytogenetic and structural foundations (combines topics 01, 02, 04)
<part-i-cytogenetic-and-structural-foundations-combines-topics-01-02-04>
=== From banded karyotypes to fluorescent probes (1916–1990)
<from-banded-karyotypes-to-fluorescent-probes-19161990>
The story of human subtelomeres begins, paradoxically, in grasshoppers.
In 1916 William Robertson catalogued V-shaped chromosomes in Tettigidae
and Acrididae, observed acrocentric–acrocentric fusions and gave the
field the word "Robertsonian" that still describes the most common human
chromosomal rearrangement involving the five acrocentric short arms
(chr13, 14, 15, 21, 22) @Robertson1916. Acrocentric short arms — and
through them, the subtelomeric problem of "what is on the chromosome
distal to a heterochromatic block, and how does it move?" — entered
human cytogenetics as soon as a banded karyotype existed.

That banded karyotype arrived with quinacrine in 1970. Caspersson, Zech,
Johansson and Modest’s quinacrine-banding paper @Caspersson1970 was the
first method that could identify all 24 human chromosomes by their
longitudinal pattern; it gave the field G-banding nomenclature and made
it possible to talk about chromosome arms (p / q), telomeric and
centromeric extremes, and recurrent rearrangements with
chromosome-specific resolution. Hsu’s 1979 monograph collected the
surrounding cytogenetic foundations @Hsu1979. Two years after
Caspersson, Henderson, Warburton and Atwood used in situ hybridisation
to localise human ribosomal DNA to the secondary constrictions of all
five acrocentric short arms @Henderson1972, proving for the first time
that an evolutionarily conserved, highly repetitive, partially
homogenised gene family was distributed across multiple non-homologous
chromosomes. This is the original
"concerted-evolution-at-chromosome-ends" observation in humans, and it
underwrites Part II’s pseudohomology concept and Part V’s interpretation
of why C5’s acrocentric clade and C8’s "concerted evolution" framing
connect at all.

Fluorescence in situ hybridisation (FISH) entered the toolset in 1986
with Pinkel, Straume and Gray’s quantitative high-sensitivity protocol
@Pinkel1986, scaled to chromosome painting in Lichter et al.’s 1990
cosmid-FISH demonstration @Lichter1990 and made systematic in Trask’s
1991 #emph[Trends in Genetics] review @Trask1991. Together with the
discovery of the canonical telomeric (TTAGGG)n repeat — first
characterised in #emph[Tetrahymena] rDNA termini @BlackburnGall1978 and
shortly after extended to human chromosome ends — these methods
transformed subtelomeres into a measurable rather than theoretical
entity. By the early 1990s a researcher could pick a cosmid from a
chromosome-end region, tag it with a fluorophore, and ask which other
chromosome ends carried homologous sequence in interphase or metaphase
nuclei. The Trask lab’s flow-sorted-chromosome FISH programme exploited
this: it identified blocks of olfactory-receptor genes duplicated
polymorphically among chromosome ends @Trask1998 and produced cosmid
f7501 (L78442), which would later become the field’s single
most-revisited subtelomeric block.

In parallel, Brown, MacKinnon, Villasanté and colleagues described
"telomere-associated DNA" — #emph[Cell] 1990 — as polymorphic, variably
arranged sequence shared across chromosome ends @Brown1990. This is the
founding paper for what the present manuscript calls TAR1
(Telomere-Associated Repeat 1). The next year, Wilkie, Higgs and
colleagues used segregation analysis at the tip of chromosome 16p to
document a stable length polymorphism of up to 260 kb, which they argued
probably arose by occasional exchange among non-homologous subtelomeres
@Wilkie1991. Together these two papers established that chromosome ends
were #emph[exchange substrates] rather than fixed reference termini, and
they are Part I’s most important pre-1995 claim-anchored finding (it
directly motivates C4 — extended interchromosomal homology at nearly all
subtelomeres).

=== The pseudoautosomal region as the prototype (1982–2025)
<the-pseudoautosomal-region-as-the-prototype-19822025>
A decade before the modern molecular era, Burgoyne had already given
mammalian cytogenetics the obligate-X/Y-crossover argument: the
pseudoautosomal region (PAR1) must recombine in male meiosis or X and Y
will fail to segregate @Burgoyne1982. The molecular era began the same
year as Pinkel’s quantitative FISH protocol. In 1986, Goodfellow,
Darling, Thomas and Goodfellow showed that #emph[MIC2] (later CD99) is a
pseudoautosomal gene with measurable recombination from the
sex-determining locus @sexchrompars_goodfellow1986, and Rouyer, Simmler
and colleagues localised the broader phenomenon as a #emph[gradient] of
sex linkage @Rouyer1986: markers very near the chromosome tip behaved as
if autosomal (high recombination; low LD with sex), while markers
proximal to the PAR boundary behaved as sex-linked. Pritchard et
al.~used long-range restriction mapping to place MIC2 close to
Y-specific sequence @sexchrompars_pritchard1987, and Ellis, Goodfellow
and colleagues then identified an Alu-repeat insertion at the proximal
PAR1 boundary on the Y chromosome
@sexchrompars_ellis1989@sexchrompars_ellis1990 — the first molecular
description of an abrupt boundary between recombining and
non-recombining sex-chromosome territory. PAR1 thus became the first
explicitly characterised "pseudohomologous" region in the human genome,
with an obligate functional role tying physical proximity at the
chromosome tip to recombination.

The pseudoautosomal frame expanded in the 1990s. Freije, Helms, Watson
and Donis-Keller used microsatellite inheritance in CEPH pedigrees to
discover a #emph[second] pseudoautosomal region near the Xq/Yq tips
@sexchrompars_freije1992; Graves, Wakefield and Toder synthesized
comparative-mammal evidence that PAR1 and PAR2 have different
evolutionary origins (PAR1 as relic added autosomal material, PAR2 as a
more recent compound event) @sexchrompars_graves1998; Lahn and Page
placed both PARs alongside the X-conserved (XCR), X-added (XAR) and
X-transposed (XTR) regions in their evolutionary-strata model
@sexchrompars_lahn1999. Filatov then demonstrated a gradient of silent
substitution rate across the human PAR consistent with the sex-averaged
recombination gradient @Filatov2004. The sequence era refined both ends:
Ciccodicola et al.~produced the first full sequence of the Xq/Yq
pseudoautosomal interval and showed it has internally distinct
gene-regulation and evolutionary subregions
@sexchrompars_ciccodicola2000; Charchar et al.~mapped PAR2 as arising
through multiple species-specific events rather than a simple single
duplication @Charchar2003; Skaletsky et al.’s Y-chromosome reference
assembly placed PAR1, PAR2, X-degenerate, X-transposed and ampliconic
territory in one mosaic @Skaletsky2003; and Ross et al.’s X-chromosome
reference completed the picture by naming PAR1 / PAR2 / XAR / XCR / XTR
in one framework @Ross2005. Mangs and Morris consolidated this into a
current-genomics review @MangsMorris2007; Flaquer, Rappold, Wienker and
Fischer wrote the genetic-epidemiology PAR primer
@sexchrompars_flaquer2008 and the sex-specific genetic-map paper that
quantified the contrast between PAR1 (extreme male recombination) and
PAR2 (near-genome-average) @sexchrompars_flaquer2009. PAR1 is
approximately 2.7 Mb and supports the obligate X–Y crossover; PAR2 is
roughly 320–400 kb, has a different origin, and undergoes much less
recombination.

From 2014 onward the PAR literature shifted from "where are the PARs?"
to "how stable are PAR boundaries, how much sequence is required for
function, and what does direct recombination mapping look like at the
molecular level?" Hinch, Altemose, Noor, Donnelly and Myers produced the
canonical fine-scale recombination map of PAR1 from pedigrees and LD,
with intense male recombination and resolvable hotspot structure
@sexchrompars_hinch2014. Mensah et al.~then identified recurrent human
PAR1 length polymorphism caused by X-to-Y translocation — the first ePAR
description @sexchrompars_mensah2014; Cotter, Brotman and Wilson Sayres
showed that population-genetic diversity does not produce a clean cliff
at the annotated PAR1 boundary, supporting boundary plasticity
@sexchrompars_cotter2016; Poriswanish et al.~characterised recombination
hotspots within an extended human pseudoautosomal domain by sperm-based
crossover analysis @sexchrompars_poriswanish2018; Acquaviva et
al.~revealed in mouse that ANKRD31 and MEI4-related axis-elongation
factors force DSB formation into the pseudoautosomal segment so that
even a small interval can carry the obligate exchange
@sexchrompars_acquaviva2020; Fukami et al.~provided pedigree evidence
that human spermatogenesis tolerates massive size reduction of PAR1
@sexchrompars_fukami2020; Monteiro, Arenas, Prata and Amorim built the
population-genetic dynamics model for PARs @sexchrompars_monteiro2021;
Bergman and Schierup placed PAR1 in great-ape population-genomic context
with recombination-associated mutation and GC-biased gene conversion
@sexchrompars_bergman2022. The T2T era arrived with Miga et al.’s
telomere-to-telomere X chromosome @sexchrompars_miga2020, Rhie et al.’s
complete Y @sexchrompars_rhie2023 and Hallast et al.’s 43 assembled
human Y chromosomes documenting extensive PAR1 length polymorphism in
the population @sexchrompars_hallast2023; and Bellott, Hughes,
Skaletsky, Owen and Page re-examined genetic and sequence evidence to
narrow the canonical PAR1 boundary against a competing 2023
assembly-based shift claim @sexchrompars_bellott2024. Most recently
Poriswanish et al.~analysed approximately 218,300 46,XY males in UK
Biobank, estimating ePAR1 incidence, recurrent origins across multiple Y
haplogroups, reciprocal X deletions and phenotypic associations
@sexchrompars_poriswanish2025 — establishing that PAR1 is not a fixed
textbook interval but a recombination #emph[system] with recurrent
structural variants creating new pseudoautosomal sequence in living
populations.

The result is a clean prototype: PAR1 and PAR2 are bona fide
pseudohomologous regions; their recombination behaviour is functionally
required for sex-chromosome segregation or for evolutionary turnover of
male-specific gene content; their boundaries are biologically plastic at
population scale; and their existence proves that interchromosomal
sequence sharing at chromosome ends is biologically tolerated — even
#emph[required] — at the human germline level. This is the empirical
anchor for C4’s "comparable in scale to the canonical PAR2" wording:
when the present manuscript reports tens-to-hundreds of kilobases of
homology at nearly all subtelomeres comparable to PAR2, the reader is
being asked to view PAR2 not as the only such region but as one named
instance of a much broader phenomenon. It is also the empirical anchor
for C5’s positive controls: the manuscript’s NJ-tree clades C14 (Xq/Yq →
PAR2) and C15 (Xp/Yp → PAR1) are recovering the two best-known
sex-chromosome subtelomeric clades. Per the topic\_04 review’s
interpretation, PAR1 supplies the #emph[strong] mechanistic anchor
(genetically mapped, sperm-typed obligate male crossover) and PAR2 the
#emph[weaker] anchor (true Xq/Yq homology but with crossover rate near
genome average); together they calibrate the manuscript’s broader
autosomal-PHR claims by giving two named instances where recombination
is #emph[already] part of the definition.

=== Subtelomeric structural biology: TAR1, internal telomeric sequences and the two-domain model (1996–2008)
<subtelomeric-structural-biology-tar1-internal-telomeric-sequences-and-the-two-domain-model-19962008>
Brown’s TAR1 @Brown1990 and Wilkie’s chr16p polymorphism @Wilkie1991
launched a structural-biology programme that asked what the subtelomeric
repeat landscape actually looks like at the sequence level. The 1996
sequencing of the chr4q35 D4Z4 macrosatellite @vanDeutekom1996 and the
1992 mapping of 4q-tip rearrangements to FSHD @Wijmenga1992 gave the
field a disease-linked exemplar of subtelomeric exchange — the chr4q ↔
chr10q D4Z4 translocation observed in roughly 20% of typed individuals
@vanDeutekom1996. Flint, Bates and colleagues then proposed the
two-domain model in 1997 @Flint1997: a distal subtelomeric domain that
shares sequence with many chromosome ends, and a more proximal domain
shared with fewer ends, separated by degenerate internal (TTAGGG)n
tracts. This model was inferred from comparison of human and yeast
subtelomeres and was characterised at chr4p, chr16p and chr22q only; it
predicts both that interchromosomal sharing should #emph[decrease] with
distance from the telomere and that internal telomeric sequence islands
should sit at the boundary. It is the single most-tested prediction of
the Flint era and is precisely what the present manuscript’s "two-domain
Flint/Mefford model" result re-tests at pangenome scale (39/48 arms
gradient; 39/41 piecewise-breakpoint; 11/19 ITS within 25 kb of the
breakpoint — a confirmation orders of magnitude broader than the
original three arms).

Internal (TTAGGG)n repeats themselves became a research object.
Ruiz-Herrera, Nergadze and colleagues reviewed the mechanisms by which
telomeric tracts end up far from chromosome ends — interstitial
telomeric sequences (ITS) — and the conserved evolutionary signal they
carry @RuizHerrera2008. By 2008 the field accepted that subtelomeres
contain a stable architecture of (i) a very distal canonical (TTAGGG)n
telomere; (ii) TAR1 (Brown 1990) immediately proximal; (iii) a
duplicated patchwork of subtelomere-specific repeats; (iv) often,
internal (TTAGGG)n islands deeper in. In the present manuscript,
chapters 02 and 04 of Andrea’s report quantify each of these layers
across 466 haplotypes (TAR1 prevalence 94.6%; absent from PAR1 at 0.5%;
18,352 internal (TTAGGG)n islands across 8,321 sequences), and these
observations are direct population-scale validation of the Flint and
Brown structural model.

=== Riethman, Ambrosini and the assembly era (2001–2014)
<riethman-ambrosini-and-the-assembly-era-20012014>
Riethman, Xiang and colleagues integrated telomere clones with the draft
human genome in 2001 @Riethman2001 and established the systematic
subtelomeric assembly programme that would dominate the subsequent
decade @Riethman2004. The conceptual upshot was that \~80% of the most
distal 100 kb of human chromosomes consists of shared duplicated blocks;
subtelomeres are not single-copy regions punctuated by telomeres but
rather #emph[patchworks] of duplications. Ambrosini, Paul, Hu and
Riethman then catalogued these patchworks in 2007 as a duplicon-block
organisation @Ambrosini2007: 11 subtelomere-specific block families, 6
subterminal block families, 6 one-copy (TTAGGG)n-adjacent regions, and a
bimodal sequence-identity distribution with peaks at \~91% and \~98%
(the two peaks reflecting older and more recent exchange events
respectively). Riethman’s 2008 review then placed copy-number variation
at subtelomeres in the broader CNV literature @Riethman2008. Stong, Deng
and colleagues’ 2014 paper extended the assembly programme into
structural and chromatin annotation @Stong2014, establishing CTCF and
cohesin organisation across haplotype-resolved subtelomeres.

The assembly era is conceptually important for two reasons. First, it
converted the FISH-era observation that subtelomeres exchange into a
#emph[catalogue] of which sequence blocks exchange between which
chromosome ends — vocabulary the present manuscript inherits when it
names communities (C1 D4Z4-bearing 4q–10q, C2 chr10p–chr18p, C3
f7501-cluster, C5 RPL23A/WASH, etc.). Second, it exposed the limits of
single-haplotype references: even with a full-length BAC tile across a
chromosome-end region, individual subtelomeres show major copy-number,
allele-organisation and gene-content variation that a reference cannot
represent. This exposure motivates C2 (the implicit pangenome graph as a
reference-free response) and C3 (466 near-complete haplotype assemblies
as the substrate that finally makes the assembly problem tractable
across populations).

=== Long-molecule and T2T eras (2018–2022)
<long-molecule-and-t2t-eras-20182022>
Optical mapping and long-read sequencing reopened subtelomere assembly
from the technology side. Shao, Zhou and colleagues @Shao2018 used
BioNano and nanopore data to detect ongoing chromosome-end extension
events — heterozygous extensions that could be modelled as chimeric
chromosome ends. Young, Abid and colleagues’ whole-genome optical
mapping of 154 individuals showed widespread large-scale subtelomeric
structural variation undersampled in standard references @Young2020.
Grigorev, Foox and colleagues used long-read telomere sequencing to
characterise telomeric repeat-motif heterogeneity across individuals
@Grigorev2021. Adam, Ranjan and Riethman’s NPGREAT pipeline integrated
linked reads with ultralong nanopore reads to assemble CHM13
subtelomeres across representative segmental-duplication structures
@Adam2022. Logsdon et al.’s complete chr8 @Logsdon2021 and Nurk et al.’s
T2T-CHM13 assembly @Nurk2022 then closed many formerly inaccessible
terminal regions; Vollger et al.~catalogued the segmental-duplication
landscape at T2T scale @Vollger2023. The HPRC v1 draft pangenome
demonstrated that haplotype-resolved diverse assemblies could be
organised into a population reference @Liao2023.

This sequence of advances is the technical bridge to C1–C3 of the
present manuscript. The 466 HPRC v2 near-complete haplotype assemblies
are the operational endpoint of three decades of subtelomeric structural
biology: each arm is now a sequenced haplotype rather than an inferred
FISH signal, each haplotype is part of a population-scale collection
rather than a single reference, and each subtelomere can be aligned
all-to-all without imposing chromosomal partitioning. Part I of this
synthesis ends here, at the moment in 2022 when reading every
subtelomere in 233 individuals became practical. Parts II–IV pick up the
conceptual, mechanistic and methodological histories that turn that data
into the manuscript’s claims.

#line(length: 100%)

== Part II — The pseudohomology concept and its variants (combines topics 03, 05, 06, 13, 14)
<part-ii-the-pseudohomology-concept-and-its-variants-combines-topics-03-05-06-13-14>
=== The pseudohomology concept proper (1990–present)
<the-pseudohomology-concept-proper-1990present>
The central conceptual unit of this manuscript is the
#strong[pseudohomologous region (PHR)] — telomere-adjacent sequence that
is paralogously shared across non-homologous chromosome ends, in
contrast to the formal #emph[pseudoautosomal] region (PAR) which is
functionally homologous between X and Y. The PHR concept emerged in
three layers, each documented in detail in the topic\_03 review
\[paper\_prep/lit\_review/topic\_03\_pseudohomologous\_regions\_concept.md\]
and summarised here.

The first layer is #emph[architectural];: Flint et al.’s two-domain
model @Flint1997 proposed that chromosome ends share a distal multi-end
domain and a proximal few-end domain separated by internal (TTAGGG)n
tracts. The second layer is #emph[evolutionary];: Trask 1998 @Trask1998,
Mefford & Trask 2002 @MeffordTrask2002, Der-Sarkissian et al.~2002
@DerSarkissian2002, Linardopoulou et al.~2005 @Linardopoulou2005 and
Ambrosini et al.~2007 @Ambrosini2007 turned the architecture into a
dynamic patchwork of duplicated blocks moved among ends at high rate.
Mefford & Trask in particular framed subtelomeres as transitions between
chromosome-specific DNA and terminal telomere repeats with recurrent
exchange and copy-number variation, and explicitly hedged that block
frequencies might carry population-history signal even though exchange
history could make them unreliable phylogenetic markers — a hedge that
the present manuscript’s chapter 04 + chapter 12 results revisit by
recovering the out-of-Africa topology from cross-arm-affinity Fst
(chapter 12 novel contribution \#19; supports C6). The third layer is
#emph[operational];: Riethman 2004 @Riethman2004, Stong 2014 @Stong2014,
Adam 2022 @Adam2022, Nurk 2022 @Nurk2022 and Liao 2023 @Liao2023 made it
possible to #emph[measure] PHRs at sequence resolution across many
haplotypes — which is the technical premise of C2 (implicit pangenome
graph) and C3 (466 haplotypes).

The PHR concept is also bounded by alternatives. Bailey, Gu and
Eichler’s "recent segmental duplications" framework @Bailey2002, updated
by Eichler 2001 @Eichler2001, Sharp 2006 @Sharp2006, Sudmant 2015
@Sudmant2015 and Vollger 2023 @Vollger2023, placed subtelomeric exchange
within a genome-wide segmental-duplication story rather than a
chromosome-end-specific phenomenon. Stankiewicz & Lupski’s NAHR /
non-allelic homologous recombination framework
@StankiewiczLupski2002@StankiewiczLupski2010 provided the mechanistic
vocabulary of #emph[how] large duplicated repeats mediate ectopic
recombination through misalignment at homologous templates. Lovett’s
misalignment-mediated mutation review @Lovett2004 supplied the molecular
detail for slipped-strand and template-switching events that produce
sub-block variation. And Rudd, Friedman et al.~2007 @Rudd2007 showed
that even sister-chromatid exchange rates are elevated at chromosome
ends, which is a cell-biological reason chromosome ends accumulate
exchange events broadly, not merely through one mechanism.

The current manuscript inherits and refines this concept. The key
operational definition is: a PHR is the telomere-adjacent interval of a
haplotype that shares high-identity sequence with other chromosome ends
in the implicit pangenome graph; communities are inferred from the
resulting similarity graph rather than imposed from cytogenetic labels
\[paper\_prep/lit\_review/topic\_03\_pseudohomologous\_regions\_concept.md,
lines 31–35\]. The novelty of C4–C5 is strongest where all three
historical layers converge: the old two-domain and patchwork models are
recovered, known clades (PAR1, PAR2, acrocentric, 10p–18p,
f7501-cluster, 4q–10q D4Z4) are reproduced, novel clades appear, and the
population-scale structure becomes measurable for the first time.

=== The acrocentric short arms and rDNA (topic 05 — refreshed v3-opus pass)
<the-acrocentric-short-arms-and-rdna-topic-05-refreshed-v3-opus-pass>
The five acrocentric short arms (chr13, 14, 15, 21, 22) are the single
best-studied "concerted-evolution at chromosome ends" example, and the
v3-opus refresh of topic\_05 substantially deepens the substrate
available to the manuscript. Henderson, Warburton and Atwood’s 1972 in
situ hybridisation @Henderson1972 showed rDNA arrays distributed across
all five short arms; Schmickel 1973
@acrocentric_rdna_robertsonian_schmickel1973 quantified human rDNA copy
number at roughly 200 per haploid genome, the first dosage estimate for
the five-site cistrons. Krystal et al.~1981
@acrocentric_rdna_robertsonian_krystal1981 used somatic-cell hybrids to
show that nonhomologous human acrocentric chromosomes carry the same
intergenic-spacer length variants — direct molecular evidence that
sequence information moves between heterologous NORs and a foundational
test of concerted evolution; Worton et al.~1988 @acrocentric_Worton1988
mapped the orientation of the rDNA tandem arrays; Therman, Susman and
Denniston 1989 @Therman1989 catalogued Robertsonian translocations and
noted the dominance of chr13–14 and chr14–21 fusions in the empirical
distribution. The molecular mapping era then narrowed the breakpoints:
Page, Shaffer and Shaffer 1996 @acrocentric_rdna_robertsonian_page1996
localised rob(13;14) and rob(14;21) breakpoints to the satellite-III /
rDNA boundary; Bandyopadhyay et al.~2001
@acrocentric_rdna_robertsonian_bandyopadhyay2001 identified seven new
satellite III subfamilies across the five p-arms; Bandyopadhyay et
al.~2002 @acrocentric_rdna_robertsonian_bandyopadhyay2002 established
that common Robertsonian translocations are predominantly dicentric,
occur preferentially in female meiosis, and may be driven by nucleolar
association. Stults et al.~2008 @acrocentric_Stults2008 then used
optical pulsed-field gel electrophoresis to measure rDNA cluster lengths
across multigenerational pedigrees, reporting a meiotic rearrangement
frequency exceeding 10% per cluster per meiosis — among the highest
documented rates of structural change in the human genome.

Floutsakou and colleagues’ 2013 paper @Floutsakou2013 then sequenced the
#emph[distal junction] (DJ) and #emph[proximal junction] (PJ) regions
flanking the rDNA on chromosome 21, found that PJ and DJ are 95% and 99%
identical across the five acrocentrics, and identified the
distal-junction long noncoding transcript that orchestrates the
surrounding chromatin. Stimpson et al.~2014
@acrocentric_rdna_robertsonian_stimpson2014 linked nucleolar
organisation, rDNA stability and acrocentric integrity to telomere
function. Jarmuz-Szymczak et al.~2014
@acrocentric_rdna_robertsonian_jarmuzSzymczak2014 further refined
Robertsonian breakpoints using FISH on de novo and familial cases.
McStay’s 2016 review consolidated these findings into the "nucleolar
organiser regions are genomic dark matter" frame @McStay2016; van Sluis
et al.~2019 extended the architecture and confirmed near-identity of
acrocentric short arms in human and chimpanzee @vanSluis2019. Most
recently, Guarracino et al.~2023 @Guarracino2023 used pangenome-graph
methods on HPRC v1 data to demonstrate ongoing recombination between
heterologous acrocentric short arms — the canonical predecessor of the
present manuscript and the paper that establishes "subtelomeric exchange
across non-homologous chromosomes is detectable in pangenome graphs." de
Lima, Guarracino et al.~2025 @deLima2025 completed the assemblies of
three Robertsonian chromosomes (rob(13q14q), rob(14q21q)) and identified
the SST1 macrosatellite array — inverted on chromosome 14 — as the
recombining homology that produces the dominant fusion product. Cechova,
Hartley et al.~2025
@acrocentric_rdna_robertsonian_cechovaHartley2025apes extended the
analysis across complete great-ape acrocentric assemblies, showing
recurrent segmental duplication, FRG1 amplification and centromere
repositioning across 25 Myr of ape evolution. Hartley et al.~2026
@acrocentric_rdna_robertsonian_hartley2026biobank developed a short-read
genotyping method for ROB carriers and applied it to UK Biobank (\~490 k
samples) plus a healthy newborn cohort, recovering the historical
0.11–0.12% prevalence and uncovering hidden distal-junction copy-number
variation in HPRC genomes.

For the present manuscript these papers are the substrate of #emph[C5’s
acrocentric clade] (chapter 01 community C7) and of the "concerted
evolution" framing in C8: the rDNA-adjacent acrocentric homogenisation
is the most extreme case of within-community paralog distance \< allelic
distance (chapter 04: C7 reverses with paralog closer p=2.0e-7, chapter
12 novel contribution \#18). Concerted evolution is #emph[the] canonical
mode of evolution for the rDNA gene family
@Dover1982@Smith1976@Liao1999@Hurles2004; the present manuscript extends
that mode beyond rDNA into the surrounding subtelomeric blocks.
Critically, de Lima 2025’s identification of the SST1 macrosatellite as
the recombining homology supplies the #emph[molecular] mechanism for the
dominant Robertsonian products — direct evidence that ongoing exchange
among the five short arms is mediated by specific homologous templates,
not by generic NOR-clustering alone.

=== D4Z4, FSHD biology and the chr4q–chr10q axis (topic 06)
<d4z4-fshd-biology-and-the-chr4qchr10q-axis-topic-06>
The D4Z4 macrosatellite repeat at the 4q35 subtelomere is the second
canonical example. Wijmenga et al.~1992 mapped 4q-tip rearrangements to
facioscapulohumeral muscular dystrophy @Wijmenga1992; van Deutekom et
al.~1996 @vanDeutekom1996 showed Southern-blot evidence of subtelomeric
exchange of 3.3 kb tandem D4Z4 units between chr4q35 and chr10q26 in
roughly 20% of typed individuals; Lemmers et al.~2002 @Lemmers2002
established that FSHD is associated with one of two 4q subtelomere
variants (4qA / 4qB); Masny et al.~2004 @Masny2004 localised 4q35.2 to
the nuclear periphery via lamin A/C; Kowaljow et al.~2007 @Kowaljow2007
showed that DUX4 within D4Z4 encodes a pro-apoptotic protein; Ottaviani
et al.~2008 @OttavianiGilson2008 identified a perinuclear positioning
element in subtelomeres that requires A-type lamins and CTCF; Ottaviani
et al.~2009 @Ottaviani2009 then established that D4Z4 itself acts as a
CTCF / A-type lamin-dependent insulator; Snider et al.~2010 @Snider2010
showed FSHD reflects incomplete suppression of a retrotransposed gene;
Lemmers et al.~2010 @Lemmers2010 gave the unifying genetic model —
permissive 4qA haplotype produces stable DUX4 mRNA via a 4q-specific
polyadenylation signal; Cabianca et al.~2012 @Cabianca2012 connected
D4Z4 copy-number to a long ncRNA-mediated polycomb/trithorax epigenetic
switch.

This dense literature is the molecular underpinning of #emph[C5’s
chr4q–chr10q DUX4 community] (chapter 01 community C1, named for the
D4Z4 sharing) and #emph[C7’s nuclear-envelope hypothesis] (the
abstract’s "facilitated by the physical proximity of subtelomeres at the
nuclear envelope" wording derives directly from the Masny–Ottaviani
D4Z4-CTCF-lamin axis). The present manuscript’s chapter 03 confirms 28
DUX4L genes in C1, chapter 07 cites the Masny–Ottaviani mechanism in the
integrated synthesis, and chapter 14 detects a chr10q ← chr4q
gene-conversion event at score 0.957 in PAN028 maternal — a
#emph[directly observed] subtelomeric exchange event in the WashU
pedigree confirming the FSHD-context translocation as ongoing rather
than only inferable from population frequency.

=== Subtelomere population genetics and out-of-Africa structure (topic 13)
<subtelomere-population-genetics-and-out-of-africa-structure-topic-13>
The pseudohomology concept’s third population-genetics dimension is
whether subtelomeric block frequencies carry phylogenetic signal.
Mefford & Trask 2002 raised this directly, hedging that "measuring the
frequency of specific subtelomeric blocks on particular chromosomes
might inform phylogenetic studies of modern humans" but cautioning that
"such frequencies could be an unreliable indicator of the relationships
of human populations" because they reflect exchange history as well as
population history @MeffordTrask2002. Der-Sarkissian et al.~2002
@DerSarkissian2002 supplied early evidence that proterminal markers vary
across populations.

The decade of population-genomics that followed (1000 Genomes Project
Consortium 2015 @ThousandGenomes2015; Sudmant et al.~2015 SV map
@Sudmant2015; Mallick et al.~2016 SGDP @Mallick2016; Bergström et
al.~2020 HGDP @Bergstrom2020) made it possible to test Mefford’s hedge:
do subtelomeric block frequencies, treated as alleles, recover the
canonical out-of-Africa topology? Andrea’s chapter 04 (Hudson Fst) and
chapter 12 (novel contribution \#19) report that they do. From cross-arm
affinity Fst at 9 arms × 5 superpopulations, AFR is the deepest split
(distance 0.73–0.81), AMR–EUR are closest (distance 0.22), the topology
mirrors out-of-Africa, and AFR/non-AFR Fst is 0.10–0.15 while non-AFR
pairs are 0.00–0.02. This is the direct quantitative answer to Mefford &
Trask’s 2002 hedge: subtelomeric blocks #emph[do] carry phylogenetic
signal at population scale, and the signal is large enough to recover
deep splits. It supports C6 (PCA + community detection on similarity
matrix → subtelomere clustering across human populations) and provides
one of the strongest novel contributions of the manuscript — chapter 12
novel contribution \#19.

=== Olfactory-receptor pseudogenisation and the f7501 / OR4F gradient (topic 14)
<olfactory-receptor-pseudogenisation-and-the-f7501-or4f-gradient-topic-14>
Olfactory receptors are the largest mammalian gene family and the
canonical example of a duplicating, pseudogenising gene family at
chromosome ends. Glusman et al.~2001 catalogued the complete human
olfactory subgenome @Glusman2001; Niimura & Nei 2003 documented
evolutionary dynamics @Niimura2003; Olender et al.~2008 published a
structured update @Olender2008; Hasin-Brumshtein, Lancet and Olender
2009 connected genomic variation to phenotypic diversity @Hasin2009;
Niimura, Matsui and Touhara 2014 placed the human OR repertoire in a
13-mammal comparative frame @Niimura2014. Trask 1998 @Trask1998 is the
founding empirical observation that OR cluster blocks are duplicated
polymorphically near chromosome ends in humans; cosmid f7501 (L78442)
was the field’s reference probe for one such block.

The Mefford & Trask 2002 review @MeffordTrask2002 documented f7501
distribution across \~52 individuals at five chromosome ends (chr3q,
chr15q, chr19p, chr16q, chr7p) and found chr16q enriched in African
populations. Andrea’s chapter 01 §"f7501 (L78442) direct alignment
validation" reproduces this #emph[quantitatively] at population scale:
chr3q (91.8%), chr19p (90.5%), chr15q (85.6%) confirmed #emph[fixed]
across the 466-haplotype dataset; chr16q AFR-enriched at p=6.6e-27 (70%
AFR carriers vs \<10% non-AFR); novel SAS-enriched chr2q (p=6.8e-11),
AMR-enriched chr6p (p=7.0e-04), and three additional AFR-enriched arms
(chr8p, chr16p, chr9q). Chapter 03 then reports that 62.1% of 5,023 OR4F
annotations are pseudogenes, with pseudogenisation rate varying from
11.1% (chr7p) to 99.8% (chr15q) — the "OR4F gradient" finding in the C5
community structure (chapter 12 novel contribution implicit in the
OR4F→C8 / DUX4L→C1 / IQSEC3→C5 / IL9RP1→C3 mapping table from chapter
03).

For the present manuscript, the OR-family literature is critical context
for two reasons. First, OR pseudogenisation gradients across communities
give the C5 cladistic structure a #emph[mechanistic] interpretation
beyond pure sequence affinity: communities differ in the rate at which
OR genes lose function as they move across chromosome ends, suggesting
different selective regimes at different chromosome ends. Second, the
f7501 reproduction in chapter 01 is the single cleanest "old discovery
confirmed at scale" result in the manuscript and should be foregrounded
in any BoG Q&A on the validity of community-detection methods.

=== Summary — pseudohomology, four ways
<summary-pseudohomology-four-ways>
By the end of Part II the field knows four kinds of pseudohomology. (i)
PAR-style obligate homology — chrX/Y PAR1 and PAR2, where exchange is
functionally required and the gradient of sex linkage
@Rouyer1986@sexchrompars_goodfellow1986@sexchrompars_freije1992@Filatov2004@Charchar2003@sexchrompars_hinch2014@sexchrompars_bellott2024
is the canonical signature, with recurrent population-scale length
polymorphism establishing that PAR1 is a recombination #emph[system]
rather than a fixed interval
@sexchrompars_mensah2014@sexchrompars_poriswanish2018@sexchrompars_poriswanish2025;
this is Part I material but is the cleanest prototype. (ii) Acrocentric
rDNA-flanking concerted evolution — the chr13/14/15/21/22 short arms,
where homogenisation extends across the rDNA and into the surrounding
subtelomere; canonical predecessor papers
@Henderson1972@Floutsakou2013@McStay2016@vanSluis2019@Guarracino2023@deLima2025;
the present manuscript’s C7 community. (iii) Block-shared dynamic
patchworks — the FISH-era / Linardopoulou / Ambrosini concept of
recurrent ectopic exchange among subtelomeric blocks at non-acrocentric
ends @MeffordTrask2002@Linardopoulou2005@Ambrosini2007; this is what
most of the present manuscript’s communities (C1, C2, C3, C5, C6)
represent. (iv) Disease-associated subtelomeric duplicons — D4Z4 / FSHD
at chr4q-chr10q @Wijmenga1992@vanDeutekom1996@Lemmers2010 — a single
named case where pseudohomology has clinical weight. The C5 cladistic
claim of the abstract names instances of all four kinds; the C8
"concerted evolution" framing is the umbrella across all four. Part III
now asks the mechanistic question: what is the cellular and molecular
machinery that creates and maintains this pseudohomology?

#line(length: 100%)

== Part III — Mechanism: meiotic 3D and recombination (combines topics 07, 08, 09)
<part-iii-mechanism-meiotic-3d-and-recombination-combines-topics-07-08-09>
=== Concerted evolution, gene conversion and NAHR (topic 07 — refreshed v3-opus pass)
<concerted-evolution-gene-conversion-and-nahr-topic-07-refreshed-v3-opus-pass>
The concerted-evolution concept enters the field with Smith 1976
@Smith1976, who showed that unequal crossover among tandem repeats
generates a homogenisation dynamic, and Brown 1972 @Brown1972 who had
already demonstrated the empirical pattern in #emph[Xenopus] rDNA. Dover
1982 @Dover1982@Dover1986 gave the field the #emph[molecular drive]
framework — a cohesive mode of species evolution in which intragenomic
processes (gene conversion, unequal crossover, transposition) homogenise
repetitive sequence families across the genome. Liao 1999 @Liao1999
consolidated the concerted-evolution literature and articulated unequal
crossover and gene conversion as the two principal molecular engines;
Hurles 2004 @Hurles2004 supplied the CMT1A paralogous-repeat
gene-conversion exemplar. Nei and Rooney 2005 @NeiRooney2005 introduced
the #emph[birth-and-death] alternative for some multigene families — an
important caveat for OR4F (where the family does turn over rather than
being homogenised) and a reason the abstract uses "concerted evolution"
only in the loose sense.

The NAHR / genomic-disorders frame entered with Lupski 1998 @Lupski1998
and was canonised by Stankiewicz & Lupski 2002 / 2010
@StankiewiczLupski2002@StankiewiczLupski2010, with the broader
genomic-disorders synthesis built up across Lupski & Stankiewicz 2005
@LupskiStankiewicz2005 and Lupski 2009 @Lupski2009. Bailey & Eichler
2006 @BaileyEichler2006 argued that primate segmental duplications are
evolutionarily young, structurally biased toward pericentromeric and
subtelomeric locations, and recombinationally active — directly relevant
to subtelomeres as a high-NAHR class. Chen, Cooper and colleagues 2007
@ChenCooper2007 reviewed gene conversion mechanisms and human disease;
Hastings et al.~2009 @Hastings2009 integrated NAHR, NHEJ and
FoSTeS/MMBIR into a single mechanism-of-CNV-change synthesis. Lovett
2004 @Lovett2004 supplied the misalignment / template-switching
mechanism at the nucleotide level. Bailey 2002 @Bailey2002 and Sharp
2006 @Sharp2006 established that recent segmental duplications are the
genome-wide arena for these mechanisms, and Sudmant 2015 @Sudmant2015 /
Vollger 2023 @Vollger2023 placed segmental duplication in the
population-scale and T2T-scale context.

The third mechanism — break-induced replication (BIR) — joined the
inventory through Llorente, Smith and Symington 2008 @Llorente2008 and
Anand, Lovett and Haber 2013 @Anand2013, which showed that long-tract
synthesis from a single homologous template can copy hundreds of
kilobases without crossover. BIR is the operationally relevant mechanism
for the kilobase-scale gene-conversion patches the present manuscript
reports: the WashU pedigree’s 538 patches range from a few kb up to
hundreds of kb, and 133 are gene-conversion-like (no crossover) at
perfect 1.000/1.000 score — a size distribution and exchange morphology
more compatible with BIR than with classical short-tract gene
conversion.

The intra-array biology of rDNA — the strict-sense concerted-evolution
case — was quantified by Eickbush & Eickbush 2007 @EickbushEickbush2007
and Ganley & Kobayashi 2007 @GanleyKobayashi2007: total repeat
homogenisation in real time, mechanism mapped onto recombination at the
ribosomal locus. Galtier et al.~2001 @Galtier2001 and Duret & Galtier
2009 @DuretGaltier2009 established that gene conversion is GC-biased,
with population-genetic and isochore consequences; gBGC is a candidate
explanation for why subtelomeric communities tend to retain GC-rich
block structure. Jeffreys & May 2004 @JeffreysMay2004 measured human
gene-conversion tract length at \~50–500 bp from sperm crossover
hotspots; Williams et al.~2015 @Williams2015 and Halldorsson et al.~2019
@Halldorsson2019 produced direct pedigree-resolved gene-conversion calls
at sequence resolution, providing the closest external comparator for
the WashU pedigree analysis.

For the present manuscript, two threads matter most. First, NAHR is the
canonical mechanism for #emph[crossover-like] PHR exchange events
(chapter 14 reports 16 such events in the WashU pedigree), and gene
conversion (with BIR as a long-tract variant) is the mechanism for the
#emph[gene-conversion-like] events (chapter 14: 133 events, 96 of which
sit at perfect 1.000/1.000 score). Second, "concerted evolution" in the
abstract title refers — per the lead author’s clarification — to the
#emph[loose] sense (ongoing recombination exchange producing
sequence-similarity homogenisation across non-homologous arms), not the
strict molecular-evolution sense (a coherent, often co-evolutionary
process across a gene family). Nei & Rooney 2005 @NeiRooney2005 is the
explicit reason for this loose-sense framing: not all subtelomeric
paralog families behave like rDNA, and the loose sense covers both
homogenising (rDNA) and turnover-prone (OR4F) systems. The pedigree
result (Part IV) is what upgrades that loose sense from inference to
direct observation.

=== The meiotic bouquet and nuclear-envelope tethering (topic 08)
<the-meiotic-bouquet-and-nuclear-envelope-tethering-topic-08>
The meiotic bouquet is the central cell-biological observation that ties
subtelomeric sequence homology to chromosome positioning. Scherthan,
Weich and colleagues 1996 @Scherthan1996 showed that centromere and
telomere movements during early meiotic prophase of mouse and man are
associated with the onset of chromosome pairing. Zickler & Kleckner’s
1998 review @ZicklerKleckner1998 consolidated the leptotene-zygotene
transition as the key window; their 2015 #emph[Cold Spring Harbor
Perspectives] article @ZicklerKleckner2015 is the canonical modern
review of recombination, pairing and synapsis with the LINC complex /
SUN-KASH machinery framing.

The molecular machinery is now well-characterised. Sosa, Rothballer and
colleagues 2012 @Sosa2012 solved the structural basis of LINC complex
formation through KASH-peptide / SUN-trimer binding; Penkner et al.~2009
@Penkner2009 identified the SUN-1 modification cycle that drives meiotic
chromosome homology search through the nuclear envelope; Hiraoka &
Dernburg 2009 @Hiraoka2009 reviewed the SUN-protein meiotic chromosome
dynamics; Chikashige et al.~2006 @Chikashige2006 identified the
fission-yeast Bqt1/Bqt2 telomere tethers; Lee et al.~2015 @Lee2015
worked out rapid mouse telomere prophase movements; Boateng, Bellani and
Camerini-Otero 2013 @Boateng2013 documented homologous pairing preceding
SPO11-mediated double-strand breaks in mice. Burgoyne 1982 @Burgoyne1982
is the original X/Y obligate-crossover argument that anchors the C7
hypothesis. The mouse-PAR mechanism extends this further: Acquaviva,
Boekhout, Karasu et al.~2020 @sexchrompars_acquaviva2020 showed that
ANKRD31, MEI4 and axis-elongation factors actively #emph[force] DSB
formation into the small pseudoautosomal interval to ensure the obligate
X-Y crossover, and Chu, Froberg, Kesner et al.~2017
@sexchrompars_chu2017 identified PAR-TERRA RNA as a pairing factor for
sex-chromosome homologous regions in mouse ESCs. Together these papers
establish that chromosome-end pairing is #emph[actively maintained] by
sequence-aware molecular machinery, not merely a passive consequence of
telomere clustering.

Two takeaways matter for C7 and C8. First, the meiotic bouquet places
telomere-adjacent regions at the nuclear envelope, in physical proximity
to one another, during a window (leptotene→zygotene) in which homology
search and double-strand-break repair are most active. This is the
biophysical basis for the C7 "facilitated by the physical proximity of
subtelomeres at the nuclear envelope" wording. Second, sequence
similarity at chromosome ends biases the homology search toward ectopic
(interchromosomal) rather than purely allelic (intrachromosomal
homologue-pair) repair, which is the proposed feedback-loop mechanism in
C8 (chapter 07 §"Causal feedback loop"; chapter 11 finding 11; chapter
12 novel contribution \#11): sequence similarity → 3D proximity at the
bouquet → ectopic exchange → increased similarity. Bandyopadhyay et
al.~2002’s nucleolar-association-as-cause argument for Robertsonian
translocation @acrocentric_rdna_robertsonian_bandyopadhyay2002 is the
corresponding mechanism for the acrocentric C7 community: the five
p-arms are co-located at a #emph[different] nuclear landmark (the
nucleolus rather than the lamina), and the spatial mechanism therefore
differs from the D4Z4 / lamin tether but the principle is identical.

=== Hi-C, Pore-C, single-cell and meiotic 3D (topic 09)
<hi-c-pore-c-single-cell-and-meiotic-3d-topic-09>
Lieberman-Aiden et al.~2009 @LiebermanAiden2009 introduced Hi-C and
demonstrated bulk principles of human genome folding (compartments,
contact frequency vs distance, principal genome partitions); Rao et
al.~2014 @Rao2014 extended this to kilobase resolution and identified
looping principles. Single-cell variants followed: Nagano et al.~2013
@Nagano2013 introduced single-cell Hi-C; Stevens et al.~2017
@Stevens2017 resolved 3D structures of individual mammalian genomes;
Ramani et al.~2017 @Ramani2017 built massively multiplex single-cell
Hi-C; Tan, Xing and colleagues 2018 introduced Dip-C @Tan2018 and
resolved GM12878 and PBMC single-cell 3D structures; Ulahannan et
al.~2019 @Ulahannan2019 gave Pore-C the long-read concatemer framing.

Meiotic 3D came of age with Patel et al.~2019 @Patel2019 (mouse meiotic
Hi-C, chromosome-end clustering during meiotic prophase) and especially
Zuo et al.~2021 @Zuo2021 (stage-resolved mouse meiotic Hi-C across
leptotene, zygotene, pachytene, diplotene; LINC complex / SUN1 force
transmission; meiotic chromosome-end alignment extends roughly 20% of
chromosome length; the SUN1 W151R mutant compresses this to roughly 5%).
Xu et al.~2025 @Xu2025 extended single-cell 3D into 20 sperm cells. For
the present manuscript these papers supply the analytical backbone of
C7: Mantel tests of sequence-similarity × Hi-C contact at six bulk Hi-C
samples, two Pore-C and CiFi datasets, GM12878 and sperm Dip-C, and
Zuo’s mouse meiotic Hi-C all converge on rho ≈ 0.66–0.72 (chapter 05;
chapter 06; chapter 08). The mouse zygotene rho=0.715 (Zuo 2021
stage-resolved Hi-C) is the strongest existing evidence that the
sequence-similarity ↔ 3D-proximity correspondence peaks at the meiotic
stage where the bouquet is most active.

Methodologically, three subsidiary innovations are also notable: (i)
Stergachis et al.~2020 @Stergachis2020 introduced Fiber-seq for
single-molecule chromatin-fibre regulatory architectures; (ii) Gershman
et al.~2022 @Gershman2022 re-aligned ENCODE CTCF ChIP-seq to T2T-CHM13,
confirming CTCF enrichment at TAR loci; (iii) Tan 2018’s Dip-C produced
per-haplotype radial-position data that the present chapter 05 §"Nuclear
lamina cross-reference" uses as an envelope proxy (C1 D4Z4 radial=0.732
peripheral; C14 PAR2 0.840). These are the analytical bricks that allow
chapter 05’s "the Mantel rho strengthens for ALL 8 datasets at ALL
resolutions when acrocentric / sex / strongest-community exclusions are
applied" — the cleanest available defence against confound-driven Hi-C
signal in subtelomeres.

=== Synthesis of mechanism
<synthesis-of-mechanism>
Three mechanistic threads converge. First, #emph[concerted-evolution
machinery] (gene conversion, NAHR, BIR, segmental-duplication-mediated
rearrangement) supplies the molecular events that produce sub-block
similarity, with the strict-sense rDNA case
@EickbushEickbush2007@GanleyKobayashi2007 as the asymptote of fast
intra-array homogenisation and the loose-sense subtelomere-block case as
the slower inter-arm version. Second, the #emph[meiotic bouquet]
supplies a cell-biological window in which interchromosomal homology
search at chromosome ends is maximal, biasing repair toward ectopic
outcomes when sequence similarity exists; nucleolar association supplies
an analogous window for the acrocentric short arms. Third, #emph[Hi-C
and single-cell 3D] methods supply the empirical access point: 3D
contact across non-homologous chromosome ends should track sequence
similarity if the bouquet hypothesis is correct, and the present
manuscript reports exactly that correspondence at multiple scales. Part
IV now asks how the methods (pangenomes, pedigrees) make this story
testable at population scale.

#line(length: 100%)

== Part IV — Methods: pangenomes and pedigrees (combines topics 10, 11, 12)
<part-iv-methods-pangenomes-and-pedigrees-combines-topics-10-11-12>
=== Pangenome graphs and IMPG transitive closure (topic 10)
<pangenome-graphs-and-impg-transitive-closure-topic-10>
The pangenome graph framework emerged from the recognition that
single-reference workflows under-represent diverse alleles. The
Computational Pan-Genomics Consortium 2018 @Computational2018
consolidated the early conceptual frame; Eizenga et al.~2020
@Eizenga2020 published the canonical "Pangenome graphs" review; Garrison
et al.~2018 introduced vg @Garrison2018; Li, Feng & Chu 2020 introduced
minigraph @Li2020minigraph for reference-based graph construction;
Garrison et al.~2024 introduced pggb (pangenome graph builder)
@Garrison2024pggb for reference-free all-vs-all graph construction;
Heumos et al.~2024 published the odgi infrastructure used to compute
pairwise Jaccard similarity on the pggb GFA @Heumos2024. The
Computational Pan-Genomics Consortium 2018 paper @Computational2018
(Marschall, Marz, Abeel et al.) is the consortium-level review.

For the present manuscript the canonical pipeline is: wfmash all-vs-all
pairwise mapping at p95 identity over the 18,827 telomere-anchored 500
kb flanks → IMPG (https:\/\/github.com/pangenome/impg) transitive
closure on the implicit pangenome graph → pggb GFA construction →
odgi-derived pairwise Jaccard matrix → Leiden community detection at
arm-level (k=15) and sequence-level (k=50). The conceptual point that
the topic\_10 review establishes is that the all-vs-all PAF set
#emph[is] the implicit pangenome graph: each PAF edge is a pairwise
mapping, the union over flanks is the graph, and IMPG `query -x`
performs transitive closure on it. The pggb-built explicit GFA is one
downstream product, not the primary object. Per Erik’s clarification
(CROSSWALK §7a), the Methods should cite IMPG and reframe the all-vs-all
PAF set as the implicit graph; the field already accepts the concept and
the Methods does not need to invent it. The Erdős-Rényi connectivity
argument (CROSSWALK §7b: n=18,827 → p\* = log(n)/n ≈ 5.21e-4; 12%
sampling is \~230× above threshold; the resulting random graph is
densely connected and transitive closure reaches everywhere) is the
formal justification for "no chromosomal partitioning."

Three canonical companion papers extend the pangenome frame. Guarracino
et al.~2023 @Guarracino2023 is the closest predecessor: pangenome-graph
evidence of recombination between heterologous human acrocentric
chromosomes, using HPRC v1 data. de Lima et al.~2025 @deLima2025 is a
follow-on pangenome-resolved analysis of human Robertsonian breakpoints.
Rautiainen et al.~2023 @Rautiainen2023 introduced Verkko, the diploid
T2T assembly pipeline whose output feeds the WashU pedigree analysis.
These three together with Liao 2023 @Liao2023 (HPRC v1) define the
pangenome-method branch of the lineage that the present manuscript
continues.

=== Pedigree-based variant and recombination detection (topic 11)
<pedigree-based-variant-and-recombination-detection-topic-11>
Topic 11 supplies the second methodological branch. The committed
topic\_11 review
\[paper\_prep/lit\_review/topic\_11\_pedigree\_based\_recombination\_detection.md\]
documents this in detail. Briefly: deCODE’s fine-scale pedigree maps
@Kong2010@Bherer2017 and Halldorsson et al.’s sequence-level genetic map
@Halldorsson2019 established the logic of using transmitted alleles to
infer recombination and de novo mutation; Sasani et al.~2019’s
three-generation Utah families @Sasani2019 showed how grandchildren
resolve postzygotic vs germline mosaicism; Chaisson et al.~2019
@Chaisson2019, Audano et al.~2019 @Audano2019, Collins et al.~2020
@Collins2020 and Ebert et al.~2021 @Ebert2021 showed that long-read or
multi-platform structural-variation discovery exposed many variants
invisible to short reads; Cheng et al.~2021 (hifiasm) @Cheng2021 and
Rautiainen et al.~2023 (Verkko) @Rautiainen2023 supplied
haplotype-resolved assembly methods; Wagner et al.~2022 @Wagner2022 and
Zook et al.~2020 @Zook2020 gave benchmarking frameworks for difficult
variants and SVs; Ahsan et al.~2023 @Ahsan2023 and Liu et al.~2024
@Liu2024 reviewed long-read SV detection algorithm tradeoffs; Nurk et
al.~2022 @Nurk2022 and Liao et al.~2023 @Liao2023 supplied the T2T-CHM13
reference and HPRC v1 pangenome.

The two papers that are most directly load-bearing for C8 are Cechova et
al.~2025 @Cechova2025 and Porubsky et al.~2025 @Porubsky2025. Cechova et
al.~2025 supplies the WashU 3-generation T2T pedigree (PAN010/PAN011 →
PAN027 → PAN028), the primary substrate for the abstract’s "538 patches
of recent recombination, 92% confined within sequence-similarity
communities" claim. Porubsky et al.~2025 supplies the CEPH1463
4-generation pedigree (28 hifiasm + 14 verkko phased assemblies), used
for cross-assembler validation: 11 robust inter-chromosomal subtelomeric
features that replicate across both pipelines (chr10/chr18 in
NA12877+NA12878 = the Linardopoulou 2005 pair; chr19/chr22 transmitted
via NA12878; chr12/chr9 in both NA12889 and NA12890). The methodological
principle the topic\_11 review documents most carefully — and that the
manuscript’s chapter 14 inherits — is that an inter-chromosomal patch in
a child is evidentiarily #emph[strongest] when it (a) sits on an
inherited parental haplotype with the expected Mendelian phase, (b) is
reproducible across independent assemblers or generations, and (c) sits
inside a Leiden community independently defined from unrelated HPRC v2
haplotypes.

=== HPRC v1 / v2 population pangenomes (topic 12)
<hprc-v1-v2-population-pangenomes-topic-12>
The HPRC literature is the dataset substrate. Nurk 2022 @Nurk2022
(T2T-CHM13) closed the reference-genome side; Liao 2023 @Liao2023 is the
HPRC v1 draft pangenome (47 individuals, 94 haplotypes); HPRC v2
(forthcoming, preprinting per lead-author clarification) extends this to
233 individuals × 2 haplotypes = 465 plus CHM13 = 466
haplotype-equivalent units in the present manuscript. The technical
story is straightforward: hifiasm @Cheng2021 and Verkko @Rautiainen2023
deliver haplotype-resolved assemblies; pggb @Garrison2024pggb / odgi
@Heumos2024 / IMPG (in submission) deliver the pangenome-graph
operations; Wagner 2022 @Wagner2022 / Zook 2020 @Zook2020 / GIAB-style
benchmarking discipline keeps the dataset rigorous. The dataset is large
enough that subtelomeric block frequencies become measurable as
population variables; this is the single technical advance that makes C5
(cladistic structure with named clades) and C6 (population clustering
across human ancestries) possible.

The methodological hand-off from topic 10 → topic 11 → topic 12 is the
architecture of the present paper’s Methods: pangenome graph ingests 466
haplotypes; community detection partitions them; pedigree analysis tests
whether the partition predicts inheritance; population-genetics analysis
tests whether the partition tracks ancestry. Each branch supplies a
different #emph[kind] of validation. The pangenome graph supplies
#emph[internal] structure (sequence-similarity communities). The
pedigree supplies #emph[transmission] validation (do new generations
inherit and recombine within those communities?). The
population-genetics analysis supplies #emph[external] validation (do
those communities mirror known population structure?). These three
together are the empirical core of C5–C8.

#line(length: 100%)

== Part V — How this paper sits in the field
<part-v-how-this-paper-sits-in-the-field>
Part V is the synthesis the manuscript itself depends on: which of the
abstract’s eight claims are old discoveries reproduced at scale, which
are new, and which depend on framing decisions that the abstract has
chosen to make explicit. This section follows the C-numbering of
`paper_prep/synthesis/CROSSWALK.md`.

=== What is recovered (pre-existing literature confirmed at population scale)
<what-is-recovered-pre-existing-literature-confirmed-at-population-scale>
C1 (HPRC v2 companion-paper framing) is editorial rather than empirical;
it commits the present manuscript to the HPRC v2 publication arc and is
correct as soon as the HPRC v2 main paper is in press. C2 (implicit
pangenome graph; \~12% pairwise sampling) reframes the standard wfmash +
pggb + odgi pipeline as a graph object and supplies an Erdős-Rényi
connectivity argument; the #emph[concept] is established in
@Eizenga2020@Computational2018 and the #emph[implementation] is
established in @Garrison2018@Garrison2024pggb@Heumos2024. C3 (466
near-complete haplotype assemblies) is direct continuation of @Liao2023
using HPRC v2; the dataset itself is the new contribution.

C4 (extended interchromosomal homology at nearly all subtelomeres
comparable to PAR2) recovers and dramatically extends Flint, Bates et
al.~1997 @Flint1997 and Riethman et al.~2004 @Riethman2004. The
two-domain model is confirmed on 39/48 arms by gradient analysis and
39/41 by piecewise breakpoint analysis (chapter 04 §"Two-domain
subtelomeric model test"); internal (TTAGGG)n islands co-localise with
the proximal-domain boundary within 25 kb on 11/19 testable arms
(chapter 04). The "comparable to PAR2" wording is a #emph[new] explicit
comparison; the Flint era did not run this comparison directly.

C5 (cladistic analysis with named clades) recovers Linardopoulou 2005
@Linardopoulou2005 (chr10p–chr18p), Mefford & Trask 2002
@MeffordTrask2002 (f7501-cluster, chr4q–chr10q D4Z4), Lemmers 2010
@Lemmers2010 (chr4qA / DUX4 polymorphism), and Trask 1998 @Trask1998 (OR
cluster polymorphism), and supplies a new {22q, 21q, 19q, 1q, 13q, 17q}
clade plus a large moderate-similarity clade. Note framing: the abstract
names "neighbour-joining trees" but Andrea’s analysis runs Leiden +
UPGMA with 12-of-15-community agreement (chapter 01); the named clades
all map to communities, so this is editorial rather than analytical. The
novel contribution is #emph[quantitative] (population-scale community
structure of 41 arms across 466 haplotypes — chapter 12 novel
contribution \#1) and #emph[cladistic resolution] (named clades that the
reference-only era could not assign confidently).

C6 (PCA / community detection on similarity matrix → subtelomere
clustering across human populations) recovers Mefford & Trask 2002’s
hedge about whether subtelomeric block frequencies could be phylogenetic
markers. The answer is positive: the out-of-Africa Fst tree from chapter
04 / chapter 12 novel contribution \#19 recovers AFR as the deepest
split with mean Fst 0.044, AFR/non-AFR Fst 0.10–0.15 vs non-AFR
0.00–0.02; Mefford’s hedge is now a confirmed positive result.

C7 (Hi-C 3D maps testing nuclear-envelope-proximity recombination
hypothesis) recovers Stout et al.~1999 @Stout1999 (interphase
chr4q/chr10q clustering by FISH), Masny 2004 @Masny2004 (4q35.2
peripheral lamin-A/C tethering), Ottaviani 2009 @Ottaviani2009
(D4Z4-CTCF-lamin), Zuo 2021 @Zuo2021 (chromosome-end alignment extends
\~20% in mouse meiosis), Patel 2019 @Patel2019 (mouse meiotic Hi-C),
Acquaviva 2020 @sexchrompars_acquaviva2020 (PAR-axis-elongation force on
DSB localisation in mouse), and Chu 2017 @sexchrompars_chu2017
(PAR-TERRA-mediated sex-chromosome pairing). The novel contribution is
#emph[direct quantitative testing] in human bulk and single-cell Hi-C
and Pore-C: Mantel rho ≈ 0.66 across 8 datasets, strengthening to 0.79
after exclusion controls (chapter 05). The "flanking-paradox" finding
(chapter 05 / chapter 07 / chapter 11 finding 9) — unique-sequence
flanking shows #emph[stronger] 3D enrichment than duplicated PHR — rules
out multi-mapping artifacts. The framing gap noted in CROSSWALK §6 is
that no direct LAD / Lamin B1 ChIP overlay exists yet; Dip-C radial
position is used as a proxy. REWRITE\_PLAN TASK-22 is the appropriate
follow-up; if unavailable, the abstract should soften "nuclear envelope"
to "peripheral telomere positioning."

C8 (synthesising thesis: "concerted evolution and unorthodox
recombination") recovers Dover 1982 @Dover1982 / Liao 1999 @Liao1999
(concerted evolution as molecular drive), Stankiewicz & Lupski
@StankiewiczLupski2002 (NAHR), Henderson 1972 @Henderson1972 /
Floutsakou 2013 @Floutsakou2013 / McStay 2016 @McStay2016 (rDNA-flanking
concerted evolution), Linardopoulou 2005 @Linardopoulou2005
(subtelomeres as hot spots of interchromosomal recombination), and
Guarracino 2023 @Guarracino2023 (pangenome-graph evidence for
inter-acrocentric recombination); the v3-opus topic\_05 refresh
explicitly adds de Lima 2025 @deLima2025 (SST1 macrosatellite as the
recombining homology for dominant Robertsonian fusions) and Stults 2008
@acrocentric_Stults2008 (\>10% per-meiosis rDNA-cluster rearrangement
frequency) as direct mechanistic anchors for the loose-sense usage. The
single most-novel contribution is #strong[direct pedigree evidence];:
538 high-quality interchromosomal patches in the WashU T2T pedigree, 92%
within Leiden communities, including 133 gene-conversion-like events at
perfect 1.000/1.000 score, 16 crossover-like events, 229 acros-like
events; the chr10q ← chr4q gene-conversion event observed at score 0.957
in PAN028 maternal. Per the lead-author clarification (CROSSWALK §5),
the pedigree IS the proof of concerted evolution in the loose sense.

=== What is genuinely new
<what-is-genuinely-new>
A short list, in priority order:

+ #strong[Direct pedigree evidence of ongoing inter-chromosomal
  exchange] (chapter 14). 538 patches, 92% within communities, 133
  gene-conversion-like at perfect score, 16 crossover-like, transmitted
  across three generations in WashU T2T plus 11
  cross-assembler-validated CEPH1463 features. This is the most-novel
  contribution and the empirical core of C8.

+ #strong[Population-scale community structure] (chapter 01). 41 arms in
  15 arm-level Leiden communities and 50 sequence-level communities
  across 466 haplotypes, the first population-scale quantification of
  the patchwork model.

+ #strong[Two-domain Flint/Mefford model confirmed at pangenome scale]
  (chapter 04). 39/48 arms gradient, 39/41 piecewise, 11/19 ITS
  co-localised, 99.7% per-haplotype. Extends Flint 1997 by a factor of
  \~13.

+ #strong[Subtelomeric Fst mirrors out-of-Africa] (chapter 04 / chapter
  12 novel contribution \#19). Subtelomeric block frequencies do carry
  phylogenetic signal; Mefford & Trask’s 2002 hedge is now a positive
  result.

+ #strong[3D-mirrors-sequence at multiple scales and resolutions]
  (chapter 05 / chapter 06 / chapter 07 / chapter 08). Mantel rho ≈
  0.66–0.79 across 8 bulk Hi-C / Pore-C / CiFi / Dip-C / sperm / mouse
  meiotic Hi-C datasets, with mouse zygotene rho=0.715.

+ #strong[Allele-vs-paralog Wilcoxon at population scale] (chapter 04).
  5,946 paired distances, 8/9 multi-arm communities allele-closer
  (overall p\<1e-300); C7 acrocentric reverses to paralog-closer
  (p=2.0e-7) as quantitative confirmation of complete homogenisation.

+ #strong[Cross-arm type discordance up to 47.5%] at chr22q (chapter
  04). Refines Mefford 2002’s chr4q/chr10q figure from \~20% to 43.4%.

+ #strong[Flanking-paradox] (chapter 05 / chapter 07 / chapter 11
  finding 9). Unique-sequence 100 kb regions centromere-ward of PHR
  boundaries show #emph[stronger] 3D enrichment than the duplicated PHR
  itself — a strong methodological defence against multi-mapping
  artifacts.

+ #strong[Three-category arm classification] (chapter 04 / chapter 11
  finding 2). Homogeneous (8/41 with 0% cross-arm), polymorphic (34/41),
  fully interchangeable (acrocentric p-arms). First quantitative
  framework extending Mefford’s qualitative patchwork model.

+ #strong[Gene-replacement scoring] (chapter 04). Complete (0.91–1.0) at
  chr13\_p/chr14\_p/chr15\_p and PAR; partial (0.0–0.72) at autosomal
  communities. First quantitative gradient of homogenisation onto
  community structure.

+ #strong[OR4F pseudogenisation gradient] (chapter 03). 62.1%
  pseudogenes across 5,023 OR4F annotations; pseudogenisation rate
  11.1%–99.8% across 16 arms.

=== Framing decisions explicit in the abstract
<framing-decisions-explicit-in-the-abstract>
Three framing decisions matter for Q&A:

- #strong["Implicit pangenome graph" terminology and "\~12% pairwise
  sampling"] (C2): the standard pipeline is wfmash all-vs-all + pggb +
  odgi; the abstract reframes the all-vs-all PAF set as the implicit
  graph, citing IMPG, and supplies an Erdős-Rényi argument for why no
  chromosomal partitioning is required. Per CROSSWALK §6 modify-task
  TASK-10 / TASK-11, this is editorial rather than analytical.
- #strong["Neighbour-joining trees"] (C5): the analysis runs Leiden +
  UPGMA, which agree on 12 of 15 communities. The named clades all map;
  the algorithm differs. Either run NJ on the existing arm-distance
  matrix (a one-line `ape::nj()` call) or relax the abstract to
  "cladistic analysis (Leiden on the arm-level Jaccard distance)". Per
  Erik’s clarification, both are legitimate.
- #strong["Concerted evolution"] (C8 title term): used in the
  #emph[loose] sense (ongoing recombination exchange producing
  similarity homogenisation), not the strict molecular-evolution sense.
  Pedigree IS the evidence. The Discussion paragraph should make this
  explicit and reference chapter 14. Nei & Rooney 2005 @NeiRooney2005 is
  the explicit alternative model that justifies the loose-sense framing
  — not all subtelomeric paralog families behave like rDNA.

A fourth, less-discussed framing decision is the #strong["facilitated by
… nuclear envelope"] wording (C7). Andrea’s chapter 05 §"Nuclear lamina
cross-reference" uses Dip-C radial position as a proxy; no genome-wide
LAD / Lamin B1 ChIP overlay exists. The abstract should either be
supplemented with a LAD overlay (REWRITE\_PLAN TASK-22) or softened to
"consistent with peripheral telomere positioning" / "consistent with
envelope association."

=== Closing
<closing>
The Guarracino & Garrison subtelomere manuscript stands at the
convergence of three decades of subtelomeric biology, six layers of
mechanism (concerted-evolution machinery, NAHR, gene conversion, BIR,
segmental duplication, the meiotic bouquet, pangenome representation),
and two new empirical access points (HPRC v2 466-haplotype pangenome
graph; T2T multigenerational pedigrees). Its single most-novel
contribution is the pedigree evidence in chapter 14: subtelomeric
exchange among non-homologous chromosome ends is #emph[currently
happening] in real human families, and it happens predominantly within
sequence-similarity communities that were defined independently from
unrelated population samples. That is the empirical proof of "ongoing
recombination" / "concerted evolution" in the abstract title. Everything
else in the manuscript — the pangenome graph, the cladistic structure,
the Hi-C 3D evidence, the population-genetics signatures — is
corroborating evidence at population scale for what the pedigree
directly observes at family scale. Read in this order, the field history
makes the manuscript’s claims feel inevitable rather than surprising.

The v2 refresh deepens this story in two specific places. Topic 05’s
v3-opus pass adds Krystal 1981, Schmickel 1973, Worton 1988,
Bandyopadhyay 2001/2002, Stults 2008, Stimpson 2014 and the de Lima 2025
SST1 mechanism — converting the acrocentric clade from a recovered
cytogenetic clade to a sequence-resolved, mechanism-anchored community
with a per-meiosis rate (\>10%, Stults 2008) that is the highest
documented in human structural biology. Topic 07’s v3-opus pass adds
Brown 1972, Lupski 1998, Lupski & Stankiewicz 2005, Bailey & Eichler
2006, Eickbush & Eickbush 2007, Ganley & Kobayashi 2007, Llorente 2008,
Hastings 2009, Anand 2013, Williams 2015 and the explicit Nei & Rooney
2005 birth-and-death alternative — converting the C8 mechanism story
from a single-vocabulary NAHR claim to a layered framework (NAHR + gene
conversion + BIR + gBGC, with birth-and-death as the explicit
alternative for OR4F-like families). Both upgrades make the manuscript’s
loose-sense "concerted evolution" usage defensible against a
strict-sense critic.

#line(length: 100%)

#emph[End of SYNTHESIS\_v2.md.]

#pagebreak()
= Appendix A — Chronology of Subtelomere Research

#strong[Author:] lit-review-v2 (agent-937), 2026-05-06. Supersedes
`CHRONOLOGY.md` (lit-review-synthesis, agent-877). #strong[Source:]
distilled from `paper_prep/synthesis/REFERENCES_v3.bib` (295 entries)
and `end-to-end-report/report/12_literature.md`. Each row corresponds to
one published paper. Claims supported (C1–C8) follow
`paper_prep/synthesis/CROSSWALK.md`. Brief description distinguishes
mechanism, dataset, method, prediction, or framework. Entries are sorted
chronologically; ties broken alphabetically by first-author surname.

#strong[Versioning note:] v2 supersedes the v1 chronology and reflects
the final `REFERENCES_v3.bib` merge that incorporates the v3-opus
refresh of topic 05 (acrocentric / rDNA / Robertsonian) and topic 07
(concerted evolution / NAHR / gene conversion).

#figure(
  align(center)[#table(
    columns: (2.47%, 10.7%, 69.14%, 7.41%, 10.29%),
    align: (auto,auto,auto,auto,auto,),
    table.header([Year], [Author(s)], [Brief
      description], [Claim(s)], [Citation key],),
    table.hline(),
    [1916], [Robertson], [Robertsonian acrocentric–acrocentric fusion
    concept (grasshoppers); eponym for human acrocentric short-arm
    rearrangements.], [C5], [Robertson1916],
    [1916], [Robertson], [XII.—Chromosome studies. I. Taxonomic
    relationships shown in the chromosomes of Tettigidae and Acrididae:
    V…], [C5], [cytofound\_Robertson1916],
    [1970], [Caspersson et al.], [Differential binding of alkylating
    fluorochromes in human
    chromosomes], [C5], [cytofound\_Caspersson1970],
    [1970], [Caspersson et al.], [Quinacrine chromosome banding; first
    method to identify all 24 human chromosomes; G-banding
    nomenclature.], [C1, C5], [Caspersson1970],
    [1971], [Caspersson et al.], [The 24 fluorescence patterns of the
    human metaphase chromosomes — distinguishing characters and
    variability], [C5], [cytofound\_Caspersson1971],
    [1972], [Brown et al.], [Xenopus rDNA homogenisation; empirical
    pattern that motivated Smith 1976’s unequal-crossover
    model.], [C8], [Brown1972],
    [1972], [Henderson et al.], [In situ localisation of rDNA arrays to
    all five acrocentric short arms; original cross-chromosome
    concerted-evolution observation.], [C5, C8], [Henderson1972],
    [1972], [Lewontin], [The Apportionment of Human
    Diversity], [C6], [subtel\_popgen\_lewontin1972],
    [1973], [Schmickel], [Quantitation of human ribosomal DNA:
    hybridization of human DNA with ribosomal RNA for quantitation and
    fra…], [C5, C8], [acrocentric\_rdna\_robertsonian\_schmickel1973],
    [1976], [Smith], [Unequal crossover among tandem repeats produces
    homogenisation; mechanistic basis of concerted evolution.], [C7,
    C8], [Smith1976],
    [1978], [Blackburn & Gall], [Original (TTAGGG)n-like canonical
    telomere repeat in Tetrahymena rDNA termini; molecular target of
    TAR1 and ITS biology.], [C4, C5], [BlackburnGall1978],
    [1979], [Hsu], [Standard historical reference for cytogenetic
    foundations (banding, FISH, satellite-DNA
    cytology).], [C5], [Hsu1979],
    [1980], [Arnheim et al.], [Molecular evidence for genetic exchanges
    among ribosomal genes on nonhomologous chromosomes in man and
    apes], [C7, C8], [concerted\_evolution\_nahr\_Arnheim1980],
    [1981], [Krystal et al.], [Human nucleolus organizers on
    nonhomologous chromosomes can share the same ribosomal gene
    variants], [C5, C8], [acrocentric\_rdna\_robertsonian\_krystal1981],
    [1982], [Burgoyne], [Obligate X/Y crossover at PAR1; founding
    statement of pseudoautosomal recombination.], [C4,
    C7], [Burgoyne1982],
    [1982], [Dover], [Concerted evolution / molecular drive theory;
    foundational framework for repetitive-sequence
    homogenisation.], [C8], [Dover1982],
    [1983], [Guichaoua et al.], [Structural basis for Robertson
    translocations in man: association of ribosomal genes in the
    nucleolar fibri…], [C5, C8], [acrocentric\_Guichaoua1983],
    [1984], [Ohta], [Some models of gene conversion for treating the
    evolution of multigene families], [C7,
    C8], [concerted\_evolution\_nahr\_Ohta1984],
    [1984], [Weir & Cockerham], [Estimating F-statistics for the
    analysis of population structure], [C6], [subtel\_popgen\_weir1984],
    [1985], [Greider & Blackburn], [Identification of a specific
    telomere terminal transferase activity in Tetrahymena
    extracts], [C5], [cytofound\_Greider1985],
    [1986], [Dover], [Molecular drive: a third mode of evolution
    alongside selection and drift.], [C8], [Dover1986],
    [1986], [Goodfellow et al.], [A pseudoautosomal gene in man], [C4,
    C7], [sexchrompars\_goodfellow1986],
    [1986], [Pinkel et al.], [Cytogenetic analysis by in situ
    hybridization with fluorescently labeled nucleic acid
    probes], [C5], [cytofound\_Pinkel1986b],
    [1986], [Pinkel et al.], [Quantitative high-sensitivity FISH; method
    enabling Trask-lab subtelomeric polymorphism
    programme.], [C5], [Pinkel1986],
    [1986], [Rouyer et al.], [Gradient of sex linkage in PAR1; original
    empirical pseudoautosomal observation.], [C4, C7], [Rouyer1986],
    [1987], [Pritchard et al.], [Mapping the limits of the human
    pseudoautosomal region and a candidate sequence for the
    male-determining gene], [C4, C7], [sexchrompars\_pritchard1987],
    [1988], [Moyzis et al.], [A highly conserved repetitive DNA
    sequence, (TTAGGG)n, present at the telomeres of human
    chromosomes], [C5], [cytofound\_Moyzis1988],
    [1988], [Worton et al.], [Human ribosomal RNA genes: orientation of
    the tandem array and conservation of the 5’ end], [C5,
    C8], [acrocentric\_Worton1988],
    [1989], [Ellis et al.], [The pseudoautosomal boundary in man is
    defined by an Alu repeat sequence inserted on the Y
    chromosome], [C4, C7], [sexchrompars\_ellis1989],
    [1989], [Greider & Blackburn], [A telomeric sequence in the RNA of
    Tetrahymena telomerase required for telomere repeat
    synthesis], [C5], [cytofound\_GreiderBlackburn1989],
    [1989], [Hastie & Allshire], [Human telomeres: fusion and
    interstitial sites], [C4, C5], [subtelstruct\_HastieAllshire1989],
    [1989], [Morin], [The human telomere terminal transferase enzyme is
    a ribonucleoprotein that synthesizes TTAGGG
    repeats], [C5], [cytofound\_Morin1989],
    [1989], [Therman et al.], [Non-random distribution of human
    acrocentric Robertsonian translocations; chr13–14 and chr14–21
    dominance.], [C5], [Therman1989],
    [1990], [Brown et al.], [Original characterisation of TAR1
    (Telomere-Associated Repeat 1) as a polymorphic, variably arranged
    subtelomeric element.], [C4, C5], [Brown1990],
    [1990], [Ellis et al.], [Population structure of the human
    pseudoautosomal boundary], [C4, C7], [sexchrompars\_ellis1990],
    [1990], [Lichter et al.], [High-resolution FISH using cosmid clones;
    chromosome painting infrastructure for subtelomere
    studies.], [C5], [Lichter1990],
    [1991], [Buck & Axel], [A novel multigene family may encode odorant
    receptors: a molecular basis for odor
    recognition], [—], [BuckAxel1991],
    [1991], [Hillis et al.], [Evidence for biased gene conversion in
    concerted evolution of ribosomal DNA], [C7,
    C8], [concerted\_evolution\_nahr\_Hillis1991],
    [1991], [Ijdo et al.], [Origin of human chromosome 2: an ancestral
    telomere-telomere fusion], [C5], [cytofound\_Ijdo1991],
    [1991], [Trask], [FISH applications in cytogenetics and gene mapping
    (review); Trask-lab toolset for subtelomeric
    polymorphism.], [C5], [Trask1991],
    [1991], [Wilkie et al.], [A truncated human chromosome 16 associated
    with alpha thalassaemia is stabilized by addition of telomeric
    r…], [C4, C5], [subtelstruct\_Wilkie1991],
    [1991], [Wilkie et al.], [Clinical features and molecular analysis
    of the $a l p h a$ thalassemia/mental retardation syndromes. II.
    Cases…], [C5], [cytofound\_Wilkie1991b],
    [1991], [Wilkie et al.], [260 kb stable length polymorphism at
    chr16p tip; argued to arise from non-homologous-end exchange.], [C4,
    C8], [Wilkie1991],
    [1992], [Freije et al.], [Identification of a second pseudoautosomal
    region near the Xq and Yq telomeres], [C4,
    C7], [sexchrompars\_freije1992],
    [1992], [Hudson et al.], [Estimation of levels of gene flow from DNA
    sequence data], [C6], [subtel\_popgen\_hudson1992],
    [1992], [Telenius et al.], [Degenerate oligonucleotide-primed PCR:
    general amplification of target DNA by a single degenerate
    primer], [C5], [cytofound\_Telenius1992],
    [1992], [Wijmenga et al.], [chr4q DNA rearrangements associated with
    FSHD; founding mapping paper for the disease/subtelomere
    link.], [C5, C8], [Wijmenga1992],
    [1994], [Charlesworth et al.], [The evolutionary dynamics of
    repetitive DNA in eukaryotes], [C7,
    C8], [concerted\_evolution\_nahr\_Charlesworth1994],
    [1994], [Chikashige et al.], [Telomere-led premeiotic chromosome
    movement in fission yeast], [C7,
    C8], [bouquet\_ChikashigeTelomere1994],
    [1994], [Hewitt et al.], [Analysis of the tandem repeat locus D4Z4
    associated with facioscapulohumeral muscular dystrophy], [C5, C7,
    C8], [dux4\_d4z4\_fshd\_hewitt1994],
    [1996], [Health & Collaboration], [A complete set of human telomeric
    probes and their clinical application], [C5], [cytofound\_NIH1996],
    [1996], [Page & Shaffer], [Identification of a common breakpoint at
    14q21q from de novo and familial Robertsonian translocations:
    brea…], [C5, C8], [acrocentric\_rdna\_robertsonian\_page1996],
    [1996], [Scherthan et al.], [Centromere and telomere movements
    during early meiotic prophase of mouse and man are associated with
    the on…], [—], [Scherthan1996],
    [1996], [van Deutekom et al.], [FSHD-associated 4q35 D4Z4
    macrosatellite sequenced; chr4q↔chr10q D4Z4 translocation observed
    in \~20% of typed individuals.], [C5, C8], [vanDeutekom1996],
    [1997], [Flint et al.], [Two-domain subtelomere model: distal
    multi-end, proximal few-end, separated by internal (TTAGGG)n; the
    field’s predictive structural model.], [C4, C5], [Flint1997],
    [1998], [Graves et al.], [The origin and evolution of the
    pseudoautosomal regions of human sex chromosomes], [C4,
    C7], [sexchrompars\_graves1998],
    [1998], [Lupski], [Genomic disorders concept; recurrent
    rearrangements with predictable SD-flanked breakpoints.], [C7,
    C8], [Lupski1998],
    [1998], [Rouquier et al.], [Distribution of olfactory receptor genes
    in the human genome], [—], [Rouquier1998],
    [1998], [Trask et al.], [OR cluster blocks duplicated
    polymorphically near chromosome ends; founding empirical
    block-polymorphism paper.], [C5], [Trask1998],
    [1998], [Zickler & Kleckner], [The leptotene-zygotene transition of
    meiosis], [—], [ZicklerKleckner1998],
    [1999], [Lahn & Page], [Four evolutionary strata on the human X
    chromosome], [C4, C7], [sexchrompars\_lahn1999],
    [1999], [Liao], [Concerted-evolution synthesis: unequal crossover
    and gene conversion as the two principal molecular engines.], [C7,
    C8], [Liao1999],
    [1999], [Stout et al.], [Interphase chr4q/chr10q clustering by FISH
    in human muscle nuclei; first 3D-proximity claim for D4Z4-bearing
    arms.], [C7], [Stout1999],
    [1999], [Zickler & Kleckner], [Meiotic chromosomes: integrating
    structure and function], [C7, C8], [bouquet\_ZicklerKleckner1999],
    [2000], [Ciccodicola et al.], [Differentially regulated and evolved
    genes in the fully sequenced Xq/Yq pseudoautosomal region], [C4,
    C7], [sexchrompars\_ciccodicola2000],
    [2001], [Bandyopadhyay et al.], [Identification and characterization
    of satellite III subfamilies to the acrocentric chromosomes], [C5,
    C8], [acrocentric\_rdna\_robertsonian\_bandyopadhyay2001],
    [2001], [Consortium], [Initial sequencing and analysis of the human
    genome], [—], [IHGSC2001],
    [2001], [de Vries et al.], [Clinical studies on submicroscopic
    subtelomeric rearrangements: a
    checklist], [C5], [cytofound\_deVries2001],
    [2001], [Eichler], [Recent segmental duplications and human genome
    evolution; founding SD framework.], [C5], [Eichler2001],
    [2001], [Eichler], [Segmental duplications: an expanding role in
    genomic instability and disease], [C7,
    C8], [concerted\_evolution\_nahr\_Eichler2001],
    [2001], [Galtier et al.], [GC-biased gene conversion in humans;
    biased-substitution effect on isochore structure.], [C7,
    C8], [Galtier2001],
    [2001], [Glusman et al.], [Complete human olfactory subgenome
    catalogue; OR family is the largest mammalian gene family and a
    chromosome-end staple.], [C5], [Glusman2001],
    [2001], [Mefford et al.], [Comparative sequencing of a multicopy
    subtelomeric region containing olfactory receptor genes reveals
    multi…], [—], [Mefford2001],
    [2001], [Riethman et al.], [Integration of telomere clones with the
    draft human genome; founding paper of the systematic subtelomeric
    assembly programme.], [C2, C4], [Riethman2001],
    [2001], [Scherthan], [A bouquet makes ends meet], [C7,
    C8], [bouquet\_Scherthan2001],
    [2001], [Zozulya et al.], [The human olfactory receptor
    repertoire], [—], [Zozulya2001],
    [2002], [Bailey et al.], [Recent segmental duplications framework;
    SDs as the genome-wide arena of NAHR-mediated
    exchange.], [C5], [Bailey2002],
    [2002], [Bandyopadhyay et al.], [Parental origin and timing of de
    novo Robertsonian translocation formation], [C5,
    C8], [acrocentric\_rdna\_robertsonian\_bandyopadhyay2002],
    [2002], [Der-Sarkissian et al.], [Subtelomeric variation of human
    proterminal markers; population-level evidence for inter-arm block
    exchange.], [C6], [DerSarkissian2002],
    [2002], [Fan et al.], [Genomic structure and evolution of the
    ancestral chromosome fusion site in 2q13-2q14.1 and paralogous
    regio…], [—], [Fan2002],
    [2002], [Lee et al.], [Multiple sequence alignment using partial
    order graphs], [C2, C3], [pangenome\_graphs\_impg\_Lee2002],
    [2002], [Lemmers et al.], [FSHD associated with one of two 4q
    subtelomere variants (4qA / 4qB); first allele structure
    model.], [C5, C8], [Lemmers2002],
    [2002], [Martin et al.], [The evolutionary origin of human
    subtelomeric homologies–or where the ends begin], [—], [Martin2002],
    [2002], [Mefford & Trask], [Subtelomere transition zones:
    chromosome-specific DNA → terminal telomere with recurrent exchange
    and CNV; 'patchwork' model.], [C5, C6], [MeffordTrask2002],
    [2002], [Rosenberg et al.], [Genetic Structure of Human
    Populations], [C6], [subtel\_popgen\_rosenberg2002],
    [2002], [Samonte & Eichler], [Segmental duplications and the
    evolution of the primate genome], [C7,
    C8], [concerted\_evolution\_nahr\_SamonteEichler2002],
    [2002], [Stankiewicz & Lupski], [NAHR / Non-Allelic Homologous
    Recombination canonical paper; flanking SD architecture as substrate
    of ectopic exchange.], [C7, C8], [StankiewiczLupski2002],
    [2002], [van Geel et al.], [Genomic analysis of human chromosome 10q
    and 4q telomeres suggests a common origin], [C5, C7,
    C8], [dux4\_d4z4\_fshd\_vanGeel2002],
    [2002], [Young et al.], [Different evolutionary processes shaped the
    mouse and human olfactory receptor gene
    families], [—], [YoungFriedman2002],
    [2003], [Charchar et al.], [PAR2 origins via multiple
    species-specific events; PAR1 vs PAR2 evolutionary
    contrast.], [C4], [Charchar2003],
    [2003], [Flint & Knight], [The use of telomere probes to investigate
    submicroscopic rearrangements associated with mental
    retardation], [C5], [cytofound\_Flint2003],
    [2003], [Gilad et al.], [Natural selection on the olfactory receptor
    gene family in humans and
    chimpanzees], [—], [GiladNatSelection2003],
    [2003], [Niimura & Nei], [Olfactory-receptor evolutionary dynamics;
    comparative genomics of OR pseudogenisation.], [C5], [Niimura2003],
    [2003], [Rozen et al.], [Abundant gene conversion between arms of
    palindromes in human and ape Y chromosomes], [C7,
    C8], [concerted\_evolution\_nahr\_Rozen2003],
    [2003], [Scherthan], [Telomere dynamics unique to meiotic prophase:
    formation and significance of the bouquet], [C7,
    C8], [bouquet\_Scherthan2003],
    [2003], [Skaletsky et al.], [Reference Y-chromosome assembly placing
    PAR1, PAR2, X-degenerate, X-transposed and ampliconic territory in
    one mosaic.], [C4], [Skaletsky2003],
    [2004], [Filatov], [Gradient of silent substitution rate across
    human PAR consistent with sex-averaged recombination
    gradient.], [C4], [Filatov2004],
    [2004], [Gilad et al.], [Loss of olfactory receptor genes coincides
    with the acquisition of full trichromatic vision in
    primates], [—], [Gilad2004],
    [2004], [Harper et al.], [A bouquet of chromosomes], [C7,
    C8], [bouquet\_HarperBouquet2004],
    [2004], [Hurles], [Gene conversion between paralogous repeats: CMT1A
    as canonical exemplar of disease-relevant ectopic
    recombination.], [C7, C8], [Hurles2004],
    [2004], [Jeffreys & May], [Human gene-conversion tract length
    \~50–500 bp from sperm crossover hotspots.], [C7,
    C8], [JeffreysMay2004],
    [2004], [Lovett], [Mismatch / template-switching mechanisms
    producing slipped-strand and sub-block variation at the nucleotide
    level.], [C7], [Lovett2004],
    [2004], [Masny et al.], [chr4q35.2 nuclear-periphery localisation
    via lamin A/C; founding observation for D4Z4 envelope
    tethering.], [C5, C7], [Masny2004],
    [2004], [Nergadze et al.], [Insertions of telomeric repeats at
    intrachromosomal break sites during primate evolution], [C4,
    C5], [subtelstruct\_NergadzeITS2007],
    [2004], [Niimura & Nei], [The human olfactory receptor gene
    family], [—], [HumanORFamily2004],
    [2004], [Riethman et al.], [Subtelomeric assembly programme: \~80%
    of distal 100 kb of human chromosomes consists of shared duplicated
    blocks.], [C2, C4], [Riethman2004],
    [2005], [Linardopoulou et al.], [Subtelomere paralogous-block
    sequence catalogue; chr10p/chr18p hot-spot pair; founding empirical
    paper for inter-arm exchange.], [C5, C7], [Linardopoulou2005],
    [2005], [Lupski & Stankiewicz], [Genomic disorders synthesis: NAHR +
    flanking SD architecture in human disease.], [C7,
    C8], [LupskiStankiewicz2005],
    [2005], [Nei & Rooney], [Birth-and-death alternative to concerted
    evolution for some multigene families; explicit caveat for
    OR4F.], [C8], [NeiRooney2005],
    [2005], [Riethman], [Human subtelomere structure and
    variation], [C6], [subtel\_popgen\_riethman2005],
    [2005], [Ross et al.], [Reference X-chromosome assembly: PAR1 / PAR2
    \/ XAR / XCR / XTR named in one framework.], [C4], [Ross2005],
    [2006], [Bailey & Eichler], [Primate SDs as crucibles of evolution
    and disease; subtelomeric / pericentromeric SD
    bias.], [C7], [BaileyEichler2006],
    [2006], [Chikashige et al.], [Fission-yeast Bqt1/Bqt2 telomere
    tethers; SUN-KASH meiotic chromosome attachment to nuclear
    envelope.], [C7], [Chikashige2006],
    [2006], [Patterson et al.], [Population Structure and
    Eigenanalysis], [C6], [subtel\_popgen\_patterson2006],
    [2006], [Sharp et al.], [Segmental duplication update; SD landscape
    across the human genome with subtelomeric/pericentromeric
    bias.], [C5], [Sharp2006],
    [2006], [Trelles-Sticken et al.], [The meiotic bouquet promotes
    homolog interactions and restricts ectopic recombination in
    textitSchizosaccha…], [C7, C8], [bouquet\_TrellesBouquetPombe2005],
    [2007], [Ambrosini et al.], [Subtelomeric duplicon-block
    organisation: 11 subtelomere-specific block families, bimodal \~91%
    \/ \~98% identity peaks.], [C4, C5], [Ambrosini2007],
    [2007], [Chen et al.], [Gene-conversion review of record: SDSA, BIR
    mechanisms, disease examples, tract-length distributions.], [C7,
    C8], [ChenCooper2007],
    [2007], [Ding et al.], [SUN1 is required for telomere attachment to
    nuclear envelope and gametogenesis in mice], [C7,
    C8], [bouquet\_DingSUN12007],
    [2007], [Eickbush & Eickbush], [rDNA total repeat homogenisation in
    real time; recombination at the ribosomal
    locus.], [C8], [EickbushEickbush2007],
    [2007], [Ganley & Kobayashi], [Whole-genome quantification of
    intra-array homogenisation rates.], [C8], [GanleyKobayashi2007],
    [2007], [Keller et al.], [Genetic variation in a human odorant
    receptor alters odour perception], [—], [Keller2007],
    [2007], [Kowaljow et al.], [DUX4 within D4Z4 encodes a pro-apoptotic
    protein; molecular underpinning of FSHD.], [C5], [Kowaljow2007],
    [2007], [Lemmers et al.], [Specific sequence variations within the
    4q35 region are associated with facioscapulohumeral muscular
    dystrophy], [C5, C7, C8], [dux4\_d4z4\_fshd\_lemmers2007],
    [2007], [Mangs & Morris], [Pseudoautosomal region review (current
    genomics).], [C4], [MangsMorris2007],
    [2007], [Nergadze et al.], [Contribution of telomerase RNA
    retrotranscription to DNA double-strand break repair during
    mammalian genome…], [C4, C5], [subtelstruct\_Nergadze2007],
    [2007], [Nergadze et al.], [Endings in the middle: current knowledge
    of interstitial telomeric sequences], [C4,
    C5], [subtelstruct\_NergadzeITSReview2007],
    [2007], [Niimura & Nei], [Extensive gains and losses of olfactory
    receptor genes in mammalian evolution], [—], [Niimura2007],
    [2007], [Rudd et al.], [Sister-chromatid exchange rates elevated at
    chromosome ends; cell-biological reason chromosome ends accumulate
    exchange events.], [C5, C7], [Rudd2007],
    [2008], [Anderson et al.], [Molecular Population Genetics of
    Drosophila Subtelomeric DNA], [C6], [subtel\_popgen\_anderson2008],
    [2008], [Flaquer et al.], [The human pseudoautosomal regions: a
    review for genetic epidemiologists], [C4,
    C7], [sexchrompars\_flaquer2008],
    [2008], [Gu et al.], [Frequency of nonallelic homologous
    recombination is correlated with length of homology: evidence that
    ectop…], [C7, C8], [concerted\_evolution\_nahr\_Gu2011],
    [2008], [Hasin et al.], [High-resolution copy-number variation map
    reflects human olfactory receptor diversity and
    evolution], [—], [Hasin2008],
    [2008], [Koszul et al.], [Meiotic chromosomes move by linkage to
    dynamic actin cables with transduction of force through the
    nuclear…], [C7], [hic3d\_koszul2008],
    [2008], [Li et al.], [Worldwide Human Relationships Inferred from
    Genome-Wide Patterns of Variation], [C6], [subtel\_popgen\_li2008],
    [2008], [Llorente et al.], [Break-induced replication (BIR) as third
    mechanism distinct from synthesis-dependent strand
    annealing.], [C7], [Llorente2008],
    [2008], [Olender et al.], [Olfactory-receptor structured update; OR
    cluster annotation across human chromosomes.], [C5], [Olender2008],
    [2008], [Riethman], [Copy-number variation at subtelomeres in the
    broader CNV literature.], [C4], [Riethman2008],
    [2008], [Ruiz-Herrera et al.], [Internal telomeric sequences (ITS)
    review: mechanisms by which telomeric tracts end up far from
    chromosome ends.], [C4, C5], [RuizHerrera2008],
    [2008], [Stults et al.], [Genomic architecture and inheritance of
    human ribosomal RNA gene clusters], [C5,
    C8], [acrocentric\_Stults2008],
    [2008], [Young et al.], [Extensive copy-number variation of the
    human olfactory receptor gene family], [—], [Young2008CNV],
    [2009], [Benovoy & Drouin], [Ectopic gene conversions in the human
    genome], [C7, C8], [concerted\_evolution\_nahr\_Benovoy2009],
    [2009], [Duret & Galtier], [GC-biased gene conversion and isochore
    structure; population-genetic consequences.], [C7,
    C8], [DuretGaltier2009],
    [2009], [Flaquer et al.], [A new sex-specific genetic map of the
    human pseudoautosomal regions (PAR1 and PAR2)], [C4,
    C7], [sexchrompars\_flaquer2009],
    [2009], [Hasin-Brumshtein et al.], [OR genomic variation linked to
    phenotypic diversity in humans.], [C5], [Hasin2009],
    [2009], [Hastings et al.], [A microhomology-mediated break-induced
    replication model for the origin of human copy number
    variation], [C7, C8], [concerted\_evolution\_nahr\_Hastings2009],
    [2009], [Hastings et al.], [Mechanisms of change in gene copy
    number: NAHR, NHEJ, FoSTeS/MMBIR synthesis.], [C7,
    C8], [Hastings2009],
    [2009], [Hiraoka & Dernburg], [SUN-protein meiotic chromosome
    dynamics review; LINC complex meiotic role.], [C7], [Hiraoka2009],
    [2009], [Lieberman-Aiden et al.], [Hi-C introduced; bulk principles
    of human genome folding (compartments, contact frequency vs
    distance).], [C7], [LiebermanAiden2009],
    [2009], [Lupski], [Genomic disorders ten years on; NAHR + clinical
    genomics synthesis.], [C7, C8], [Lupski2009],
    [2009], [Niimura], [On the origin and evolution of vertebrate
    olfactory receptor genes: comparative genome analysis among 23
    ch…], [—], [Niimura2009],
    [2009], [Ottaviani et al.], [Perinuclear positioning element in
    subtelomeres requires A-type lamins and CTCF.], [C5,
    C7], [OttavianiGilson2008],
    [2009], [Ottaviani et al.], [D4Z4 acts as a CTCF / A-type
    lamin-dependent insulator; mechanistic envelope-tether
    anchor.], [C5, C7], [Ottaviani2009],
    [2009], [Penkner et al.], [SUN-1 modification cycle drives meiotic
    chromosome homology search through nuclear
    envelope.], [C7], [Penkner2009],
    [2009], [Rudd et al.], [Comparative sequence analysis of primate
    subtelomeres originating from a chromosome fission
    event], [—], [Rudd2009],
    [2009], [Snider et al.], [RNA transcripts, miRNA-sized fragments and
    proteins produced from D4Z4 units: new candidates for the
    pathop…], [C5, C7, C8], [dux4\_d4z4\_fshd\_snider2009],
    [2010], [Consortium], [A map of human genome variation from
    population-scale sequencing], [C6], [subtel\_popgen\_1000g2010],
    [2010], [de Greef et al.], [Common epigenetic changes of D4Z4 in
    contraction-dependent and contraction-independent FSHD], [C5, C7,
    C8], [dux4\_d4z4\_fshd\_degreef2010],
    [2010], [Kong et al.], [Fine-scale recombination map from deCODE
    pedigrees; sequence-level genetic maps of Iceland
    trios.], [C8], [Kong2010],
    [2010], [Lemmers et al.], [Unifying FSHD genetic model: permissive
    4qA haplotype produces stable DUX4 mRNA via 4q-specific
    polyadenylation signal.], [C5, C8], [Lemmers2010],
    [2010], [Lemmers et al.], [Worldwide population analysis of the 4q
    and 10q subtelomeres identifies only four discrete
    interchromosomal…], [C5, C7,
    C8], [dux4\_d4z4\_fshd\_lemmers2010worldwide],
    [2010], [Myers et al.], [Drive against hotspot motifs in primates
    implicates the PRDM9 gene in meiotic recombination], [C7,
    C8], [concerted\_evolution\_nahr\_Myers2010],
    [2010], [Snider et al.], [FSHD reflects incomplete suppression of a
    retrotransposed gene.], [C5], [Snider2010],
    [2010], [Stankiewicz & Lupski], [NAHR vocabulary, refined; recurrent
    rearrangements with predictable breakpoints at large SDs.], [C7,
    C8], [StankiewiczLupski2010],
    [2012], [Cabianca et al.], [D4Z4 copy-number controls a long
    ncRNA-mediated polycomb/trithorax epigenetic
    switch.], [C5], [Cabianca2012],
    [2012], [Consortium et al.], [An integrated map of genetic variation
    from 1,092 human genomes], [C3, C6], [hprc\_1000g2012],
    [2012], [Dixon et al.], [Topological domains in mammalian genomes
    identified by analysis of chromatin
    interactions], [C7], [hic3d\_dixon2012],
    [2012], [Geng et al.], [DUX4 activates germline genes,
    retroelements, and immune mediators: implications for
    facioscapulohumeral dy…], [C5, C7,
    C8], [dux4\_d4z4\_fshd\_geng2012],
    [2012], [Imakaev et al.], [Iterative correction of Hi-C data reveals
    hallmarks of chromosome organization], [C7], [hic3d\_imakaev2012],
    [2012], [Lemmers et al.], [Digenic inheritance of an SMCHD1 mutation
    and an FSHD-permissive D4Z4 allele causes facioscapulohumeral
    mus…], [C5, C7, C8], [dux4\_d4z4\_fshd\_lemmers2012],
    [2012], [Morimoto et al.], [A conserved KASH domain protein
    associates with telomeres, SUN1, and dynactin during mammalian
    meiosis], [C7, C8], [bouquet\_MorimotoKASH2012],
    [2012], [Pickrell & Pritchard], [Inference of Population Splits and
    Mixtures from Genome-Wide Allele Frequency
    Data], [C6], [subtel\_popgen\_pickrell2012],
    [2012], [Sosa et al.], [Structural basis of LINC complex formation
    through KASH-peptide / SUN-trimer binding.], [C7], [Sosa2012],
    [2013], [Anand et al.], [Break-induced replication characterised in
    yeast; long-tract synthesis from single homologous
    template.], [C7], [Anand2013],
    [2013], [Bhatia et al.], [Estimating and interpreting FST: The
    impact of rare variants], [C6], [subtel\_popgen\_bhatia2013],
    [2013], [Boateng et al.], [Homologous pairing precedes
    SPO11-mediated double-strand breaks in mice.], [C7], [Boateng2013],
    [2013], [Chen et al.], [Signals of historical interlocus gene
    conversion in human segmental duplications], [C7,
    C8], [concerted\_evolution\_nahr\_Chen2013],
    [2013], [Floutsakou et al.], [Acrocentric distal/proximal junctions
    sequenced; PJ and DJ are 95–99% identical across the five
    acrocentrics.], [C5, C8], [Floutsakou2013],
    [2013], [Horn et al.], [A mammalian KASH domain protein coupling
    meiotic chromosomes to the cytoskeleton], [C7,
    C8], [bouquet\_HornKASH52013],
    [2013], [Nagano et al.], [Single-cell Hi-C; chromosome-conformation
    variability across individual cells.], [C7], [Nagano2013],
    [2014], [Hinch et al.], [Recombination in the human pseudoautosomal
    region PAR1], [C4, C7], [sexchrompars\_hinch2014],
    [2014], [Jarmuz-Szymczak et al.], [Narrowing the localization of the
    region breakpoint in most frequent Robertsonian
    translocations], [C5,
    C8], [acrocentric\_rdna\_robertsonian\_jarmuzSzymczak2014],
    [2014], [Mainland et al.], [The missense of smell: functional
    variability in the human odorant receptor
    repertoire], [—], [Mainland2014],
    [2014], [Mensah et al.], [Pseudoautosomal region 1 length
    polymorphism in the human population], [C4,
    C7], [sexchrompars\_mensah2014],
    [2014], [Niimura et al.], [Human OR repertoire in 13-mammal
    comparative frame; canonical OR-pseudogenisation gradient
    paper.], [C5], [Niimura2014],
    [2014], [Rao et al.], [Hi-C at kilobase resolution; chromatin loops,
    contact domains, CTCF anchoring.], [C7], [Rao2014],
    [2014], [Stimpson et al.], [Nucleolar organization, ribosomal DNA
    array stability, and acrocentric chromosome integrity are linked to
    t…], [C5, C8], [acrocentric\_rdna\_robertsonian\_stimpson2014],
    [2014], [Stong et al.], [Subtelomeric assembly extended into
    structural and chromatin annotation; CTCF/cohesin organisation
    across haplotype-resolved subtelomeres.], [C2, C4], [Stong2014],
    [2015], [Consortium], [1000 Genomes Project main paper;
    population-genetic substrate for subtelomere allele-frequency
    analyses.], [C6], [ThousandGenomes2015],
    [2015], [Lee et al.], [Rapid mouse telomere prophase movements
    measured by live imaging.], [C7], [Lee2015],
    [2015], [Shibuya et al.], [Mechanism and regulation of rapid
    telomere prophase movements in mouse meiotic chromosomes], [C7,
    C8], [bouquet\_ShibuyaRPMs2015],
    [2015], [Sudmant et al.], [1000 Genomes structural-variation map;
    population-scale CNV catalogue.], [C5, C6], [Sudmant2015],
    [2015], [Williams et al.], [Direct pedigree-resolved gene conversion
    calls at sequence resolution.], [C7, C8], [Williams2015],
    [2015], [Zickler & Kleckner], [Recombination, pairing and synapsis
    review; LINC complex / SUN-KASH machinery
    framing.], [C7], [ZicklerKleckner2015],
    [2016], [Cotter et al.], [Genetic diversity on the human X
    chromosome does not support a strict pseudoautosomal boundary], [C4,
    C7], [sexchrompars\_cotter2016],
    [2016], [Mallick et al.], [Simons Genome Diversity Project; 300
    deeply sequenced genomes from 142 worldwide
    populations.], [C6], [Mallick2016],
    [2016], [McStay], [Nucleolar organiser regions as 'genomic dark
    matter'; chromosomal context as determinant of NOR stability.], [C5,
    C8], [McStay2016],
    [2017], [Bh\'erer et al.], [Sex-specific recombination map;
    large-scale pedigree analysis.], [C8], [Bherer2017],
    [2017], [Chu et al.], [PAR-TERRA directs homologous sex chromosome
    pairing], [C4, C7], [sexchrompars\_chu2017],
    [2017], [Hendrickson et al.], [Conserved roles of mouse DUX and
    human DUX4 in activating cleavage-stage genes and MERVL/HERVL
    retrotranspo…], [C5, C7, C8], [dux4\_d4z4\_fshd\_hendrickson2017],
    [2017], [Ramani et al.], [Massively multiplex single-cell Hi-C;
    combinatorial barcoding for population-scale single-cell
    3D.], [C7], [Ramani2017],
    [2017], [Stevens et al.], [3D structures of individual mammalian
    genomes resolved at single-cell level.], [C7], [Stevens2017],
    [2018], [Consortium], [Computational pan-genomics consortium review;
    consolidated early conceptual frame.], [C2], [Computational2018],
    [2018], [Garrison et al.], [vg variation-graph toolkit; genotyping
    and read-mapping on graphs.], [C2], [Garrison2018],
    [2018], [Jain et al.], [A fast adaptive algorithm for computing
    whole-genome homology maps], [C2,
    C3], [pangenome\_graphs\_impg\_Jain2018],
    [2018], [Poriswanish et al.], [Recombination hotspots in an extended
    human pseudoautosomal domain predicted from double-strand break
    maps…], [C4, C7], [sexchrompars\_poriswanish2018],
    [2018], [Shao et al.], [BioNano + nanopore optical mapping detects
    ongoing chromosome-end extension events.], [C2, C4], [Shao2018],
    [2018], [Tan et al.], [Dip-C; per-haplotype 3D structures of GM12878
    and PBMC single cells.], [C7], [Tan2018],
    [2018], [Wolff et al.], [Galaxy HiCExplorer: a web server for
    reproducible Hi-C data analysis, quality control and
    visualization], [C7], [hic3d\_wolff2018],
    [2019], [Alavattam et al.], [Attenuated chromatin
    compartmentalization in meiosis and its maturation in sperm
    development], [C7], [hic3d\_alavattam2019],
    [2019], [Audano et al.], [Long-read structural-variation discovery
    exposes many variants invisible to short
    reads.], [C8], [Audano2019],
    [2019], [Blokhina et al.], [The telomere bouquet is a hub where
    meiotic double-strand breaks, synapsis, and stable homolog
    juxtapositio…], [C7, C8], [bouquet\_BlokhinaZebrafish2019],
    [2019], [Chaisson et al.], [Multi-platform discovery of structural
    variation in 15 human genomes; long-read SV
    catalogue.], [C8], [Chaisson2019],
    [2019], [Halld\'orsson et al.], [Sequence-level genetic map from
    deCODE pedigrees; direct gene-conversion calls.], [C7,
    C8], [Halldorsson2019],
    [2019], [Levy-Sakin et al.], [Genome maps across 26 human
    populations reveal population-specific patterns of structural
    variation], [C6], [subtel\_popgen\_levysakin2019],
    [2019], [Patel et al.], [Mouse meiotic Hi-C; chromosome-end
    clustering during meiotic prophase.], [C7], [Patel2019],
    [2019], [Qiang et al.], [The meiotic TERB1–TERB2–MAJIN complex
    tethers telomeres to the nuclear envelope], [C7,
    C8], [bouquet\_QiangTERB2019],
    [2019], [Sasani et al.], [Three-generation Utah families;
    grandchildren resolve postzygotic vs germline
    mosaicism.], [C8], [Sasani2019],
    [2019], [Ulahannan et al.], [Pore-C long-read concatemer Hi-C
    variant; multi-way contact resolution.], [C7], [Ulahannan2019],
    [2019], [van Sluis et al.], [Acrocentric short-arm architecture
    extended; near-identity of acrocentric short arms in human and
    chimpanzee.], [C5, C8], [vanSluis2019],
    [2019], [Vara et al.], [Three-Dimensional Genomic Structure and
    Cohesin Occupancy Correlate with Transcriptional Activity during
    Sp…], [C7], [hic3d\_vara2019],
    [2019], [Wang et al.], [Structural basis of meiotic telomere
    attachment to the nuclear envelope by MAJIN–TERB2–TERB1], [C7,
    C8], [bouquet\_WangTERB2019],
    [2020], [Acquaviva et al.], [Ensuring meiotic DNA break formation in
    the mouse pseudoautosomal region], [C4,
    C7], [sexchrompars\_acquaviva2020],
    [2020], [Bell et al.], [Insights into variation in meiosis from
    31,228 human sperm genomes], [—], [Bell2020],
    [2020], [Bergstr\"om et al.], [Human Genome Diversity Project; 929
    genomes from 54 populations.], [C6], [Bergstrom2020],
    [2020], [Collins et al.], [Structural-variation reference at
    population scale; gnomAD-SV.], [C8], [Collins2020],
    [2020], [Eizenga et al.], [Pangenome graphs canonical
    review.], [C2], [Eizenga2020],
    [2020], [Fukami et al.], [Human spermatogenesis tolerates massive
    size reduction of the pseudoautosomal region], [C4,
    C7], [sexchrompars\_fukami2020],
    [2020], [Kota & Bhatt], [Tethering of Telomeres to the Nuclear
    Envelope Is Mediated by SUN1–MAJIN and Possibly Promoted by
    SPDYA–C…], [C7, C8], [bouquet\_KotaSUN1MAJIN2020],
    [2020], [Li et al.], [minigraph; reference-based pangenome graph
    construction.], [C2], [Li2020minigraph],
    [2020], [Miga et al.], [Telomere-to-telomere assembly of a complete
    human X chromosome], [C4, C7], [sexchrompars\_miga2020],
    [2020], [Stergachis et al.], [Fiber-seq for single-molecule
    chromatin-fibre regulatory architectures.], [C7], [Stergachis2020],
    [2020], [West & others], [The TERB1–TERB2–MAJIN complex of mouse
    meiotic telomeres dates back to the common ancestor of
    metazoans], [C7, C8], [bouquet\_BhattTERBEvolution2020],
    [2020], [Young et al.], [Whole-genome optical mapping of 154
    individuals; widespread large-scale subtelomeric SV.], [C2,
    C4], [Young2020],
    [2020], [Zook et al.], [Genome-in-a-Bottle benchmarking framework
    for difficult variants and SVs.], [C8], [Zook2020],
    [2021], [Bhatt & others], [Control of meiotic chromosomal bouquet
    and germ cell morphogenesis by the zygotene cilium], [C7,
    C8], [bouquet\_ZygoteneCilium2021],
    [2021], [Cheng et al.], [hifiasm haplotype-resolved assembly; one of
    the two HPRC-grade assemblers.], [C8], [Cheng2021],
    [2021], [Ebert et al.], [Haplotype-resolved diverse human genomes;
    long-read SV discovery at population scale.], [C8], [Ebert2021],
    [2021], [Grigorev et al.], [Long-read telomere sequencing; telomeric
    repeat-motif heterogeneity across individuals.], [C2,
    C4], [Grigorev2021],
    [2021], [Logsdon et al.], [Complete chr8 telomere-to-telomere
    assembly; closed many formerly inaccessible
    regions.], [C2], [Logsdon2021],
    [2021], [Monteiro et al.], [Evolutionary dynamics of the human
    pseudoautosomal regions], [C4, C7], [sexchrompars\_monteiro2021],
    [2021], [Sirén et al.], [Pangenomics enables genotyping of known
    structural variants in 5202 diverse genomes], [C2,
    C3], [pangenome\_graphs\_impg\_Siren2021],
    [2021], [Xu et al.], [The SUN1–SPDYA interaction plays an essential
    role in meiosis prophase I], [C7, C8], [bouquet\_XuSPDYA2021],
    [2021], [Zuo et al.], [Stage-resolved mouse meiotic Hi-C; SUN1 W151R
    compresses chromosome-end alignment from \~20% to
    \~5%.], [C7], [Zuo2021],
    [2022], [Adam et al.], [NPGREAT pipeline integrates linked reads
    with ultralong nanopore reads to assemble CHM13 subtelomeres.], [C2,
    C4], [Adam2022],
    [2022], [Altemose et al.], [Complete genomic and epigenetic maps of
    human centromeres], [C5, C8], [acrocentric\_Altemose2022],
    [2022], [Bergman & Schierup], [Evolutionary dynamics of
    pseudoautosomal region 1 in humans and great apes], [C4,
    C7], [sexchrompars\_bergman2022],
    [2022], [Deshpande et al.], [Identifying synergistic high-order 3D
    chromatin conformations from genome-scale nanopore concatemer
    sequencing], [C7], [hic3d\_deshpande2022],
    [2022], [Ebler et al.], [Pangenome-based genome inference allows
    efficient and accurate genotyping across a wide spectrum of
    variant…], [C3, C6], [hprc\_ebler2022],
    [2022], [Gershman et al.], [ENCODE CTCF ChIP-seq re-aligned to
    T2T-CHM13; CTCF enrichment at TAR loci.], [C7], [Gershman2022],
    [2022], [Guarracino et al.], [ODGI: understanding pangenome
    graphs], [C2, C3], [pangenome\_graphs\_impg\_GuarracinoHeumos2022],
    [2022], [Logsdon et al.], [Short arms of human acrocentric
    chromosomes and the completion of the human genome sequence], [C5,
    C8], [acrocentric\_ShortArms2022],
    [2022], [Nurk et al.], [T2T-CHM13 telomere-to-telomere reference
    assembly.], [C2, C3], [Nurk2022],
    [2022], [Sholes et al.], [Single-molecule architecture and
    heterogeneity of human telomeric DNA and chromatin], [C4,
    C5], [subtelstruct\_Sholes2022],
    [2022], [Sirén et al.], [GBZ file format for pangenome graphs], [C2,
    C3], [pangenome\_graphs\_impg\_Siren2022],
    [2022], [Wagner et al.], [Benchmarking framework for difficult
    variants and SVs (GIAB).], [C8], [Wagner2022],
    [2022], [Wang et al.], [The Human Pangenome Project: a global
    resource to map genomic diversity], [C3, C6], [hprc\_wang2022],
    [2023], [Ahsan et al.], [Long-read SV detection algorithm tradeoff
    review.], [C8], [Ahsan2023],
    [2023], [Garrison & Guarracino], [Unbiased pangenome graphs], [C2,
    C3], [pangenome\_graphs\_impg\_GarrisonGuarracino2023],
    [2023], [Guarracino et al.], [IMPG: implicit pangenome graph], [C2,
    C3], [pangenome\_graphs\_impg\_IMPG2023],
    [2023], [Guarracino et al.], [Pangenome-graph evidence of
    recombination between heterologous human acrocentric chromosomes;
    closest predecessor.], [C5, C8], [Guarracino2023],
    [2023], [Hallast et al.], [Assembly of 43 human Y chromosomes
    reveals extensive complexity and variation], [C4,
    C7], [sexchrompars\_hallast2023],
    [2023], [Kucuk et al.], [Comprehensive de novo mutation discovery
    with HiFi long-read sequencing], [—], [Kucuk2023],
    [2023], [Liao et al.], [HPRC v1 draft pangenome; 47 individuals × 2
    haplotypes.], [C2, C3], [Liao2023],
    [2023], [Rautiainen et al.], [Verkko diploid T2T assembly
    pipeline.], [C8], [Rautiainen2023],
    [2023], [Rhie et al.], [The complete sequence of a human Y
    chromosome], [C4, C7], [sexchrompars\_rhie2023],
    [2023], [Tan et al.], [scNanoHi-C: a single-cell long-read
    concatemer sequencing method to reveal high-order chromatin
    structures…], [C7], [hic3d\_scnanoHiC2023],
    [2023], [Vollger et al.], [Increased mutation and gene conversion
    within human segmental duplications], [C7,
    C8], [concerted\_evolution\_nahr\_Vollger2023],
    [2023], [Vollger et al.], [Segmental-duplication landscape at T2T
    scale.], [C5], [Vollger2023],
    [2024], [Bellott et al.], [Where is the boundary of the human
    pseudoautosomal region?], [C4, C7], [sexchrompars\_bellott2024],
    [2024], [de Lima & Gerton], [A working model for the formation of
    Robertson chromosomes], [C5, C8], [acrocentric\_workingmodel2024],
    [2024], [Garrison et al.], [pggb (pangenome graph builder);
    reference-free all-vs-all graph
    construction.], [C2], [Garrison2024pggb],
    [2024], [Heumos et al.], [odgi infrastructure; pairwise Jaccard
    similarity computation on pggb GFA.], [C2], [Heumos2024],
    [2024], [Hickey et al.], [Pangenome graph construction from genome
    alignments with Minigraph-Cactus], [C2,
    C3], [pangenome\_graphs\_impg\_Hickey2024],
    [2024], [Liu et al.], [Long-read SV detection benchmarking;
    algorithm tradeoffs.], [C8], [Liu2024],
    [2024], [Logsdon et al.], [T2T-resolution subtelomere assemblies;
    sequence-resolution structural maps.], [C2], [Logsdon2024],
    [2024], [Ma et al.], [Exchange of subtelomeric regions between
    chromosomes 4q and 10q reverts the FSHD genotype and
    phenotype], [—], [Ma2024],
    [2024], [Smolka et al.], [Detection of mosaic and population-level
    structural variants with Sniffles2], [—], [Smolka2024],
    [2025], [Cechova et al.], [Origin and evolution of acrocentric
    chromosomes in human and great apes], [C5,
    C8], [acrocentric\_rdna\_robertsonian\_cechovaHartley2025apes],
    [2025], [Cechova & others], [WashU 3-generation T2T pedigree
    (PAN010/PAN011 → PAN027 → PAN028); primary substrate for chapter 14
    inter-chromosomal patches.], [C8], [Cechova2025],
    [2025], [Consortium], [HPRC Data Release 2: high-quality phased
    genome assemblies from 232 diverse individuals], [C3,
    C6], [hprc\_hprcv2\_2025],
    [2025], [Cramer & others], [CiFi: accurate long-read chromosome
    conformation capture with low-input
    requirements], [C7], [hic3d\_cifi2025],
    [2025], [de Lima et al.], [Pangenome-resolved Robertsonian
    breakpoint analysis; SST1 macrosatellite as recombining
    homology.], [C5, C8], [deLima2025],
    [2025], [de Lima et al.], [The formation and propagation of human
    Robertsonian chromosomes], [C5,
    C8], [acrocentric\_rdna\_robertsonian\_deLima2025],
    [2025], [Francis & others], [Long-read pedigree-based variant
    detection in HPRC-style assemblies.], [C8], [Francis2025],
    [2025], [Guarracino & others], [Origin and evolution of acrocentric
    chromosomes in human and great apes], [C5,
    C8], [acrocentric\_Guarracino2025ape],
    [2025], [Lalli et al.], [HPRC-related haplotype-resolution
    analysis.], [C8], [Lalli2025],
    [2025], [Poriswanish et al.], [Multiple origins and phenotypic
    implications of an extended human pseudoautosomal region shown by
    analysis…], [C4, C7], [sexchrompars\_poriswanish2025],
    [2025], [Porubsky & others], [Human acrocentric chromosome short arm
    de novo mutation and recombination], [C5,
    C8], [acrocentric\_Porubsky2025denovo],
    [2025], [Porubsky & others], [CEPH1463 4-generation pedigree;
    cross-assembler validation of inter-chromosomal subtelomeric
    features.], [C8], [Porubsky2025],
    [2025], [Sir\'en et al.], [Structural variation in 1,019 diverse
    humans based on long-read sequencing], [C3, C6], [hprc\_siren2025],
    [2025], [Tan & others], [Single-cell long-read Hi-C, scNanoHi-C2,
    details 3D genome reorganization in embryonic-stage germ
    cells], [C7], [hic3d\_scnanoHiC2\_2025],
    [2025], [Xu & others], [Single-cell 3D structure of 20 human sperm
    cells.], [C7], [Xu2025],
    [2026], [Cechova & others], [Biobank-scale genotyping of Robertson
    translocations reveals hidden structural variation on the human
    acroc…], [C5, C8], [acrocentric\_Cechova2026],
    [2026], [Hartley et al.], [Biobank-scale genotyping of Robertsonian
    translocations reveals hidden structural variation on the human
    ac…], [C5,
    C8], [acrocentric\_rdna\_robertsonian\_hartley2026biobank],
  )]
  , kind: table
  )

#emph[Total: 295 entries.]

#emph[End of CHRONOLOGY\_v2.md.]

#pagebreak()
= Appendix B — Open Questions and Gaps

#strong[Author:] lit-review-v2 (agent-937), 2026-05-06. Supersedes
`GAPS.md` (lit-review-synthesis, agent-877). #strong[Companion files:]
`SYNTHESIS_v2.md`, `CHRONOLOGY_v2.md`,
`paper_prep/synthesis/REFERENCES_v3.bib`,
`paper_prep/synthesis/CROSSWALK.md`.

This is a one-page note flagging anything the v2 synthesis pass surfaced
as missing, mis-cited, or open. It is not a follow-up plan; only a list
of items that future passes (or revisions to the manuscript) should
resolve. The v1 GAPS.md is retained as historical record; almost all of
its entries carry forward unchanged.

== Substrate provenance
<substrate-provenance>
All 14 topical reviews are committed under `paper_prep/lit_review/`:
`topic_01` cytogenetic foundations (21 bib entries), `topic_02` TAR1 /
TTAGGG / ITS structure (16), `topic_03` pseudohomologous regions concept
(14), `topic_04` sex-chromosome pseudoautosomal regions (28), `topic_05`
acrocentric / rDNA / Robertsonian (22 — v3-opus refresh), `topic_06`
D4Z4 / DUX4 / FSHD (13), `topic_07` concerted evolution / NAHR / gene
conversion (24 — v3-opus refresh), `topic_08` meiotic bouquet /
nuclear-envelope tethering (19), `topic_09` Hi-C / Pore-C / single-cell
\/ meiotic 3D (16), `topic_10` pangenome graphs / IMPG (15), `topic_11`
pedigree-based recombination detection (20), `topic_12` HPRC v1 / v2
population pangenomes (17), `topic_13` subtelomere population genetics /
Fst / out-of-Africa (12), `topic_14` OR / OR4F gradient (17). Total raw
substrate: 254 entries across the 14 topic bibs, with 24 in the original
`REFERENCES.bib`.

`REFERENCES_v3.bib` merges `REFERENCES.bib` (24), the prior synth-pass
merge `REFERENCES_v2.bib` (270; carried forward for canonical-key
continuity), and all 14 topic\_NN\_\*.bib (254 raw, 240 unique after
intra-topic dedup) into 295 deduplicated entries. Dedup is by DOI (when
present) and by (first-author surname, year, normalized-title prefix)
with the canonical / prefix-free key preferred on collision so existing
`SYNTHESIS_v2.md` citations resolve unchanged. All 168 distinct
citations in `SYNTHESIS_v2.md` resolve in `REFERENCES_v3.bib` (verified
2026-05-06).

The substantive change from v1 to v2 is the v3-opus refresh of two
topical reviews: - #strong[topic\_05] grew from \~16 entries to 22,
adding Schmickel 1973, Krystal 1981, Worton 1988, Page 1996,
Bandyopadhyay 2001/2002, Stults 2008, Stimpson 2014, Cechova/Hartley
2025 (great-ape acrocentric atlas) and Hartley 2026 (UK Biobank
ROB-carrier genotyping). - #strong[topic\_07] grew from \~18 entries to
24, adding Brown 1972, Lupski 1998, Lupski & Stankiewicz 2005, Bailey &
Eichler 2006, Eickbush & Eickbush 2007, Ganley & Kobayashi 2007,
Llorente 2008, Hastings 2009, Anand 2013, Williams 2015, Logsdon 2024
and the explicit Nei & Rooney 2005 birth-and-death alternative.

Both refreshes are integrated into `SYNTHESIS_v2.md` Parts II
(acrocentric/rDNA section) and III (concerted evolution / NAHR section).
The closing paragraph of `SYNTHESIS_v2.md` Part V calls out the specific
upgrades.

== Bibliography-level gaps (mostly carried forward from v1)
<bibliography-level-gaps-mostly-carried-forward-from-v1>
- #strong[`Cechova2025`] is cited as "(in press)" with no DOI; the WashU
  pedigree paper is not yet on bioRxiv per the lead author’s
  clarification. Update DOI once available.
- #strong[`Porubsky2025`] is cited as "(in press)" with no DOI; the
  CEPH1463 paper was published in #emph[Nature] (April 2025) per the
  lead author’s recent comment. Update with the published DOI before
  submission.
- #strong[`deLima2025`] is "(in submission)"; verify status before
  submission. (Note: topic\_05 v3-opus refresh assigns this paper a
  sequence-resolved SST1-mechanism finding that is now load-bearing for
  C8.)
- #strong[`Xu2025`] lists "Hanbo Xu and others" — replace with full
  author list once the manuscript is in its published form.
- #strong[`Francis2025`] has full DOI (`10.1038/s41588-025-02358-1`) but
  author list is currently truncated to "Brittany A. Francis and others"
  — replace with full author list before submission.
- #strong[Stergachis Fiber-seq telomere companion paper.]
  `REFERENCES_v3.bib` cites `Stergachis2020` with a note that the
  specific Fiber-seq-at-39-of-46-telomeres paper is the relevant
  downstream work — track that paper down before submission.
- #strong[`StankiewiczLupski2002` and `StankiewiczLupski2010`] are both
  cited; verify the 2002 #emph[Trends in Genetics] paper and the 2010
  #emph[Annual Review of Medicine] paper resolve to the intended NAHR
  vocabulary in context.
- #strong[`Computational2018` (consortium-author key)] is the canonical
  citation for the #emph[Computational pan-genomics] review (Marschall,
  Marz, Abeel et al., #emph[Briefings in Bioinformatics] 2018,
  doi:10.1093/bib/bbw089). Duplicate author-led keys are removed in
  `REFERENCES_v3.bib`.
- #strong[Topic-prefix vs canonical-key duplication.] Several topical
  reviews use prefixed keys (e.g.~`sexchrompars_charchar2003`,
  `acrocentric_rdna_robertsonian_henderson1972`,
  `dux4_d4z4_fshd_lemmers2010`) for papers that the existing
  `REFERENCES.bib` / prior synthesis pass cite as canonical PascalCase
  (`Charchar2003`, `Henderson1972`, `Lemmers2010`). The merge into
  `REFERENCES_v3.bib` keeps the canonical key on each duplicate (per the
  dedup rule: prefer prefix-free) so `SYNTHESIS_v2.md` and
  `CHRONOLOGY_v2.md` citations continue to resolve. The `topic_NN_*.md`
  files (read-only) still use their prefix keys and must therefore be
  compiled against `topic_NN_*.bib` directly, not against
  `REFERENCES_v3.bib` alone. The unique-to-topic prefix keys survive in
  `REFERENCES_v3.bib` unchanged.
- #strong[One renamed citation in v2.] `SYNTHESIS_v2.md` cites
  `@ChenCooper2007` rather than the v1 `@Chen2007`; the entry is the
  same paper (Chen, Cooper, Krawczak, #emph[Hum. Mutat.] 2007
  gene-conversion review). The v3 dedup preferred the longer
  `ChenCooper2007` key on the (Chen, 2007, gene-conversion) collision
  because both keys are non-prefixed.

== Substantive gaps
<substantive-gaps>
- #strong[C7 — direct LAD / Lamin B1 ChIP overlay.] The abstract’s
  "facilitated by the physical proximity of subtelomeres at the nuclear
  envelope" wording rests on Dip-C radial position as a proxy and on the
  Masny 2004 / Ottaviani 2009 mechanism for D4Z4 — chapter 05 §"Nuclear
  lamina cross-reference" does not run a genome-wide LAD intersection.
  `REWRITE_PLAN.md` TASK-22 is the correct follow-up. Fallback per
  `CROSSWALK §6`: soften to "consistent with peripheral telomere
  positioning."
- #strong[C5 — neighbour-joining tree] is named in the abstract but
  Andrea uses Leiden + UPGMA (12 of 15 communities agree exactly).
  One-line fix is `ape::nj()` on the existing arm-distance matrix;
  alternatively, relax the abstract wording to "cladistic analysis
  (Leiden / UPGMA)." Per `REWRITE_PLAN.md` TASK-01 / TASK-13.
- #strong[C2 — "\~12% pairwise sampling" computation.] Methods must
  compute the wfmash k-mer-evaluation rate from the on-disk PAFs (the
  realised value, not asserted) and write the Erdős-Rényi argument
  explicitly (n=18,827; p\* = log(n)/n ≈ 5.21e-4; 12% is \~230× above).
  Per `REWRITE_PLAN.md` TASK-10 / TASK-11. The argument is in
  `CROSSWALK §7b`.
- #strong[C2 / C3 dataset count.] Andrea reports 465 throughout; the
  abstract reports 466. Resolve canonically as 466 = 233 individuals × 2
  haplotypes + CHM13 reference (per `CROSSWALK §7c`).
- #strong[C8 — "concerted evolution" loose-sense framing.] The
  Discussion paragraph should explicitly state the loose-sense usage and
  tie it to the pedigree (chapter 14) as the empirical anchor. Per
  `CROSSWALK §3` and `REWRITE_PLAN.md` TASK-19. The v3-opus topic\_07
  refresh adds Nei & Rooney 2005 @NeiRooney2005 as the explicit
  alternative model: not all subtelomeric paralog families behave like
  rDNA, and the loose-sense framing is what allows the abstract to
  encompass both homogenising (rDNA) and turnover-prone (OR4F) systems.

== New gaps surfaced by the v2 substrate
<new-gaps-surfaced-by-the-v2-substrate>
- #strong[C5 — acrocentric per-meiosis rate.] Stults et al.~2008
  @acrocentric_Stults2008 reported \>10% rDNA-cluster rearrangement per
  meiosis using optical pulsed-field gel electrophoresis. Andrea’s
  chapter 14 reports 133 gene-conversion-like patches across the WashU
  pedigree but does not separate intra-NOR conversion from inter-arm
  transfer at sequence resolution. A direct comparison of Stults’ bulk
  rate vs the pedigree’s per-arm rate would either confirm or revise the
  \>10% number — useful for the Discussion.
- #strong[C5 — rare Robertsonian products (rob(13;15), rob(14;15),
  rob(15;22)) lack sequence-resolved breakpoints.] de Lima/Guarracino
  2025 @deLima2025 resolved rob(13q14q) and rob(14q21q) via the SST1
  mechanism. Whether the rarer products use the same SST1 mechanism or a
  different homology block (rDNA itself, satellite III subfamilies,
  distal junction) is open; the Hartley 2026
  @acrocentric_rdna_robertsonian_hartley2026biobank UK Biobank cohort
  may eventually carry enough rare carriers to address this.
- #strong[C8 — relative contribution of NAHR vs gene conversion vs BIR.]
  The 133 gene-conversion-like + 16 crossover-like + 229 acros-like
  patches in chapter 14 imply gene conversion dominates by an order of
  magnitude, but BIR-mediated tracts of 100+ kb may be miscategorised as
  gene conversion under the current sandwich-pattern heuristic. Anand et
  al.~2013 @Anand2013 flags exactly this miscategorisation risk in
  yeast. A future pass could add a tract-length-aware classifier that
  distinguishes classical gene conversion (50–500 bp, per Jeffreys & May
  2004 @JeffreysMay2004) from BIR-style long tracts (10–100+ kb, per
  Llorente 2008 @Llorente2008).
- #strong[C8 — gBGC fingerprint at PHR boundaries.] Galtier 2001
  @Galtier2001 / Duret & Galtier 2009 @DuretGaltier2009 establish that
  gene conversion is GC-biased. A per-PHR-boundary GC-content analysis
  on the existing pangenome graph is testable from current data and
  would either confirm or reject gBGC as a contributor to subtelomeric
  base composition.

== Open empirical questions surfaced by the synthesis (carried forward)
<open-empirical-questions-surfaced-by-the-synthesis-carried-forward>
- #strong[Crossover-rate ↔ cross-arm-affinity correlation honest
  negative.] End-to-end-report chapter 12 testable prediction \#7
  reports `rho = −0.43, p = 0.006` across all 39 arms but `rho = 0.00`
  after excluding 7 short-read-confounded arms. Open question: does a
  higher-quality recombination map (deeper sequencing, longer reads at
  acrocentric / PAR arms) recover the predicted negative correlation?
  Probably yes, but cannot be answered without new data.
- #strong[LINC-complex / SUN1-mutant test of meiotic alignment at PHR
  scale.] Zuo 2021 wild-type vs SUN1-W151R zygotene Hi-C is available
  (GEO: GSE155142, GSE155638, GSE155967). The PHR-median 105 kb scale is
  well below the \~5% mutant tip zone; tip contacts may be maintained
  without LINC-mediated force transmission, so the predicted effect is
  uncertain. Per chapter 12 testable prediction \#1.
- #strong[Somatic vs germline exchange.] HPRC v2 sample-source metadata
  (LCL vs blood) is not currently available; the manuscript cannot
  directly distinguish meiotic from somatic-LCL exchange in chapter 14
  patches. The pedigree’s #emph[inheritance] across three generations
  does rule out purely somatic origin for the transmitted patches but
  does not directly establish meiotic timing in cells.
- #strong[Subtelomeric Fiber-seq / nucleosome footprint at PHR
  boundaries.] The Stergachis 2020 framework supports single-molecule
  CTCF / nucleosome footprinting; a per-PHR-boundary analysis at 39 / 46
  telomeres would either confirm or reject the CTCF-cohesin boundary
  prediction (chapter 12 testable prediction \#4). Not currently in
  scope but a clean follow-up.

== What this synthesis does #emph[not] claim
<what-this-synthesis-does-not-claim>
This document does not extend any of the abstract’s claims beyond what
`paper_prep/synthesis/CROSSWALK.md`,
`paper_prep/synthesis/REFERENCES.bib`, the 14 committed topical reviews
(`topic_01_*` through `topic_14_*`), the v3-opus topic\_05 and topic\_07
refreshes, and `end-to-end-report/report/12_literature.md` already
establish. The prose in `SYNTHESIS_v2.md` operates at the
field-history-and-framework level; where deeper per-paper analysis is
needed, the corresponding `topic_NN_*.md` provides it. No new analytical
claim is made; this is a literature synthesis only.

#emph[End of GAPS\_v2.md.]

#pagebreak()
= References

#bibliography("../synthesis/REFERENCES_v3.bib", style: "chicago-author-date", title: none)
