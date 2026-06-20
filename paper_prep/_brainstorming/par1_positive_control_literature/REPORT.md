# PAR1 positive-control interpretation for Fig5 pedigree/sweepGA analysis

## executive_summary

Yes: a paternal chrX/chrY event whose breakpoint and inheritance are confined to PAR1 is a defensible positive control for a pedigree/assembly method that aims to detect inter-chromosomal subtelomeric recombination. The key reason is not that PAR1 is a PHR-like autosomal subtelomere; it is that human male meiosis normally requires an X-Y crossover in the short-arm pseudoautosomal region for sex-chromosome pairing, chiasma formation, and proper segregation. Human PAR1 is therefore one of the best-established biological settings in which a paternal inter-chromosomal subtelomeric exchange is expected.

The strongest wording is: "the PAN027 chrX/chrY PAR1 call provides an internal positive-control class for paternal inter-chromosomal recombination detection." This is stronger than a generic plausibility argument, because PAR1 recombination has been measured by pedigree mapping, sperm typing, and population-genetic recombination analyses, and reduced or absent pseudoautosomal recombination is associated with X-Y nondisjunction and male meiotic defects.

The language should remain cautious. The method still needs the same alignment, phasing, inheritance, and breakpoint checks required for autosomal PHR candidates. PAR1 should not be framed as evidence that autosomal PHR recombination is mechanistically analogous to the normal sex-chromosome crossover. PAR1 is a homologous pseudoautosomal interval with a dedicated meiotic role; autosomal PHR candidates involve subtelomeric homology blocks whose meiotic behavior is the biological question.

PAR2 should not be used interchangeably with PAR1 in this framing. Human PAR2 is a much smaller Xq/Yq pseudoautosomal interval that can recombine and contains recombination hotspots, but the obligate male sex-chromosome crossover is conventionally associated with PAR1, not PAR2. PAR2 is best described as pseudoautosomal but not the canonical site of the required X-Y crossover.

## evidence_by_topic

### 1. Human PAR1 is the canonical site of near-obligatory male X-Y recombination

Human males have one X and one Y chromosome. Most of the X and Y are nonhomologous, so the homologous pseudoautosomal regions provide the main substrate for meiotic pairing and exchange. PAR1 lies at the short-arm tips, is approximately 2.6-2.7 Mb in humans, and carries the canonical male X-Y crossover. PAR1 recombination is therefore unusually intense per base pair because at least one crossover must usually be accommodated in a physically short interval.

Primary evidence:

- Rouyer et al. 1986 identified a gradient of sex linkage across the human pseudoautosomal region, establishing the short-arm X/Y pseudoautosomal region as a recombining interval rather than ordinary sex-linked sequence. DOI: <https://doi.org/10.1038/319291a0>.
- Schmitt et al. 1994 used single-sperm typing across PAR1 and treated the region as the male X-Y recombination interval. URL: <https://pubmed.ncbi.nlm.nih.gov/7915881/>.
- May et al. 2002 mapped extremely localized crossover activity and rapid linkage disequilibrium decay around the PAR1 gene SHOX. DOI: <https://doi.org/10.1038/ng918>.
- Hinch et al. 2014 integrated sperm-typing, pedigree maps, linkage disequilibrium, GC content, and PRDM9 ChIP-seq. They state directly that a PAR crossover is essential for proper X/Y disjunction in male meiosis and that PAR1 is the human "obligatory crossover" PAR, with an approximately 17-fold elevation in male crossover rate relative to the genome-wide average. DOI: <https://doi.org/10.1371/journal.pgen.1004503>.
- Bergman and Schierup 2022 used population-genomic and comparative data and summarized PAR1 as a 2.7-Mb telomeric interval with a crucial role in proper X-Y segregation during male meiosis and extreme recombination-associated evolution. DOI: <https://doi.org/10.1186/s13059-022-02784-x>.

Interpretation for Fig5:

- A paternal chrX/chrY PAR1 recombinant haplotype is expected biology, not a novel claim about an unusual PHR mechanism.
- It tests whether the pedigree/assembly procedure can recover a real inter-chromosomal subtelomeric exchange when such exchange has a known meiotic context.
- It is most defensible as a "positive-control class" or "internal positive-control event" rather than as a formal spike-in control, because the event is observed in the same dataset rather than experimentally introduced.

### 2. PAR1 versus PAR2

PAR1 and PAR2 are both pseudoautosomal in the sense that X and Y share homologous sequence and can exchange genetic material. They are not equivalent for the positive-control argument.

PAR1:

- Location: Xp/Yp telomeric region.
- Size: roughly 2.6-2.7 Mb in standard references, with a boundary that has been refined and may vary in extended-PAR haplotypes.
- Biological framing: canonical site of the required male X-Y crossover.
- Recombination: very high male recombination rate per Mb; numerous pedigree/sperm/population studies.

PAR2:

- Location: Xq/Yq telomeric region.
- Size: roughly 320 kb in standard references.
- Biological framing: pseudoautosomal and capable of recombination, but not conventionally the required X-Y crossover site.
- Recombination: recombination and hotspots have been reported, but it is not the main basis for sex-chromosome segregation assurance in humans.

Primary evidence:

- Freije et al. 1992 identified the second pseudoautosomal region near the Xq/Yq telomeres. DOI: <https://doi.org/10.1126/science.1465614>.
- Charchar et al. 2003 analyzed the fully sequenced Xq/Yq PAR2 and described complex evolutionary and functional structure, reinforcing that PAR2 has a distinct origin/behavior from PAR1. DOI: <https://doi.org/10.1101/gr.390503>.
- Sarbajna et al. 2012 identified a major recombination hotspot in the Xq/Yq PAR2 region and characterized crossover/gene-conversion products in sperm. DOI: <https://doi.org/10.1093/hmg/dds019>.
- Monteiro et al. 2021 compared PAR1 and PAR2 evolutionary dynamics using 1000 Genomes data and reported higher recombination frequency in PAR1 than PAR2 with uneven recombination in both regions. DOI: <https://doi.org/10.1371/journal.pgen.1009532>.
- Flaquer et al. 2009 built sex-specific maps for PAR1 and PAR2 and estimated much larger male PAR1 genetic lengths than female PAR1 lengths; their work is useful for distinguishing sex- and region-specific pseudoautosomal maps. DOI: <https://doi.org/10.1159/000224639>.

Interpretation for Fig5:

- If PAN027 is chrX/chrY PAR1, call it PAR1 specifically.
- Do not write "PAR recombination is obligate" without specifying male PAR1, because PAR2 can recombine but is not the canonical obligate X-Y crossover interval.
- If a candidate were in PAR2, it would still be biologically plausible pseudoautosomal recombination, but it would be weaker as a positive control for the obligatory sex-chromosome crossover.

### 3. Male versus female meiosis

The positive-control logic is male-specific.

In male meiosis, X and Y are heteromorphic and share homology mainly in pseudoautosomal regions. An X-Y crossover produces the chiasma needed for normal segregation of the sex chromosomes. PAR1 is small, so the male crossover rate per base pair is very high.

In female meiosis, the two X chromosomes are homologous across most of their length. PAR1 is not uniquely required to make X chromosomes pair, synapse, and segregate. Female recombination can occur in PAR1, and some studies report elevated recombination near the Xp telomere, but this does not carry the same "obligate X-Y crossover" interpretation.

Primary evidence:

- Henke et al. 1993 reported high female recombination at the Xp telomere in PAR1, showing that PAR1 is not male-exclusive. DOI: <https://doi.org/10.1016/S0888-7543(11)80003-0>.
- Flaquer et al. 2009 produced sex-specific maps and summarized published male PAR1 map-length estimates in the 26-54 cM range versus female estimates around 2.8-6 cM. DOI: <https://doi.org/10.1159/000224639>.
- Flaquer et al. 2008 is a review, not primary evidence, but it concisely states that PARs pair and recombine like autosomes while PAR1 recombination activity differs strongly between sexes. DOI: <https://doi.org/10.1038/ejhg.2008.63>.

Interpretation for Fig5:

- The PAN027 event is especially useful if the parent of origin is paternal.
- A maternal PAR1 crossover would be real X-X recombination but would not be the same positive-control class for the male X-Y segregation crossover.

### 4. Mechanism: pairing, synapsis, homologous exchange, crossover assurance, and segregation

The standard mechanism is:

1. X and Y share sequence homology in PAR1.
2. During male meiotic prophase, the X and Y must pair/synapse at this homologous region.
3. Recombination machinery creates and repairs double-strand breaks; at least one repair outcome must mature as a crossover.
4. The crossover creates a chiasma linking the X and Y bivalent.
5. The chiasma supports proper orientation and segregation at meiosis I.

Human-specific data are strong for the phenotype and the recombination landscape. Some mechanistic details come from mouse because meiotic prophase can be studied more directly there; those should be labeled as model-organism evidence when used.

Primary evidence in humans:

- Gabriel-Robez et al. 1990 described infertile men with X;Y translocations deleting the pseudoautosomal region and lacking sex-chromosome pairing at pachytene. DOI: <https://doi.org/10.1159/000132951>.
- Mohandas et al. 1992 studied a man with a distal Xp deletion including the pseudoautosomal region and performed meiotic analyses supporting a role for PAR in male sex-chromosome pairing. URL: <https://pubmed.ncbi.nlm.nih.gov/1496984/> and <https://pmc.ncbi.nlm.nih.gov/articles/PMC1682713/>.
- Hassold et al. 1991 found that human XY nondisjunction was associated with diminished recombination in the pseudoautosomal region. URL: <https://pubmed.ncbi.nlm.nih.gov/1867189/> and <https://pmc.ncbi.nlm.nih.gov/articles/PMC1683286/>.
- Shi et al. 2001 used single-sperm typing and found reduced recombination associated with aneuploid 24,XY sperm. DOI: <https://doi.org/10.1002/1096-8628(20010215)99:1%3C34::AID-AJMG1106%3E3.0.CO;2-D>.

Primary mechanistic/model evidence:

- Kauppi et al. 2011 showed in mouse that the XY pseudoautosomal region has distinct recombination properties crucial for male meiosis, including mechanisms that help ensure double-strand break formation in the small PAR. This is mouse evidence but supports the general crossover-assurance concept. DOI: <https://doi.org/10.1126/science.1195774>.
- Acquaviva et al. 2020 showed mechanisms ensuring meiotic DNA-break formation in the mouse PAR. Again, this is mouse evidence and should not be over-translated to human without qualification. DOI: <https://doi.org/10.1038/s41586-020-2327-4>.

Interpretation for Fig5:

- Mechanistic language can mention "X-Y pairing and segregation" confidently for human male PAR1.
- More detailed statements about repeat arrays, axis structure, or exact crossover-assurance machinery should be attributed to mouse/model studies unless a human paper is specifically cited.

### 5. PRDM9, hotspots, and recombination-rate elevation

PAR1 is not just broadly recombinogenic; it contains hotspots, and multiple studies connect human PAR1 activity to PRDM9 or PRDM9-associated hotspot biology.

Primary evidence:

- May et al. 2002 showed crossover clustering and rapid LD decay in the PAR1 SHOX region. DOI: <https://doi.org/10.1038/ng918>.
- Hinch et al. 2014 found evidence that PRDM9 localizes recombination peaks in human PAR1, in contrast to the mouse PAR where PRDM9-independent mechanisms have a stronger role. DOI: <https://doi.org/10.1371/journal.pgen.1004503>.
- Poriswanish et al. 2018 studied recombination hotspots in an extended human pseudoautosomal domain using double-strand break maps and sperm-based crossover assays. They found PRDM9-associated DSB clusters and directly measured sperm crossovers in an ePAR interval, concluding that ePAR likely contributes to the critical PAR1 crossover function. DOI: <https://doi.org/10.1371/journal.pgen.1007680>.
- Baudat et al. 2010 and Myers et al. 2010 are general primary sources establishing PRDM9 as a major determinant of human recombination hotspots and hotspot motif evolution, not PAR1-specific evidence. Baudat DOI: <https://doi.org/10.1126/science.1183439>. Myers DOI: <https://doi.org/10.1126/science.1182363>.

Interpretation for Fig5:

- It is safe to say that human PAR1 has a high male recombination rate and contains hotspots.
- It is safe to say human PAR1 hotspot localization is at least partly PRDM9-associated, citing Hinch et al. and Poriswanish et al.
- Avoid saying PAR1 is PRDM9-independent in humans. That would be misleading; mouse PAR biology differs from human PAR1 in this respect.

### 6. Failure modes and clinical consequences

The most relevant failure mode for this report is not SHOX-related growth disease per se; it is failure of sex-chromosome pairing/recombination leading to male meiotic arrest/infertility or XY nondisjunction.

Human evidence:

- Gabriel-Robez et al. 1990: X;Y translocation carriers with PAR deletion lacked sex-chromosome pairing at pachytene and were infertile. DOI: <https://doi.org/10.1159/000132951>.
- Mohandas et al. 1992: distal Xp/PAR deletion case with meiotic studies supported the role of PAR in sex-chromosome pairing. URL: <https://pubmed.ncbi.nlm.nih.gov/1496984/>.
- Hassold et al. 1991: XY nondisjunction was associated with diminished recombination in the pseudoautosomal region. URL: <https://pubmed.ncbi.nlm.nih.gov/1867189/>.
- Shi et al. 2001: reduced X-Y recombination was associated with production of aneuploid 24,XY sperm. DOI: <https://doi.org/10.1002/1096-8628(20010215)99:1%3C34::AID-AJMG1106%3E3.0.CO;2-D>.

Secondary/review context:

- Hall, Hunt, and Hassold 2006 reviewed how meiotic errors cause sex-chromosome aneuploidy. This is useful background but should not be used as primary evidence for the specific PAN027 framing. DOI: <https://doi.org/10.1016/j.gde.2006.04.011>.
- Blaschke and Rappold 2006 reviewed PARs, SHOX, and disease; useful for separating PAR1 disease genetics from the recombination-positive-control argument. DOI: <https://doi.org/10.1016/j.gde.2006.04.004>.

Interpretation for Fig5:

- Clinical/failure evidence strengthens the claim that male PAR1 recombination is functionally important.
- Do not imply that a normal PAR1 crossover is itself pathological. The clinical literature concerns absent/reduced/abnormal recombination or rearranged PARs, not ordinary paternal PAR1 exchange.

### 7. What counts as a legitimate positive-control signal in pedigree/assembly alignments

A PAR1 event is a legitimate positive-control signal only if the observed pattern is consistent with a real paternal X-Y pseudoautosomal crossover rather than an alignment, phasing, sample, or assembly artifact.

Minimum evidence to treat as positive control:

- Parent of origin: paternal transmission is supported by the pedigree and haplotype pattern.
- Coordinates: the breakpoint and exchanged tract lie within PAR1 or a well-supported extended PAR1 haplotype, not in X-specific/Y-specific sequence mislabeled as pseudoautosomal.
- Reciprocal logic: the inherited haplotype switches between paternal X-like and Y-like ancestry across PAR1 in a way compatible with one meiotic exchange.
- Assembly support: contig-level alignments span the relevant interval without relying only on short, repetitive, or telomere-proximal ambiguous alignments.
- Phasing support: flanking markers and sample relationships rule out sample swap, phase switch, paralogous mapping, or contamination.
- Sex-chromosome copy-state support: the event is not better explained by sex-chromosome aneuploidy, large deletion/duplication, or assembly collapse.
- Boundary awareness: use current PAR1 boundary definitions where possible, and report if the inferred breakpoint falls near a debated/polymorphic pseudoautosomal boundary.

Useful recent boundary/variation context:

- Mensah et al. 2014 showed PAR1 length polymorphism in the human population, meaning the pseudoautosomal boundary is not perfectly static across all haplotypes. DOI: <https://doi.org/10.1371/journal.pgen.1004578>.
- Bellott et al. 2024 refined the human pseudoautosomal boundary using complete X and Y assemblies. This is current boundary context, not necessary for the general biological argument, but useful if the event is near the PAR boundary. DOI: <https://doi.org/10.1016/j.ajhg.2024.09.005>.
- Poriswanish et al. 2018 showed recombination in extended PAR1 haplotypes, reinforcing that boundary polymorphism can matter. DOI: <https://doi.org/10.1371/journal.pgen.1007680>.

Interpretation for Fig5:

- If the PAN027 event satisfies these checks, it can be described as "an expected paternal PAR1 recombination event recovered by the pipeline."
- If support is weaker or breakpoint resolution is broad, use "PAR1-consistent" or "candidate paternal PAR1 exchange" rather than definitive language.

### 8. Relationship to autosomal PHR candidate events

The logic should be asymmetric:

- PAR1: known positive-control class for paternal inter-chromosomal subtelomeric exchange because male X-Y recombination in PAR1 is normal and functionally required.
- Autosomal PHR candidates: discovery/target class where inter-chromosomal subtelomeric homology may mediate recombination-like inheritance patterns, but the mechanism, rate, and recurrence are not established by the PAR1 observation.

Good framing:

- "PAR1 validates the detection logic for a known inter-chromosomal subtelomeric recombination context; autosomal PHR calls remain candidate events requiring independent evidence."
- "The PAR1 event is a control for assay sensitivity/specificity to a known class of recombination, not a mechanistic proxy for autosomal PHR exchange."

Less good framing:

- "Because PAR1 recombines, autosomal PHRs also recombine."
- "PAR1 proves subtelomeric inter-chromosomal recombination is generally common."
- "The PAR1 event validates all candidate events."

## recommended_manuscript_language

1. "As an internal positive-control class, the pedigree analysis recovered a paternal chrX/chrY exchange in PAR1, the pseudoautosomal interval in which the male X and Y normally recombine to support sex-chromosome segregation."

2. "The PAR1 call is biologically expected, because human male meiosis concentrates the required X-Y crossover in the short-arm pseudoautosomal region; we therefore used it to calibrate interpretation of candidate autosomal PHR exchanges rather than as evidence for a shared mechanism."

3. "A paternal PAR1 recombination event provides a positive-control benchmark for detecting inter-chromosomal subtelomeric exchange in the pedigree data, while autosomal PHR events remain candidate observations requiring locus-specific support."

4. Caption phrase: "Positive-control paternal chrX/chrY PAR1 exchange, shown alongside candidate autosomal PHR recombination events."

## do_not_overclaim

- Do not say "PAR recombination is obligate" without specifying male PAR1. PAR2 can recombine, but the canonical obligate sex-chromosome crossover is PAR1.
- Do not imply that every male meiosis must have an observable PAR1 crossover in the child-level assembly data. The biological expectation is at the meiotic level; detectability depends on marker informativeness, assembly continuity, breakpoint location, and phasing.
- Do not use PAR1 as mechanistic proof for autosomal PHR recombination. PAR1 is homologous pseudoautosomal sequence with a dedicated sex-chromosome segregation role.
- Do not say "PRDM9-independent PAR1 recombination" for humans. Human PAR1 hotspot evidence supports PRDM9-associated localization; PRDM9-independent mechanisms are better established in mouse PAR biology.
- Do not ignore PAR1 boundary polymorphism if the candidate lies near the pseudoautosomal boundary. Use "within PAR1 under [reference/boundary definition]" or "PAR1/ePAR-consistent" when necessary.
- Do not frame a normal paternal PAR1 crossover as a disease event. Disease/failure studies support functional importance by showing consequences of absent/reduced/abnormal recombination.
- Do not call the PAN027 event a "spike-in" or "experimental positive control." It is an observed internal biological positive-control class.

## bottom_line_for_fig5

Use PAR1 as a clear but carefully bounded positive control. The best manuscript logic is:

"The paternal chrX/chrY PAR1 event is not part of the autosomal PHR discovery claim. It is a known, expected inter-chromosomal subtelomeric recombination context that shows the pedigree/assembly analysis can recover a biologically canonical exchange. We therefore present it as an internal positive-control benchmark and interpret autosomal PHR calls separately as candidate events."
