## Appendix: Files and tools used

### External tools

| Tool | Version | Path |
|------|---------|------|
| wfmash | v0.23.0-41-gb5f0ff1c | `/moosefs/guarracino/pggb_wfmash023/wfmash/build/bin/wfmash` |
| impg | git build (commit 5b96025) | `/moosefs/guarracino/git/impg/target/release/impg` |
| pggb | matching wfmash 0.23 | `/moosefs/guarracino/pggb_wfmash023/pggb/pggb` |
| odgi | bundled with pggb | `/moosefs/guarracino/pggb_wfmash023/smoothxg/deps/odgi/bin/odgi` |
| samtools | conda hicexplorer 3.7.4 | `/moosefs/guarracino/condatools/hicexplorer/3.7.4/bin/samtools` |
| bedtools | conda hicexplorer 3.7.4 | `/moosefs/guarracino/condatools/hicexplorer/3.7.4/bin/bedtools` |
| bgzip | conda hicexplorer 3.7.4 | `/moosefs/guarracino/condatools/hicexplorer/3.7.4/bin/bgzip` |
| pigz | system | `/usr/bin/pigz` |
| seqtk | git build | `/moosefs/guarracino/git/seqtk/seqtk` |
| Rscript | guix | `/home/guarracino/.guix-profile/bin/Rscript` |
| python3 | conda hicexplorer 3.7.4 | `/moosefs/guarracino/condatools/hicexplorer/3.7.4/bin/python3` |

*Scripts, input data, intermediate files, and output directories are listed in the "Files and scripts" section at the end of each sub-report file.*

---

## References

- Ambrosini A, Paul S, Hu S, Riethman H (2007). Human subtelomeric duplicon structure and organization. *Genome Biology* 8:R151. — Duplicon module classification; bimodal identity distribution (91%/98% peaks); subterminal vs subtelomere-only block distinction; one-copy (TTAGGG)n-adjacent regions at chr7_q, chr8_q, chr11_q, chr12_q, chr18_q, chrX_p/chrY_p.
- Brown WRA, MacKinnon PJ, Villasanté A, et al. (1990). Structure and polymorphism of human telomere-associated DNA. *Cell* 63:119–132. — Original characterization of TAR1 (Telomere-Associated Repeat 1) as a subtelomeric repeat element.
- Flint J, Bates GP, Clark K, et al. (1997). Sequence comparison of human and yeast telomeres identifies structurally distinct subtelomeric domains. *Human Molecular Genetics* 6:1305–1313. — Proposed two-domain subtelomeric model: distal domain (short blocks shared with many chromosome ends) and proximal domain (longer blocks shared with few ends), separated by degenerate (TTAGGG)n tracts; characterized at chr4p, chr16p, chr22q.
- Francis BA, et al. (2025). Complete genome assemblies of two mouse subspecies reveal structural diversity of telomeres and centromeres. *Nature Genetics* 57:2852–2862. — First T2T assemblies for C57BL/6J and CAST/EiJ; B6 adds 208 Mb, CAST 247 Mb vs GRCm39; conserved L1-LINE + TLC subtelomeric architecture in B6; heterogeneous CAST subtelomeres.
- Gershman A, Sauria MEG, Guitart X, et al. (2022). Epigenetic patterns in a complete human genome. *Science* 376:eabj5089. — ENCODE CTCF ChIP-seq realignment to T2T-CHM13; CTCF enrichment at TAR loci.
- Gonzalez IL & Sylvester JE (1995). Complete sequence of the 43-kb human ribosomal DNA repeat: analysis of the intergenic spacer. *Genomics* 27:320–328. — rDNA sequence variants spread across all five acrocentric chromosomes via inter-chromosomal gene conversion.
- Lemmers RJLF, van der Vliet PJ, Klooster R, et al. (2010). A unifying genetic model for facioscapulohumeral muscular dystrophy. *Science* 329:1650–1653. — Permissive 4qA haplotype produces stable DUX4 mRNA via polyadenylation signal.
- Ottaviani A, Rival-Gervier S, Boussouar A, et al. (2009). The D4Z4 macrosatellite repeat acts as a CTCF and A-type lamins-dependent insulator in facio-scapulo-humeral dystrophy. *PLoS Genetics* 5:e1000394. — CTCF binding within D4Z4 repeat units; lamin A/C-dependent insulator function.
- Linardopoulou EV, Williams EM, Fan Y, et al. (2005). Human subtelomeres are hot spots of interchromosomal recombination and segmental duplication. *Nature* 437:94–100. — Paralogy map of 41 blocks across 33 subtelomeres; subtelomeric interchromosomal duplication/transfer rate >60-fold higher than point mutation or retrotransposon insertion rates.
- Masny PS, Bengtsson U, Chung SA, et al. (2004). Localization of 4q35.2 to the nuclear periphery: is FSHD a nuclear envelope disease? *Human Molecular Genetics* 13:1857–1871. — 4q35.2 peripheral localization via lamin A/C; sequences proximal to D4Z4 mediate nuclear envelope interaction.
- Mefford HC & Trask BJ (2002). The complex structure and dynamic evolution of human subtelomeres. *Nature Reviews Genetics* 3:91–102. — Foundational review; f7501 block distribution; ~20% chr4_q/chr10_q translocation prevalence; patchwork duplicon architecture; telomere position effect.
- Patel L, Kang R, Rosenberg SC, et al. (2019). Dynamic reorganization of the genome shapes the recombination landscape in meiotic prophase. *Nature Structural & Molecular Biology* 26:164–174. — Mouse meiotic Hi-C; chromosome end clustering.
- Riethman H, Ambrosini A, Castaneda C, et al. (2004). Mapping and initial analysis of human subtelomeric sequence assemblies. *Genome Research* 14:18–28. — Subtelomeric assembly and characterization; 80% of the most distal 100 kb consists of shared duplicated blocks.
- Rouyer F, Simmler MC, Johnsson C, et al. (1986). A gradient of sex linkage in the pseudoautosomal region of the human sex chromosomes. *Nature* 319:291–295. — Obligate crossover in PAR1 during male meiosis.
- Stout K, van der Maarel S, Frants RR, et al. (1999). Somatic pairing between subtelomeric chromosome regions: implications for human genetic disease? *Chromosome Research* 7:323–329. — Interphase chr4_q/chr10_q clustering by FISH.
- Tan L, Xing D, Chang CH, et al. (2018). Three-dimensional genome structures of single diploid human cells. *Science* 361:924–928. — Dip-C method; GM12878 and PBMC single-cell 3D structures.
- Traag VA, Waltman L, van Eck NJ (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports* 9:5233. — Leiden community detection algorithm.
- Lalli JL, Bortvin AN, McCoy RC, Werling DM (2025). A T2T-CHM13 recombination map and globally diverse haplotype reference panel improves phasing and imputation. *bioRxiv* 2025.02.24.639687. — T2T-CHM13 recombination maps (preprint).
- Xu H, et al. (2025). Three-dimensional genome organization of human sperm at single-cell resolution. *Nature Communications* 16:3805. — Sperm single-cell 3D genome structures (20 cells).
- Zuo W, Chen G, Gao Z, et al. (2021). Stage-resolved Hi-C analyses reveal meiotic chromosome organizational features influencing homolog alignment. *Nature Communications* 12:5827. — Meiotic chromosome end alignment extends up to 20% of chromosome length; LINC complex force transmission.
