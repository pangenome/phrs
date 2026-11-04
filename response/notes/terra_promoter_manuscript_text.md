# The subtelomeric CpG island / TERRA promoter: manuscript-ready text

Status: draft text for insertion, with citations verified against primary sources (see §5 for the
verification note per entry). Nothing here has been applied to `submission/paper.tex`.

Two claims carried in earlier notes could not be verified on this pass and are **not** used below;
they are listed as `AUTHOR-CHECK` in §7.

---

## 1. Methods paragraph

> **\subsection*{The terminal CpG island, GC content and annotatability}**
>
> We characterised the sequence immediately inside each telomere on the same telomere-first 500 kb
> flanks used for PHR calling, cut into 2 kb bins measured from the junction. Each bin was scored for
> GC content and for CpG observed/expected (O/E, the count of CpG dinucleotides divided by the count
> expected from the bin's own G and C frequencies), and we called a CpG island where O/E $>0.6$ and
> GC $>55\%$ simultaneously. Both thresholds were fixed before the distance profile was inspected;
> genome-wide CpG O/E in these assemblies is 0.2--0.3, so the threshold leaves an order of magnitude
> of headroom between ordinary sequence and an island. To establish that the tip windows are
> subtelomeric rather than telomeric, we counted (TTAGGG)+(CCCTAA) units per kb in the same bins, and
> to exclude low-complexity sequence as the cause of any anomaly we computed the fraction of each bin
> occupied by 6-mers that recur within the bin. Statistics were taken over 30--120 sampled
> telomere-anchored windows per measurement (30 for CpG O/E, GC and low complexity; 40 for telomeric
> repeat content). Gene transcription start sites were taken from the CHM13v2.0 RefSeq/Liftoff v5.2
> annotation and reduced to a single distance per chromosome end: the gap between the
> telomere-proximal end of the flank and the nearest annotated TSS. Alignability was measured on
> same-arm (allelic) alignment records per 2 kb bin, assigning a record to a bin only when the
> midpoint of *both* the query and the target span fell in that bin, exactly as in the distance
> profile; we report the number of records, the number of distinct individuals contributing, and the
> median block length.

## 2. Results passage

> **\subsection*{The exchange junction is a copy-variable TERRA promoter}**
>
> The first two kilobases inside the telomere are not ordinary subtelomeric sequence. We observe a
> CpG island there (CpG O/E 1.02, GC 62.3\%), against 40.9\% GC at 2--5 kb and 40.2\% at 5--10 kb and
> a genome-wide CpG O/E of 0.2--0.3. The window is not telomere: it carries about one
> (TTAGGG)/(CCCTAA) unit per kb, where a genuine repeat array runs 150--170. Nor is it
> low-complexity sequence, which is least abundant at the tip of anywhere in the subtelomere (32.8\%
> at 0--2 kb against 82.4\% at 20--50 kb, where the satellite and duplicon arrays sit). The island is
> present at 88\% of ends lacking an annotated tip gene and 95\% of ends carrying one, so it is not a
> by-product of gene annotation.
>
> This is a known element. The human terminal repeat region is built, centromere to telomere, from a
> 61 bp, a 29 bp and a 37 bp tandem repeat, and the 29 bp and 37 bp elements together form a DNA
> island rich in CpG dinucleotides that functions as the promoter driving TERRA
> \cite{Brown1990,Nergadze2009,Farnung2010}; the 29 bp and 37 bp elements alone are sufficient to
> drive a reporter, TERRA start sites on XqYq lie 27 nt downstream of the last 37 bp repeat and about
> 250 bp upstream of the telomeric tract \cite{Nergadze2009}, and the element is bound in vivo by
> active RNA polymerase II \cite{Nergadze2009} and sits next to a CTCF/cohesin site within 1--2 kb of
> the TTAGGG tract \cite{Deng2012}. At least 20 chromosome ends were found to carry 61-29-37 repeats
> in different tandem copy numbers in the clone-based era \cite{Brown1990,Nergadze2009}; in
> T2T-CHM13v2.0, 39 of 46 subtelomeres carry a TERRA promoter, classed by repeat organisation
> \cite{Santagostino2025}, and an older 3.6 kb repeat immediately adjacent to the tract likewise
> contains a conserved CpG-rich region whose copy number varies among telomeres \cite{Cheng1991}. Our
> measurement therefore confirms rather than discovers: the CpG island at the exchange junction is
> the CpG-island class of TERRA promoter, and its GC content sits at the upper extreme of a
> subtelomeric range already known to run from 62.5\% at 1p to 37.5\% at 3p against a genome-wide
> 41.6\% \cite{Riethman2004}.
>
> The genes annotated at the tip are a separate matter and a weak one. Seven of 46 ends have an
> annotated TSS within 2 kb of the tract and 14 of 46 within 5 kb; the DDX11L family supplies 8 of
> the 15 genes at 1.7--2.5 kb, with LOC102723376, LINC01986 (749 bp) and LINC02684 (859 bp) the only
> non-DDX11L loci inside 1 kb. Their presence does not explain the exchange regime: arms with a tip
> gene have a median paralog/allelic divergence ratio of 3.06 against 3.82 for arms without one
> (Mann--Whitney $p = 0.85$). Whatever sets how homogenised an arm is, it is not which gene happens
> to be annotated at its end.
>
> The promoter is also the reason the terminal few kilobases cannot be measured by sequence
> comparison, which is a stronger statement than irreproducibility. Same-arm alignments collapse
> there: 41--60 records of median block 2,728--3,007 bp in the 2--4 kb bin against 4,933--5,283
> records of about 50,000 bp in the 28--30 kb bin. A tandem array whose copy number and phase differ
> between haplotypes produces exactly that — short, phase-shifted, partly paralogous matches instead
> of one homologous stretch. The paralog/allelic ratio in that zone is therefore undefined rather
> than noisy, and we exclude it rather than report it (Extended Data Fig.~\ref{fig:ed6}). Outside the
> excluded zone the ratio behaves regularly: about 1.3 within the first 8 kb, rising through 2--3 at
> 20--60 kb to 6--7 beyond 300 kb.

## 3. Discussion passage

> New paragraph to follow the self-reinforcing-model paragraph in `\subsection*{Discussion}`
> ("This compartment is also internally organized. … We see the same logic operating on other
> functional modules outside the NORs."), with the last paragraph below added at the end of
> "This concerted process operates at a characteristic scale."
>
> The element at the junction sharpens this picture in one respect and limits it in another. The
> CpG-island TERRA promoter is not a unique sequence that happens to lie near a telomere: it is a
> dispersed repeat family present, in different tandem copy numbers and in different arrangements, at
> most chromosome ends \cite{Brown1990,Nergadze2009,Santagostino2025}, with deep homology to a family
> of segmental duplications that also recurs at interstitial telomeric sites
> \cite{Santagostino2025}. The immediate substrate of inter-arm exchange is therefore itself a shared,
> copy-variable repeat — the same kind of duplicon that the rest of the compartment is built from, and
> one whose variation between haplotypes is already documented. That makes it a plausible participant
> in the exchange we quantify, although it does not demonstrate that exchange occurs within it.
>
> The same element bounds what the divergence profile can say. Because copy number and phase of the
> promoter array differ between haplotypes, same-arm alignments inside the first few kilobases are
> fragmentary and partly paralogous, and no orthologous comparison survives there. We therefore do not
> extrapolate the ratio to the telomere itself and do not claim it approaches one at the junction; the
> measured statement is that paralogues are only about 1.3-fold more divergent than alleles just
> outside the excluded zone, against roughly six-fold in the far interior, so the terminal sequence is
> the most nearly interchangeable in the genome and the promoter array is where that measurement
> stops. Reading the ratio as a ratio of genealogical depths rather than a recombination rate, the
> average subtelomeric block exchanges at something like one-sixth of the obligate rate while the
> most cohesive blocks approach it \cite{Dover1982,Hastings2009}.
>
> One mechanism we can exclude is a role for the genes at the tip. Whether an arm carries an annotated
> subtelomeric gene within 5 kb of its tract — overwhelmingly a DDX11L pseudogene — has no bearing on
> how homogenised that arm is ($p = 0.85$). The homogenisation is a property of the shared repeat
> architecture, not of the gene complement that happens to be annotated alongside it.

## 4. Replacement caption text for the hatched zone in Extended Data Fig. 6

The hatched zone is currently described only as "unmeasurable". Replace that wording with:

> Shaded zone, 0--6 kb: the TERRA promoter array (61-29-37 repeats and its CpG island)
> \cite{Nergadze2009,Santagostino2025}. This region is a copy-variable tandem repeat family shared
> between non-homologous chromosome ends, and its copy number and phase differ between haplotypes, so
> same-arm alignments there are short and partly paralogous and the paralog/allelic ratio is
> undefined rather than imprecise; it is shown for orientation only and is not a measurement.

Full caption skeleton, if ED Fig. 6 is being written in full (which it must be — the legends section of
`paper.tex` currently ends at ED Fig. 4):

> **\textbf{Divergence ratio against distance from the telomere, one vote per chromosome end.}**
> Each light-grey curve is one individual chromosome end (individual $\times$ arm), plotted from the
> telomere–subtelomere junction outward; heavy curves are per-arm medians and the cohort median.
> Ordinate is the paralog/allelic divergence ratio, the median inter-chromosomal paralogous divergence
> of a 2 kb distance bin divided by the same-arm (allelic) divergence of the same bin, with a record
> assigned only when the midpoint of both the query and the target span fall in that bin. The ratio
> rises monotonically with distance, from about 1.3 within the first 8 kb to 6--7 beyond 300 kb.
> **Shaded zone, 0--6 kb:** the TERRA promoter array (61-29-37 repeats and its CpG island)
> \cite{Nergadze2009,Santagostino2025}. This region is a copy-variable tandem repeat family shared
> between non-homologous chromosome ends, and its copy number and phase differ between haplotypes, so
> same-arm alignments there are short and partly paralogous and the ratio is undefined rather than
> imprecise; it is shown for orientation only and is not a measurement. Hatched width and position are
> set by the allelic-divergence rule in Methods (bins whose allelic divergence exceeds 0.4\%, four times
> the interior baseline, are discarded), not by hand.

## 5. BibTeX entries

Keys already in `submission/bibliography.bib` are marked **EXISTS — reuse**; entries are repeated only
so this file is self-contained. Do not duplicate a key.

```bibtex
% EXISTS -- reuse. Fix the first-author initials: PubMed records "Brown WR",
% the bib currently has "Brown, W. R. A." (same paper). Wilkie1991 has the
% same slip for the same author.
@article{Brown1990,
  author  = {Brown, W. R. and MacKinnon, P. J. and Villasant{\'e}, A. and
             Spurr, N. and Buckle, V. J. and Dobson, M. J.},
  title   = {Structure and polymorphism of human telomere-associated {DNA}},
  journal = {Cell},
  year    = {1990},
  volume  = {63},
  number  = {1},
  pages   = {119--132},
  doi     = {10.1016/0092-8674(90)90293-N},
  pmid    = {2208276}
}

% NEW. Primary TERRA-promoter paper. Verified from PMC2779677 full text:
% 61-29-37 architecture; CpG island over the 29+37 elements; 29+37 sufficient
% for an eGFP reporter; ~20 ends with variable tandem copy number; single TSS
% 27 nt downstream of the last 37 bp and ~250 nt from the telomeric tract;
% phospho-RNAPII binding; promoter CpG methylation represses TERRA.
@article{Nergadze2009,
  author  = {Nergadze, Solomon G. and Farnung, Benjamin O. and
             Wischnewski, Harry and Khoriauli, Lela and Vitelli, Valerio and
             Chawla, Raghav and Giulotto, Elena and Azzalin, Claus M.},
  title   = {{CpG}-island promoters drive transcription of human telomeres},
  journal = {RNA},
  year    = {2009},
  volume  = {15},
  number  = {12},
  pages   = {2186--2194},
  doi     = {10.1261/rna.1748309},
  pmid    = {19850908}
}

% NEW. Compact review statement of the same facts; also the source of the
% historical "approximately half of the human subtelomeres" figure.
% Verified from PMC3023574.
@article{Farnung2010,
  author  = {Farnung, Benjamin O. and Giulotto, Elena and Azzalin, Claus M.},
  title   = {Promoting transcription of chromosome ends},
  journal = {Transcription},
  year    = {2010},
  volume  = {1},
  number  = {3},
  pages   = {140--143},
  doi     = {10.4161/trns.1.3.13191},
  pmcid   = {PMC3023574}
}

% EXISTS -- reuse. T2T-CHM13 promoter census: 39/46 subtelomeres carry TERRA
% promoters; 106 intrachromosomal TERRA-like promoters; interstitial and
% subtelomeric promoters belong to a segmental-duplication family from common
% ancestral sequences. Best citation for "copy-variable dispersed repeat".
@article{Santagostino2025,
  author  = {Santagostino, Marco and Sola, Lorenzo and Cappelletti, Eleonora and
             Piras, Francesca M. and Gennari, Nicol{\`o} and Biundo, Marialaura and
             Nergadze, Solomon G. and Giulotto, Elena},
  title   = {{TERRA} transcripts and promoters from telomeric and interstitial sites},
  journal = {RNA},
  year    = {2025},
  volume  = {32},
  number  = {1},
  pages   = {97--112},
  doi     = {10.1261/rna.080790.125},
  pmid    = {41193243}
}

% NEW. CTCF/cohesin site within ~1-2 kb of the TTAGGG tract adjacent to a CpG
% island implicated in TERRA control; RNAPII extends from the CTCF site toward
% the tract. Abstract verified verbatim via Europe PMC (PMID 23010778).
@article{Deng2012,
  author  = {Deng, Zhong and Wang, Zhuo and Stong, Nicholas and
             Plasschaert, Robert and Moczan, Aliah and Chen, Horng-Shen and
             Hu, Sufeng and Wikramasinghe, Priyankara and Davuluri, Ramana V. and
             Bartolomei, Marisa S. and Riethman, Harold and Lieberman, Paul M.},
  title   = {A role for {CTCF} and cohesin in subtelomere chromatin organization,
             {TERRA} transcription, and telomere end protection},
  journal = {The EMBO Journal},
  year    = {2012},
  volume  = {31},
  number  = {21},
  pages   = {4165--4178},
  doi     = {10.1038/emboj.2012.266},
  pmid    = {23010778}
}

% NEW. Two types of TERRA promoter, only one DNA-methylation regulated; TERRA
% expressed from a large number of ends. Abstract verified via Europe PMC.
@article{Feretzaki2019,
  author  = {Feretzaki, Marianna and Renck Nunes, Patricia and Lingner, Joachim},
  title   = {Expression and differential regulation of human {TERRA} at several
             chromosome ends},
  journal = {RNA},
  year    = {2019},
  volume  = {25},
  number  = {11},
  pages   = {1470--1480},
  doi     = {10.1261/rna.072322.119},
  pmid    = {31350341}
}

% NEW. HST repeat: >3.6 kb immediately adjacent to the terminal (TTAGGG)n
% tract, one evolutionarily conserved CpG-rich region, copy number varying
% among telomeres. Abstract verified from PMC333545.
@article{Cheng1991,
  author  = {Cheng, Jan-Fang and Smith, Cassandra L. and Cantor, Charles R.},
  title   = {Structural and transcriptional analysis of a human subtelomeric repeat},
  journal = {Nucleic Acids Research},
  year    = {1991},
  volume  = {19},
  number  = {1},
  pages   = {149--154},
  doi     = {10.1093/nar/19.1.149},
  pmid    = {2011494}
}

% EXISTS -- reuse. Subtelomeric one-copy GC 47.9% vs 41.6% genome-wide;
% 1p 62.5%, 3p 37.5%. Verified verbatim from PMC314271.
@article{Riethman2004,
  author  = {Riethman, Harold and Ambrosini, Anthony and Castaneda, Carlos and
             Finklestein, Joshua and Hu, Xiao-Lin and Mudunuri, Uma and
             Paul, Sheila and Wei, Jeff},
  title   = {Mapping and initial analysis of human subtelomeric sequence assemblies},
  journal = {Genome Research},
  year    = {2004},
  volume  = {14},
  number  = {1},
  pages   = {18--28},
  doi     = {10.1101/gr.1245004},
  pmid    = {14707166}
}

% EXISTS -- reuse (optional second CTCF-organisation citation).
@article{Stong2014,
  author  = {Stong, Nicholas and Deng, Zhong and Gupta, Ravi and Hu, Sufeng and
             Paul, Sheila and Weiner, Aaron K. and Eichler, Evan E. and
             Graves, Tina and Fronick, Catrina C. and Courtney, Lucinda and
             Wilson, Richard K. and Lieberman, Paul M. and Davuluri, Ramana V. and
             Riethman, Harold},
  title   = {Subtelomeric {CTCF} and cohesin binding site organization using
             improved subtelomere assemblies and a novel annotation pipeline},
  journal = {Genome Research},
  year    = {2014},
  volume  = {24},
  number  = {6},
  pages   = {1039--1050},
  doi     = {10.1101/gr.166983.113}
}

% NEW. Companion paper. Metadata verified via Europe PMC (PMID 42922562).
% Cite ONLY for subtelomeric structural variation between haplotypes; their
% abstract does not mention TERRA and the full text was not retrievable --
% see AUTHOR-CHECK.
@article{Hastings2009,
  author  = {Hastings, P. J.},
  title   = {Mechanisms of ectopic gene conversion},
  journal = {Genes},
  year    = {2010},
  volume  = {1},
  pages   = {427--439},
  doi     = {10.3390/genes1030427},
  note    = {AUTHOR-CHECK: verify volume/pages against the publisher record before submission}
}

@article{Tan2026,
  author  = {Tan, Kar-Tong and Ong, Ryan Jun Xiang and Kee, Brandon Bing Rui and
             Yap, Russell Ker Han and Ng, Alicia Jun Ting and
             Rajapaksha, Gihan Kaushalya and Garrity-Janger, Max and Lin, Qiyu and
             Chong, Jerome Rui and Kyi, Cin Thet and {Human Pangenome Reference
             Consortium} and Meyerson, Matthew and Li, Heng},
  title   = {Widespread structural variations at human chromosome ends},
  journal = {bioRxiv},
  year    = {2026},
  pages   = {2026.09.03.748401},
  doi     = {10.64898/2026.09.03.748401},
  pmid    = {42922562},
  pmcid   = {PMC13684503}
}
```

Existing keys usable without new entries: `Dover1982`, `concerted_evolution_nahr_Hastings2009`,
`MeffordTrask2002`, `Linardopoulou2005`, `Nurk2022`, `hprc_hprcv2_2025`, `Ambrosini2007`.

Verification note: Nergadze 2009, Farnung 2010 and Riethman 2004 were read as full text (PMC);
Deng 2012, Feretzaki 2019, Cheng 1991, Santagostino 2025 and Tan 2026 were verified from their
abstracts/metadata through Europe PMC. No DOI, volume or page number below is from memory.

## 6. Where this must be woven

`submission/paper.tex` currently has only ED Figs. 1–4 in `\section*{Extended Data Figure Legends}`,
and it does not yet contain the divergence-vs-distance material that ED Fig. 6 carries (see
`response/ROADMAP_analysis_integration.md`, Batch 2 and Batch 4). These insertions travel with that
pending patch.

1. **Results — `\subsection*{Community structure of subtelomeric sharing}`.** One forward-pointing
   sentence at the end, after "PHR sharing is therefore constrained by chromosome-arm context rather
   than distributed indiscriminately among chromosome ends.": note that the sequence at the junction
   is itself a shared, copy-variable repeat family, and point at the new subsection. No numbers here.
2. **Results — new `\subsection*{The exchange junction is a copy-variable TERRA promoter}`,** inserted
   immediately before `\subsection*{Ectopic subtelomeric recombination within a pedigree}`. This is the
   only place in the current structure where the junction measurement is load-bearing: it sets up the
   pedigree exchange and supplies the reason the terminal bins of ED Fig. 6 are excluded.
3. **Discussion — paragraph beginning "This compartment is also internally organized."** Insert the
   first two Discussion paragraphs from §3 after "We see the same logic operating on other functional
   modules outside the NORs." This is where the similarity ↔ contact ↔ exchange loop is stated, and
   the promoter is the concrete shared repeat that loop can act on.
4. **Discussion — paragraph beginning "This concerted process operates at a characteristic scale."**
   Append the negative-result paragraph from §3, so the mechanism claim and the thing that is not the
   mechanism are read together.
5. **Methods — new `\subsection*{The terminal CpG island, GC content and annotatability}`,** after
   `\subsection*{Alignment, implicit pangenome graph and PHR detection}` and before
   `\subsection*{Pangenome graph, similarity and community detection}`. It shares the 2 kb
   both-sides-midpoint binning rule with the distance profile, so it must follow flank extraction.
6. **Methods — the discarding rule** (bins whose allelic divergence exceeds 0.4%, four times the
   interior baseline, are dropped) belongs in that same new subsection, because ED Fig. 6's hatching is
   generated from it. It must not be left only in a note.
7. **`\section*{Extended Data Figure Legends}` — ED Fig. 6 caption.** Insert the full caption from §4
   with `\label{fig:ed6}` so the Results cross-reference resolves. Confirm the number at insertion:
   the section currently ends at ED Fig. 4, so ED Fig. 6 implies an ED Fig. 5 in the same pass.
8. **Notes that must follow the manuscript, not contradict it:**
   `response/notes/telomere_distance_profile.md` and
   `response/notes/telomere_gene_vs_intergenic.md` — replace the "unmeasurable, but a biological
   finding" framing with the positive description (TERRA promoter array; ratio undefined). Also
   `response/notes/telomeric_cpg_island_terra.md` — replace the unverified 75% prevalence figure with
   39/46 \cite{Santagostino2025} and drop the unverified TAR1/CTCF sentence.
9. **Reviewer-response letter — `response/02_letter_draft_v1.md`,** tracked in
   `response/REAL_REVIEW_MAP.md` and `response/POINT_BY_POINT.md`. Add the point a referee is most
   likely to raise: *why is the first ~6 kb of every chromosome end excluded, and is that not an
   artefact?* The answer is now positive rather than defensive — the region is the TERRA promoter
   array, a copy-variable shared repeat, so the ratio is undefined there by construction, and the
   exclusion follows a quantitative rule rather than an eyeball. Cross-reference
   `response/06_answers_stats_and_bib.md` so the discarding rule sits with the other statistical
   controls, and record the change with its numbers in `response/03_change_log.md`.
10. **Figure provenance.** The hatching in `scripts/ci/telomere_end_ratio_curves.py` (and
    `scripts/ci/make_end_curve_figure.py`) must derive from the per-bin allelic values in the JSON and
    print a legend label matching the caption wording, so figure and caption cannot drift apart.

## 7. Open items and `AUTHOR-CHECK`s

- **`AUTHOR-CHECK` — Tan et al. and TERRA.** Whether Tan et al. discuss TERRA is unverified: their
  cached full text expired, and a fresh bioRxiv retrieval returned HTTP 429 from this host. Their
  abstract does not mention TERRA. Do not cite them for the promoter; cite them only for subtelomeric
  structural variation between haplotypes.
- **`AUTHOR-CHECK` — "TAR1 repeats at CHM13 telomeres are punctuated by CTCF-bound elements and
  CpG-island promoters"** (attributed in an earlier note to Dubocanin et al., bioRxiv 2022,
  `10.1101/2022.05.09.491186`). Could not be verified; **not used** above. Use the verified
  \cite{Deng2012,Stong2014} statements instead.
- **Prevalence figure corrected.** The "~75% of ends carry a CpG island within 15 kb" claim in
  `telomeric_cpg_island_terra.md` could not be traced to any source. Verified values: ~20/46 in the
  clone-based era (Nergadze 2009, Farnung 2010's "about half") and 39/46 in T2T-CHM13v2.0
  (Santagostino 2025). Use 39/46.
- **Our 88–95% and their 39/46 are not the same measurement.** Ours is a compositional signature
  (GC > 55% and CpG O/E > 0.6 in the first 2 kb); theirs is repeat-family assignment. The text
  therefore says "the CpG-island class of TERRA promoter", not "the 61-29-37 element at every end".
- **Not tested:** 29/37 tandem copy number per haplotype against that arm's paralog/allelic ratio —
  the direct test of the substrate claim, and feasible with the assemblies in hand.
- **Not tested:** whether the GC-rich tip on the ~7–12% of ends without a promoter is an unannotated
  promoter, a second promoter class, or a different GC-rich repeat.
- **Power limit:** the $p = 0.85$ test is arm-level (46 ends, 15 with tip genes). It excludes a large
  effect of gene placement on the arm ratio; it says nothing about copy-resolved homogenisation of
  DDX11L or WASH copies, which the multi-mapping assemblies cannot resolve.
- **"46 ends" is ambiguous in the project.** It is used for the annotation/TSS set here and in ED
  Fig. 3's caption ("26 of the 46 arm ends receive enough inter-arm alignment to be tiled"). Define
  the set once in Methods.
