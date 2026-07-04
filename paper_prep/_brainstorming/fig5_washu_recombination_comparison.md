# Fig. 5 WashU Recombination Comparison Note

**Status checked:** 2026-07-04.

This note records the provenance boundary for the Fig. 5 comparison between our
whole-genome PHR-exchange scan and the WashU pedigree recombination analysis.
It is intended for manuscript, caption, and downstream validation-planning
passes.

## Current Citation Status

The WashU pedigree assembly paper can now be cited directly as the bioRxiv
preprint `Cechova2025`:

- Monika Cechova et al., "Complete genomes of a multi-generational pedigree to
  expand studies of genetic and epigenetic inheritance", bioRxiv v1,
  DOI `10.64898/2025.12.14.693655`.
- The preprint is indexed by PubMed as PMID `41473289` and PMCID
  `PMC12746033`.

This supersedes older internal notes that said the WashU pedigree manuscript was
not yet on bioRxiv. Those notes reflected earlier author communication before
the preprint became public.

Use the preprint citation for:

- the existence and provenance of the PAN010/PAN011/PAN027/PAN028 T2T pedigree
  assemblies;
- general statements that this pedigree is a public telomere-to-telomere
  multigeneration benchmark;
- methods that consume those assemblies as source data.

Do not use the preprint alone for the Fig. 5 lineage-specific negative
recombination statement unless a future manuscript version or supplement exposes
the corresponding recombination table/call.

## Fig. 5 Comparison Result

Our Fig. 5 manuscript analysis compares the PAN027 paternal haplotype against a
joint two-haplotype PAN011 father target. The whole-genome 2 kb
IMPG-similarity scan highlights a chr9q/chr3q candidate in which non-homologous
chr3q donor windows outscore the same-chromosome or homologous father match
across adjacent windows. The same chr3q donor locus is reproduced by an
independent wfmash alignment at greater than 99.8% identity.

Manual comparison with the WashU pedigree recombination annotation, as
communicated by M. Cechova and T. Marschall, agrees that the relevant PAN011 to
PAN027 paternal lineage does not contain an ordinary chr3 crossover or standard
homologous chr3 recombination event that would explain the chr9q/chr3q signal.

Therefore, the current manuscript statement should remain:

- cite `Cechova2025` for the pedigree assemblies and preprint;
- cite `M. Cechova and T. Marschall, personal communication` for the specific
  negative comparison to their recombination annotation;
- describe the chr9q/chr3q result as a putative PHR-exchange-compatible
  candidate, not a breakpoint-resolved event call.

## Claim Boundary

The comparison rules out a simple explanation in which our chr9q/chr3q signal is
just an ordinary homologous chr3 crossover already called in the WashU pedigree
recombination map. It does not prove that the candidate is a clean
interchromosomal crossover. The Fig. 5 wording should keep the current caution:
the autosomal signal is recombination-compatible and concordant across aligners,
but remains putative until tested by descendant transmission, orthogonal
cytogenetics such as FISH, or long-read single-gamete recombination data.

## Update Trigger

If the WashU authors post a revised preprint, accepted manuscript, or supplement
that includes the relevant recombination annotation, update this note and the
Fig. 5 manuscript text as follows:

1. Replace the personal-communication parenthetical with a direct citation to
   the table, figure, or supplement containing the chr3-negative call.
2. Keep `Cechova2025` for assembly provenance unless a journal version receives
   a new DOI.
3. Re-check whether the manuscript's phrase "ordinary chr3 crossover" should be
   broadened to the exact term used by the WashU annotation, for example
   "crossover", "recombination", "haplotype switch", or "switch error".
