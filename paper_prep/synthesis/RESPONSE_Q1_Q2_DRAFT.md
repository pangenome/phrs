# Reviewer-response draft: D-PeerQ1 and D-PeerQ2

Response-letter material for the two deferred concerns logged in
`OPEN_REVIEWER_CONCERNS.md` as D-PeerQ1 and D-PeerQ2. This is not manuscript
prose.

## D-PeerQ1 - MAPQ0 multi-mapping

The Hi-C, Pore-C and CiFi maps were built with mapping-quality filters disabled:
Hi-C with HiC-Pro 3.0.0 and Bowtie2 2.5.3 under `MIN_MAPQ=0` and `RM_MULTI=0`,
Pore-C and CiFi with no minimum-MAPQ filter in pairtools 1.1.3. We removed that
filter because a subtelomeric PHR read that maps equally well to several
paralogous ends is the signal under study, and a standard MAPQ filter deletes it
outright. The pipeline does not, however, distribute multi-mapping mass across
the equally scoring loci. Each read pair is placed at one position only. A
multi-mapped pair contributes a single count to a single bin pair, not a
fractional count to every candidate bin pair. Two consequences follow. First, the
misplacement error is symmetric with respect to community membership: a read
placed at the wrong paralog is as likely to land on a between-community partner
as on a within-community partner, so placement adds noise but no directional
enrichment. Second, an inflated within-community signal would require placement
to favour within-community partners, which single-locus placement does not do.

The existing falsification is the unique-sequence flanking control. The 100 kb
immediately centromere-ward of every PHR boundary contains no inter-chromosomal
duplication, so multi-mapping is essentially absent there. In HG002 Hi-C the
between/within-community contact ratio (B/W, lower is stronger) is 0.027 at PHR
coordinates and 0.002 in the flanking unique-sequence windows at the same 50 kb
resolution, a 13.5-fold strengthening of within-community enrichment when the
paralogous sequence is removed from the window (10 kb: 0.029 versus 0.001, 29-fold;
100 kb: 0.021 versus 0.003, 7-fold). A multi-mapping artefact predicts the opposite
direction.

There is also an a fortiori argument. A PHR contains many near-identical copies,
but only a few copies can be physically engaged in any given contact. If only a
few are involved, single-locus placement can dilute a real signal, because a read
from an engaged copy may be placed at an unengaged paralog; it cannot manufacture
within-community enrichment, because manufacturing requires placement to prefer
within-community partners. Dilution and fabrication are not the same error, and
only the latter would rescue the reviewer's alternative.

The closure experiment is a MAPQ-strict re-run: re-derive contact matrices from
the upstream pairs behind the existing mcool files with MAPQ0 pairs dropped,
recompute B/W at PHR and flanking coordinates, and report the strict-versus-
random placement ratio per sample and region. Estimated effort is 1-2 days.

## D-PeerQ2 - RNA-mediated proximity (formerly 'RNA soup')

This stays in the Discussion as a hypothesis, and it is no longer speculation
without precedent. Liang et al., "Complementary Alu sequences mediate
enhancer-promoter selectivity", Nature 619 (2023),
doi:10.1038/s41586-023-06323-x, show that Alu elements in enhancers and promoters
pair through complementary Alu-derived RNA and thereby determine which enhancer
contacts which promoter. That is repeat-derived RNA driving sequence-specific
long-range spatial proximity in living cells, documented experimentally.

The argument in the Discussion has three steps. (a) The mechanism is documented
for Alu repeats in enhancer-promoter pairing. (b) Our observation is an aggregate
tendency of similar-sequence subtelomeric regions to sit closer in the nucleus,
which is the same class of phenomenon at a different scale: repeat-derived
sequence similarity converting into spatial proximity. (c) An aggregate signal is
compatible with only a small fraction of repeats being engaged at any moment.
Sparse participation dilutes an aggregate; it does not create one. "Only a few
copies participate" is therefore not an objection to an aggregate observation,
and we do not treat it as one.

The aggregate contact matrix cannot resolve which copies are engaged, so the
hypothesis is stated at exactly the level the data support: an aggregate
tendency for similar-sequence regions to sit closer, with a documented mechanism
that makes such an aggregate physically plausible.

We keep the paragraph because the alternative is worse. Omitting a documented
repeat-RNA proximity mechanism from a discussion of repeat-driven spatial
proximity would misrepresent the state of the field. The paragraph is labelled as
a hypothesis, it is placed in the Discussion rather than the Results, and it
carries the Liang et al. citation as its precedent.

## Open items

- The MAPQ-strict re-run has not been done; the ratio table is not yet produced.
- No ChIA-PET or Micro-C data are available for the same individual, so
  cross-validation by an independent proximity assay is not possible.
