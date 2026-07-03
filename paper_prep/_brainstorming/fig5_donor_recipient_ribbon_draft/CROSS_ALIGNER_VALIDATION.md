# Fig5 cross-aligner validation of the two autosomal candidates

Validation of the two autosomal PHR-exchange candidates shown in the Fig5
donor-recipient ribbons (built from the SweepGA/FastGA frequency-32 10:10 PAF).
Question: does the figure's f32 call reproduce under an independent aligner, and
at which coordinates?

Raw PAFs queried on octopus01 (moosefs), same PAN027-paternal vs PAN011-father
joint target:

- f32: `/moosefs/erikg/phrs/.wg-worktrees/agent-2727/.../pedigree_whole_genome_sweepga_fastga_frequency32/raw_paf/PAN027pat_vs_PAN011_joint.sweepga_frequency32_many_many_j0.paf.gz`
- wfmash v0.24.2-12-ge040aa10: `/moosefs/erikg/phrs/.wg-worktrees/agent-2636/.../pedigree_whole_genome_wfmash_p95_updated_bin/raw_paf/updated_bin_v0.24.2-12-ge040aa10/PAN027pat_vs_PAN011_joint.literal_p95.wfmash-v0.24.2-12-ge040aa10.paf.gz`

## chr9q/chr3q (panel C) — CONCORDANT across aligners

Recipient window PAN027#2#chr9 ~135.7-136.21 Mb. Donor chr3q, all forward strand:

| aligner | main donor block | identity |
|---|---|---|
| f32 (figure) | h1_chr3 202,417,382-202,549,099; h2_chr3 205,444,079-205,578,940 | 99.3-99.9% |
| f16 | identical to f32 (byte-for-byte) | 99.3-99.9% |
| wfmash -p95 | h1_chr3 202,461,289-202,535,329 | 99.82-99.93% |

The Fig5 ribbon's displayed donor (h2_chr3 ~205,538,813-205,558,808, the 10:10
filtered core) lies inside the raw f32 h2 block. h1 (~202.5 Mb) and h2 (~205.5 Mb)
are the two PAN011 chr3q haplotypes (same locus, ~3 Mb assembly offset). Three
aligner configurations agree at matching coordinates. Coordinates in
`../fig5_updated_binary_direct_alignment_review/validated_chr3_donor_coordinates.tsv`.

## chr5q/chr1p (panel B) — NOT reproduced by wfmash

Recipient window PAN027#2#chr5 ~181.6-182.15 Mb.

- f32/FastGA (figure): best inter-chromosomal match is chr1p, but it beats the
  same-chromosome (chr5) match by only a narrow margin (mean inter identity
  0.998 vs same-chromosome 0.995; mean delta ~0.003 on the main 28 kb block).
- wfmash -p95: the chr5q window maps to chr5 only (14 rows, all target chr5;
  zero chr1, zero other inter-chromosomal donor above the p95 threshold).

wfmash does report the chr9q->chr3q inter-chromosomal donor in the same run, so
its failure to report chr1 for chr5q is method-meaningful, not a blanket
inability to call subtelomeric inter-chromosomal homology. chr5q/chr1p is
therefore a weaker, aligner-dependent candidate.

## Manuscript consequence

`submission/paper.tex` (results, Methods, Fig5 caption) was updated to state that
only chr9q/chr3q is concordant across aligners and to mark chr5q/chr1p as the
weaker, aligner-dependent candidate. No coordinates in the figure were changed;
the f32 call is faithful, just not independently reproduced for chr5q/chr1p.
