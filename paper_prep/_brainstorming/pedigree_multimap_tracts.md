# WashU Untangle Multimap-Aware Tract Calls

Script: `scripts/pedigree/untangle_multimap_tracts.py`

Primary outputs:

- `scripts/pedigree/untangle_multimap_tracts.tsv`
- `scripts/pedigree/untangle_multimap_tract_summary.tsv`
- `paper_prep/_brainstorming/pedigree_multimap_tracts/representative_segments.tsv`
- `paper_prep/_brainstorming/pedigree_multimap_tracts/representative_segments.svg`
- `paper_prep/_brainstorming/pedigree_multimap_tracts/plot_representative_segments.py`

## What Was Tried

The previous lower-threshold check merged only `nth.best == 1` intervals by a
single donor haplotype. This pass instead groups every exact child/query
interval in the WashU `odgi untangle` BEDs, keeps top-N hits that are within a
score delta of the best hit, and treats those hits as the interval's donor
equivalence class. Adjacent intervals are merged only when those equivalence
classes remain compatible under the selected bridge mode.

The default sensitivity runs are:

| setting | top-N | score delta | min score | max bridge gap | bridge mode |
|---|---:|---:|---:|---:|---|
| m1000_default | 4 | 0.001 | 0.8 | 0 | arm |
| m1000_bridge1kb | 8 | 0.002 | 0.8 | 1000 | community |

The committed run processes the available m1000 BEDs directly. The script also
accepts lower-threshold BEDs through the `merge_m` field in `--run` specs, for
example `--run lower_m0:0:8:0.002:0.8:1000:community`, but those large m0 BED
intermediates are not committed in this worktree. For that reason the summary
table carries forward the existing `m0/n1` first-best lower-merge denominator
from `scripts/pedigree/patch_tract_lower_merge_summary.tsv` as a comparison
rather than silently treating it as multimap-aware evidence.

## sweepga Check

`sweepga --help` shows that sweepga accepts FASTA inputs or a single PAF and
then applies plane-sweep filtering with `--num-mappings`. The WashU evidence
available here is already the downstream `odgi untangle` BED output, with
child/query coordinates, reference path coordinates, scores, and `nth.best`.
There is no PAF-equivalent retained in the repository or the WashU untangle
directory. Re-running sweepga would therefore require reconstructing a separate
FASTA/PAF alignment path rather than preserving odgi's graph-specific interval
calls. For this task, the appropriate equivalent is an interval sweep over the
untangle BED rows themselves.

Relevant sweepga interface excerpt:

```text
Input files: FASTA (1+) or PAF (1 only), auto-detected
            1 FASTA: align to self and filter
            2+ FASTA: align all pairs and filter
            1 PAF: filter alignments
          Aligner to use for FASTA input
  -n, --num-mappings <NUM_MAPPINGS>
          Output PAF format instead of default .1aln (text instead of binary)
```

## Results

Under the recommended default `m1000_default`, the multimap-aware caller
returns 373 interchromosomal candidate tracts. The matched
first-best projection contains 585 tracts, so
the multimap-aware representation merges 212 and
splits 0 relative to the arbitrary first-best
view. Ambiguity is explicit: 201 tracts are not
unique donor haplotypes, including 22
unique-arm/multiple-haplotype tracts, 165
same-community ambiguous tracts, and
14 cross-community ambiguous tracts.

With the more permissive `m1000_bridge1kb` sensitivity, the caller returns
141 interchromosomal candidate tracts, merges
444 relative to first-best, and records
0 tracts with non-zero bridged sequence. This confirms
that consecutive m1000 intervals can be merged through compatible donor
equivalence classes. In this m1000 run the compatible joins are mostly zero-gap
adjacency rather than non-zero missing sequence; some intervals remain genuinely
unresolved or cross-community ambiguous.

Length distributions against the primate literature windows are summarized in
`scripts/pedigree/untangle_multimap_tract_summary.tsv`. For the recommended
default, median tract length is 3567 bp (IQR
16429 bp; min 493 bp; max 500000 bp).
Counts in the descriptive windows are:

- 22-95 bp: 0/373
  (0.000000)
- 318-688 bp: 1/373
  (0.002681)
- 159-1376 bp: 125/373
  (0.335121)

## Interpretation

The analysis no longer treats `nth.best == 1` as uniquely true when tied or
near-tied donor paths exist. It strengthens the tract-length audit by showing
which m1000 intervals are mergeable under explicit donor equivalence classes,
and by identifying same-arm, same-community, and cross-community ambiguity.

The result remains a supportive compatibility analysis. It does not by itself
prove conversion or crossover mechanisms, because repeated subtelomeric
haplotypes often make the exact donor unresolved. The useful claim is cautious:
some WashU candidate regions are compatible with merged tract interpretations
once multimapping is represented, while a substantial fraction should remain
ambiguous rather than be forced into a first-best donor.
