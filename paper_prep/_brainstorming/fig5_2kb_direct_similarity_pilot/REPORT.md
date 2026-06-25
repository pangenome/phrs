# Fig5 2 kb Direct-Similarity Pilot

## Purpose

This pilot evaluates a safer non-graph alternative for Fig5 whole-genome-style
support. It avoids `seqwish`/ODGI construction and instead summarizes existing
PAF evidence in 2 kb query-space bins around candidate and control regions.

Inputs:

- Candidate whole-genome WFMASH windows:
  `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/config/candidate_windows.tsv`
- Existing WFMASH 2 kb query-grid filtered support:
  `paper_prep/_brainstorming/pedigree_whole_genome_wfmash_p95_updated_bin/summaries/query_grid_filter_candidate_window_support.tsv`
- Existing joint-parent SweepGA direct PAFs:
  `paper_prep/_brainstorming/pedigree_direct_sweepga_joint_parent/filtered_paf/*.paf.gz`
- Pilot window manifest:
  `config/evaluation_windows.tsv`

No new all-vs-all alignment or graph construction was run for this pilot. The
script streams already-computed compressed PAFs and slices the small candidate
and control windows.

## Outputs

- `summaries/bin_target_support.tsv`: per-method, per-event, per-2 kb bin
  support by target chromosome/arm, with weighted identity and match-distance
  where PAF fields allow it.
- `summaries/window_target_summary.tsv`: per-method, per-event aggregate support
  by target chromosome and class.
- `summaries/method_comparison.tsv`: one row per event/method with top target,
  expected-target support, same-chromosome support, and control-class support.
- `scripts/run_direct_similarity_pilot.py`: reproducible PAF-only summarizer.
- `run.sh`: rerun wrapper.

## Candidate Result

The direct 2 kb summaries recover only weak chr3q evidence relative to dominant
same-chromosome chr9 support.

| event | method | chr3 support | same-chrom chr9 support | interpretation |
| --- | --- | ---: | ---: | --- |
| PAN027 chr9q candidate | WFMASH 2 kb 1:1 | 2,000 bp in 1 bin | 435,398 bp in 223 bins | chr3 is detected, but as a single conservative 2 kb bin. |
| PAN027 chr9q candidate | SweepGA many/four-many | 0 bp | 619,724 bp in 250 bins | Direct subtelomeric PAF support does not recover chr3 for PAN027. |
| PAN027 chr9q candidate | SweepGA one-one/simple | 0 bp | 499,208 bp in 250 bins | Strict filters retain only same-chromosome support. |
| PAN028 chr9q candidate | WFMASH 2 kb 1:1 | 10,000 bp in 5 bins | 385,875 bp in 197 bins | chr3 is detected in the conservative WFMASH grid. |
| PAN028 chr9q candidate | SweepGA many/four-many | 7,765 bp in 5 bins | 784,316 bp in 250 bins | Direct subtelomeric PAF also recovers chr3 for PAN028, but same-chromosome support dominates. |
| PAN028 chr9q candidate | SweepGA one-one/simple | 0 bp | 681,106 bp in 250 bins for one-one; 0 bp for simple threshold | Strict SweepGA filters do not retain chr3. |

Conclusion: the 2 kb direct-query approach is useful as a sanity check, but it
does not replace the graph/untangle story as-is. PAN028 has reproducible chr3
support in both WFMASH and many/four-many SweepGA direct summaries. PAN027 has
only a very small WFMASH 2 kb chr3 signal and no chr3 support in the direct
SweepGA subtelomeric PAF filters evaluated here.

## Control Result

PAR and acrocentric controls behave as expected for the available target content.

- PAN027 PAR1 positive control (`PAN027pat_vs_PAN011_joint`) recovers chrY
  support in SweepGA many/four/one-one filters: 144,103 bp across 73 bins,
  alongside stronger same-chromosome chrX support.
- PAN028 PAR control (`PAN028mat_vs_PAN027_joint`) is not an X-to-Y test because
  the PAN027 target set has no chrY sequence. It correctly recovers chrX support
  instead: 868,720 bp across 232 bins in the many-many filter.
- PAN027 acrocentric p-arm controls are clean same-chromosome positives:
  chr13p, chr14p, chr15p, chr21p, and chr22p each recover approximately the full
  500 kb expected chromosome in many/four/one-one/simple filters.
- PAN028 acrocentric p-arm controls recover expected acrocentric signal but
  often top-rank cross-acrocentric targets in many/four-many filters, consistent
  with shared acrocentric p-arm sequence rather than a chr3-specific artifact.

## Caveats

- WFMASH 2 kb filtered rows are available here only for the two Fig5 chr9q
  candidates, not for the PAR/acros controls.
- SweepGA direct PAFs are subtelomeric 500 kb window PAFs with native coordinates
  embedded in sequence names. The pilot summarizes them in local 0-500 kb query
  coordinates to match the direct-query-bin question.
- Aggregated aligned bp can exceed the 500 kb window length under many:many or
  four:many filters because multiple target mappings can overlap the same query
  bin. Use `bins_with_support` and target class alongside bp totals.
- Match-distance is computed as `1 - matches/alignment_block_length` from PAF
  columns where available. It is a direct PAF proxy, not a full re-alignment of
  extracted 2 kb FASTA windows.

## Reproduction

From the repository root:

```bash
bash paper_prep/_brainstorming/fig5_2kb_direct_similarity_pilot/run.sh
```
