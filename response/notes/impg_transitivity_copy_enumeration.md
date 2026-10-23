# IMPG transitivity copy-enumeration of core sequences

Status: IN PROGRESS
Author: sxgc.deepseek-flash-writer (child agent)
Date: (see changelog)

## 0. Changelog
- (t0) skeleton created

## 1. Question
Placeholder.

## 2. Method
Placeholder.

## 3. Cores
Placeholder.

## 4. Results
Placeholder.

## 5. Timing / cost
Placeholder.

## 6. Cross-check vs FastGA
Placeholder.

## 7. What this establishes / does not
Placeholder.

## 8. Failures and anomalies
Placeholder.

## 0a. Changelog (append-only)
- t0: skeleton created.
- t1: resources probed; **the provided 37 GB IMPG index is not queryable** (see §8 failure F1). Pivoting to a self-built index over self-generated FastGA PAFs. Decision recorded, no numbers estimated.

## 8. Failures and anomalies (running log)
### F1 - `hprc25272.bidirectional.impg` cannot resolve CIGAR data (orphaned index)
Evidence:
- `impg query -i /moosefs/.../hprc25272.bidirectional.impg -a <placeholder> -r CHM13#0#chr10:38000-53000 -d 1k -o bed`
  -> `thread 'main' panicked at src/impg.rs:1217: index out of bounds: the len is 1 but the index is 3048`
  (the index references >=3049 alignment files; `-i` still requires an `--alignment-*` argument, and
  `alignment_file_index` is an index into that list.)
- With the shipped list `all-vs-all.1Mb.p95.paf.list` (18827 entries, all under `/scratch/...`, `/scratch` copies deleted):
  `Failed to open compressed file '/scratch/all-vs-all.1Mb.p95/hprcv2.1Mb...paf.gz'` (impg.rs:623).
- With a rewritten list pointing at the surviving copies in
  `/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95/` (18827 `*.paf.gz` present) and at
  `all-vs-all.p95/` (18226 `*.paf.gz`), the query dies in BGZF seek:
  `slice index starts at 22159 but ends at 0` (noodles-bgzf data.rs:45).
- Timestamps explain it: index mtime 2026-01-19T22:28:06Z; `all-vs-all.p95/*.paf.gz` 2026-01-21;
  `all-vs-all.1Mb.p95/*.paf.gz` 2026-03-20. Both surviving PAF sets are **newer** than the index, so the
  BGZF virtual offsets stored in the index (data_offset per PAF record) do not address these files.
Conclusion: the pre-built index is orphaned from its source alignments; every CIGAR-dependent query path
(bed/bedpe/paf, and any identity filtering) is unusable. Reported, not worked around silently.
