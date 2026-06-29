# Fig5 Raw Many:Many IMPG 2kb Top-N Depth-Capped Scan

This supersedes the cancelled all-hit raw many:many IMPG scans.

Run settings:

- Existing raw unfiltered PAFs only; no new WFMASH/SweepGA alignment.
- 2 kb target windows from the previous full-genome BEDs.
- CHM13 centromere intervals removed before `impg similarity`; no flanking pad.
- `impg similarity` uses `--no-merge --num-mappings many:many --scaffold-jump 0 --max-depth 1`.
- Output stream keeps only interchromosomal top 20 rows per target window.
- Windows with more than 500 raw candidate rows are skipped and recorded under `outputs/skipped_windows/`.

Main outputs:

- `outputs/shards/`: bounded top-N shard outputs.
- `outputs/skipped_windows/`: per-shard over-dense skipped-window reports.
- `outputs/assembled/`: final compressed assembled top-N outputs after the dependency finalizer runs.
- `summaries/`: best-per-window plotting summaries after finalization.
