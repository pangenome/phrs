# Fig5 Raw 2 kb IMPG Similarity, Plane-Sweep 4:4

This is a bounded rerun of the full-genome 2 kb IMPG similarity scan using the
existing raw WFMASH and SweepGA/FastGA PAFs. It applies IMPG/SweepGA plane-sweep
mapping filtering with `--num-mappings 4:4` before local similarity.

Run settings:

- Existing raw unfiltered PAFs only; no new WFMASH/SweepGA alignment.
- 2 kb target windows from the previous full-genome BEDs.
- CHM13 centromere intervals removed before `impg similarity`; no flanking pad.
- `impg similarity` uses `--no-merge --num-mappings 4:4 --scaffold-jump 0`.
- It is non-transitive: the command does not pass `--transitive`.
- Output stream keeps only interchromosomal top 20 rows per target window.
- Windows with more than 500 raw candidate rows are skipped and recorded under `outputs/skipped_windows/`.
