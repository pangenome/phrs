# pafchop-rs

Streaming PAF row splitter for the Fig5 whole-genome direct-alignment checks.

Default behavior is intentionally small and safe for the sweepGA filter step:

- split every PAF row into `<= 10000 bp` query-axis fragments;
- linearly interpolate target coordinates plus match/block counts;
- drop inherited optional tags by default, avoiding huge duplicated CIGAR/cs tags;
- append concise `zp/zc/zl/zo/zs/ze/zts/zte` provenance tags;
- stream plain PAF on stdin/stdout so decompression/compression can be parallelized
  with `pigz`.

Build and test:

```bash
cargo test --manifest-path paper_prep/_brainstorming/pafchop-rs/Cargo.toml
cargo build --release --manifest-path paper_prep/_brainstorming/pafchop-rs/Cargo.toml
```

Run one file:

```bash
paper_prep/_brainstorming/pafchop-rs/scripts/chop_one.sh \
  raw.paf.gz chopped.10kb.paf.gz summary.tsv comparison_id 10000 8
```

Run a package in parallel, one process per comparison:

```bash
PAFCHOP_JOBS=3 PAFCHOP_THREADS_PER_JOB=8 \
  paper_prep/_brainstorming/pafchop-rs/scripts/chop_package_parallel.sh \
  paper_prep/_brainstorming/pedigree_whole_genome_sweepga_updated_bin 10000
```
