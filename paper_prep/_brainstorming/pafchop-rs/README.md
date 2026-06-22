# pafchop-rs

Streaming PAF row splitter for the Fig5 whole-genome direct-alignment checks.

Default behavior is intentionally strict and safe for identity-sensitive
sweepGA filter steps:

- split every PAF row into `<= 10000 bp` query-axis fragments;
- require `cg:Z` CIGAR input so target coordinates plus match/block counts can
  be recomputed from alignment operations instead of interpolated;
- recompute clipped `cg:Z`, clipped `cs:Z` when present, and identity-relevant
  tags (`NM:i`, `dv:f`, `de:f`, `df:i`) when those tags are present upstream;
- drop inherited optional tags by default, and never copy stale
  alignment-derived tags even with `--keep-tags`;
- append concise `zp/zc/zl/zo/zs/ze/zts/zte` provenance tags;
- stream plain PAF on stdin/stdout so decompression/compression can be parallelized
  with `pigz`.

See `PAF_SEMANTICS_VALIDATION.md` for the PAF column/tag contract and the
classification of the older f16 chopped outputs.

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
