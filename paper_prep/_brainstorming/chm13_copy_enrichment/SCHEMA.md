# CHM13 term-map schemas

All tables are tab-separated UTF-8. Compressed tables use deterministic gzip
headers (`mtime=0`). Coordinates are GFF3 1-based, end-inclusive coordinates.

- `outputs/copy_universe.tsv.gz`: one row per CHM13 physical `gene` feature.
  `copy_id` is the primary key. `gene_name_display_only` is never a join key.
- `outputs/copy_to_term.tsv.gz`: unique long-table key
  (`copy_id`, `source`, `term_id`). `support_count` counts distinct evidence
  codes after exact duplicate collapse; it is not a copy count.
- `outputs/term_metadata.tsv.gz`: primary key (`source`, `term_id`). GO
  `namespace` contains the ontology aspect.
- `outputs/term_hierarchy.tsv.gz`: unique child-parent ontology/pathway edges.
  GO rows include `is_a` and the relationships present in `go-basic`.
- `outputs/copy_mapping_diagnostics.tsv.gz`: exactly one row per `copy_id`.
- `outputs/unmapped_copies.tsv.gz`: one row per copy and functional source for
  which no term was assigned. Biotype is intentionally absent because every
  locus has exactly one biotype term.
- `outputs/coverage_by_biotype.tsv`: physical-locus denominators for every
  source and all 23 GFF biotypes, plus `ALL` totals.
- `outputs/hand_audit_examples.tsv`: coordinate-level examples selected by an
  exact display-name allowlist solely for human review.
- `outputs/INPUT_MANIFEST.tsv`: local audited inputs and SHA-256 values.
- `outputs/validation_summary.json`: machine-readable invariants and counts.
- `sources/SOURCE_MANIFEST.tsv`: upstream object provenance and frozen
  snapshot checksums.

`validate_outputs.py` reads the frozen outputs independently of the build
loop, compares every copy to the audited GFF, recomputes coverage, and rejects
duplicate keys, missing metadata, out-of-universe rows, or symbol routes.

Term sources are `biotype`, `GO`, `HGNC_group`, and `Reactome`. Biotype IDs use
the `BIOTYPE:` prefix, HGNC groups use `HGNC_GROUP:`, and GO/Reactome retain
their native stable IDs.
