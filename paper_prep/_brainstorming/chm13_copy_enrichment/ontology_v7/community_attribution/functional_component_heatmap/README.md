# Functional component heatmap

Review artifact for visualizing the V7 CHM13 copy-number-aware ontology results
by chromosome end.

The plotted unit is coordinate-distinct PHR copy burden. The six rows are
post-inference display classes over exact supported GO/Reactome term rows; they
are not additional enrichment tests and they do not treat redundant ontology
terms as independent biology.

Run from the repository root:

```sh
Rscript paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v7/community_attribution/functional_component_heatmap/make_functional_component_heatmap.R
```

Primary review outputs:

- `functional_component_arm_heatmap.community_order.pdf/png` - chromosome ends
  grouped by C1-C15 community, with tree order within each community.
- `functional_component_arm_heatmap.nj_order.pdf/png` - chromosome ends ordered
  by the NJ tree tip layout, with no-signal arms appended.
- `functional_component_arm_heatmap.review_pages.pdf` - both orderings in one
  two-page PDF.
- `functional_component_arm_matrix.tsv` - the exact six-class arm copy matrix.

Audit outputs:

- `exact_supported_term_arm_matrix.tsv.gz` - all exact primary-supported V7
  term rows by chromosome end.
- `exact_supported_term_copy_pattern_summary.tsv` - groups exact terms that
  share identical chromosome-end copy patterns, useful for deciding whether a
  larger exact-term heatmap is too redundant for a figure.
