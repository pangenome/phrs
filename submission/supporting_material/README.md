# Supporting material

`Fig5_pedigree_runs.tsv` lists the 82 drawn non-homologous runs of the three WashU transmissions
(father to daughter, mother to daughter, daughter to granddaughter), scanned on the assemblies we
corrected with each person's own reads. It is written by
`scripts/pedigree/fig5_v131_check/pedigree_run_stats.py` from the runs of `fig5_runs.sh`.
Columns: transmission, run, category, query and donor intervals (0-based, end excluded, in the
corrected assemblies), length, query and donor arm and Leiden community, whether both arms share a
community, mean same-chromosome and non-homologous identity, and Merqury error k-mers observed and
expected in the query and donor intervals.

## Extended Data Figure 2 copy-number ontology tables

`ED_Fig2_supported_ontology_terms.tsv` contains exactly the 209 primary
midpoint-supported exact V7 `(collection, relation, term_id)` rows. It is a
field-preserving subset of
`paper_prep/_brainstorming/chm13_copy_enrichment/ontology_v7/TERM_RESULTS.tsv.gz`,
ordered by the frozen integer `hypothesis_index`. The table retains the exact
physical-copy 2-by-2 counts, eligible-universe sizes, copy fractions, fold and
odds-ratio summaries, exact upper-tail P value, within-collection BH q value,
global Holm value, and support/background status. Each count is a
coordinate-distinct CHM13 GFF3 physical copy; records are never collapsed by
gene name, functional source, family, locus family, duplicated segment,
community, or ontology provenance.

`ED_Fig2_arm_functional_components.tsv` is the complete six-class by 48-end
physical-copy matrix plotted in Extended Data Fig. 2 (288 data rows). Rows are
ordered by display class and then by the figure's `community_order_position`:
C1 through C15, similarity-tree order within each community, followed by the
seven no-signal arms. These six classes are post-inference display summaries of
the exact supported terms. They add no hypotheses, P values, or support calls.

Both tables are rebuilt deterministically from retained V7 and community
attribution artifacts by running:

```bash
python3 submission/scripts/make_ed2_supporting_material.py
```

Stable repository paths:

- <https://github.com/pangenome/phrs/blob/main/submission/supporting_material/ED_Fig2_supported_ontology_terms.tsv>
- <https://github.com/pangenome/phrs/blob/main/submission/supporting_material/ED_Fig2_arm_functional_components.tsv>
- <https://github.com/pangenome/phrs/blob/main/submission/scripts/make_ed2_supporting_material.py>
