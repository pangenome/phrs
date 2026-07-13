# Frozen CHM13 functional annotation term maps

## Outcome

This directory freezes a stable-identifier-only mapping from **61,312 physical CHM13v2.0 gene loci** to GO, HGNC gene groups, Reactome pathways, and the 26 explicit GFF3 biotypes. The audited PHR target contains 412 loci and is an exact subset of that universe. Every output annotation row uses the coordinate-anchored `copy_id` from `copy_universe.tsv.gz`; no gene symbol, alias, edit distance, case folding, or pseudogene-parent relation is used to assign a term.

## Identity model and joins

The physical key is `CHM13v2.0|seqid:start-end|strand|GFF_ID`. It remains unique for all loci, including 45S rDNA features that omit Liftoff `copy_num_ID`; the original value is retained separately. GFF3 `db_xref` and `Dbxref` are both parsed. Exact `GeneID` and `HGNC` joins are preferred. Exact OMIM, miRBase, and IMGT identifiers may bridge to one HGNC record through the frozen HGNC fields `omim_id`, `mirbase`, and `imgt`. Withdrawn HGNC IDs are redirected only through an exact, single-target HGNC merged report; withdrawn splits are ambiguous and rejected. A bridge is accepted only when it is one-to-one and concordant with every other supplied stable identifier. Symbols are display-only.

Pseudogenes are ordinary physical loci. A pseudogene gets only its own stable-ID annotations. The pipeline contains no parent-gene parser or inheritance operation. Loci without a usable stable identifier still receive their own biotype term and remain explicit in diagnostics.

Pre-existing enrichment spreadsheets, reports, and symbol-based draft notes under `_brainstorming/` are not inputs to this build. Biological claims from those drafts were neither copied nor used to fill mapping gaps.

Mapping status | Loci
--- | ---:
mapped_direct_geneid_only | 16,419
mapped_unique_hgnc | 43,992
unmapped_no_stable_xref | 876
unmapped_stable_xref_not_in_hgnc | 25

Accepted exact route | Loci
--- | ---:
direct_geneid | 16,419
exact_gene_to_hgnc | 1,224
exact_gene_to_hgnc|exact_hgnc_to_hgnc | 1
exact_gene_to_hgnc|exact_hgnc_to_hgnc|exact_omim_to_hgnc | 3
exact_hgnc_to_hgnc | 23,745
exact_imgt_to_hgnc | 588
exact_mirbase_to_hgnc | 2,045
exact_omim_to_hgnc | 16,299
exact_withdrawn_hgnc_redirect | 87
none | 901

## Coverage

Coverage counts physical loci, not unique symbols or genes. `HGNC_group` means at least one HGNC gene-group term; stable-ID coverage can therefore exceed functional-source coverage.

Source | All loci mapped | Coverage
--- | ---: | ---:
stable_id | 60,411/61,312 | 98.530%
biotype | 61,312/61,312 | 100.000%
GO | 20,363/61,312 | 33.212%
HGNC_group | 26,913/61,312 | 43.895%
Reactome | 11,365/61,312 | 18.536%

Full source-by-biotype coverage is in `outputs/coverage_by_biotype.tsv`. Required edge categories are summarized below.

Biotype | Loci | Stable ID | GO | HGNC group | Reactome
--- | ---: | ---: | ---: | ---: | ---:
C_region_pseudogene | 5 | 5/5 (100.000%) | 1/5 (20.000%) | 5/5 (100.000%) | 0/5 (0.000%)
J_segment_pseudogene | 6 | 6/6 (100.000%) | 0/6 (0.000%) | 6/6 (100.000%) | 0/6 (0.000%)
V_segment_pseudogene | 209 | 209/209 (100.000%) | 11/209 (5.263%) | 209/209 (100.000%) | 0/209 (0.000%)
ncRNA_pseudogene | 1 | 1/1 (100.000%) | 0/1 (0.000%) | 1/1 (100.000%) | 0/1 (0.000%)
pseudogene | 16,018 | 16014/16018 (99.975%) | 117/16018 (0.730%) | 1308/16018 (8.166%) | 0/16018 (0.000%)
transcribed_pseudogene | 1,262 | 1261/1262 (99.921%) | 185/1262 (14.659%) | 347/1262 (27.496%) | 9/1262 (0.713%)
C_region | 23 | 23/23 (100.000%) | 22/23 (95.652%) | 23/23 (100.000%) | 0/23 (0.000%)
J_segment | 79 | 79/79 (100.000%) | 14/79 (17.722%) | 79/79 (100.000%) | 0/79 (0.000%)
RNase_MRP_RNA | 1 | 1/1 (100.000%) | 0/1 (0.000%) | 1/1 (100.000%) | 0/1 (0.000%)
RNase_P_RNA | 1 | 1/1 (100.000%) | 1/1 (100.000%) | 1/1 (100.000%) | 0/1 (0.000%)
V_segment | 245 | 245/245 (100.000%) | 218/245 (88.980%) | 232/245 (94.694%) | 0/245 (0.000%)
Y_RNA | 7 | 7/7 (100.000%) | 0/7 (0.000%) | 7/7 (100.000%) | 0/7 (0.000%)
antisense_RNA | 19 | 19/19 (100.000%) | 3/19 (15.789%) | 19/19 (100.000%) | 0/19 (0.000%)
lncRNA | 18,389 | 18380/18389 (99.951%) | 188/18389 (1.022%) | 5737/18389 (31.198%) | 1/18389 (0.005%)
miRNA | 2,046 | 2045/2046 (99.951%) | 645/2046 (31.525%) | 2045/2046 (99.951%) | 0/2046 (0.000%)
misc_RNA | 37 | 37/37 (100.000%) | 7/37 (18.919%) | 6/37 (16.216%) | 1/37 (2.703%)
ncRNA | 49 | 49/49 (100.000%) | 4/49 (8.163%) | 27/49 (55.102%) | 0/49 (0.000%)
other | 13 | 11/13 (84.615%) | 2/13 (15.385%) | 10/13 (76.923%) | 2/13 (15.385%)
rRNA | 982 | 106/982 (10.794%) | 41/982 (4.175%) | 106/982 (10.794%) | 0/982 (0.000%)
scRNA | 4 | 4/4 (100.000%) | 3/4 (75.000%) | 3/4 (75.000%) | 0/4 (0.000%)
snRNA | 192 | 192/192 (100.000%) | 1/192 (0.521%) | 118/192 (61.458%) | 0/192 (0.000%)
snoRNA | 1,188 | 1187/1188 (99.916%) | 32/1188 (2.694%) | 532/1188 (44.781%) | 0/1188 (0.000%)
tRNA | 523 | 523/523 (100.000%) | 0/523 (0.000%) | 515/523 (98.470%) | 0/523 (0.000%)
telomerase_RNA | 1 | 1/1 (100.000%) | 0/1 (0.000%) | 1/1 (100.000%) | 0/1 (0.000%)
vault_RNA | 4 | 4/4 (100.000%) | 0/4 (0.000%) | 4/4 (100.000%) | 0/4 (0.000%)
protein_coding | 20,008 | 20001/20008 (99.965%) | 18868/20008 (94.302%) | 15571/20008 (77.824%) | 11352/20008 (56.737%)

## Term products

`copy_to_term.tsv.gz` has 532,209 unique `(copy_id, source, term_id)` rows. Duplicate source assertions are collapsed into sorted evidence sets; `support_count` counts distinct evidence codes, never join multiplicity. `term_metadata.tsv.gz` has 53,349 terms. `term_hierarchy.tsv.gz` has 74,794 GO/Reactome child-parent edges. GO metadata includes `biological_process`, `molecular_function`, or `cellular_component` namespace and obsolete status. Only direct GO annotations are assigned to copies; hierarchy edges are frozen separately so a downstream method must opt in explicitly to ancestor propagation.

GO assertions with a `NOT` qualifier are excluded. GO identifiers absent from the matching `go-basic` snapshot are also rejected: `none`.

## Unmapped-copy and inflation diagnostics

`copy_mapping_diagnostics.tsv.gz` contains exactly one row per physical locus, with stable mapping route and per-source term counts. `unmapped_copies.tsv.gz` emits one reason per missing functional source. The build fails on duplicate physical keys, a PHR copy outside the universe, ambiguous stable-ID indices, conflicting cross-references, HGNC group ID/name mismatches, missing term metadata, non-human Reactome rows, hierarchy endpoints absent from metadata, duplicate long-table keys, or output copy IDs outside the universe.

## Hand-audited duplicated families

Selection uses the exact, frozen display-name lists in `HAND_AUDIT`; it is not an annotation join. Up to four coordinate-distinct copies per listed name are written to `outputs/hand_audit_examples.tsv`. This deliberately includes pseudogenes and the unmapped rDNA control.

Family | Rows | Distinct copies | Stable-mapped | Notes
--- | ---: | ---: | ---: | ---
OR4F subtelomeric receptor family | 6 | 6 | 6 | coordinate-distinct loci retain independent copies; terms follow each locus's own ID
TUBB8/TUBB8B family with pseudogenes | 8 | 8 | 8 | coordinate-distinct loci retain independent copies; terms follow each locus's own ID
DUX4-like multicopy family | 12 | 12 | 12 | coordinate-distinct loci retain independent copies; terms follow each locus's own ID
45S rDNA arrays without GFF stable xrefs | 12 | 12 | 0 | 45S loci lack stable GFF xrefs and correctly receive biotype only

## Frozen provenance

The checked-in `sources/*.gz` objects are the analysis inputs; `sources/SOURCE_MANIFEST.tsv` records complete upstream-object and snapshot SHA-256 values, HTTP metadata, release/object identifiers, filters, retrieval time, and licenses. `fetch_sources.py --fetch` is an explicit renewal workflow and refuses to overwrite snapshots. Plain `fetch_sources.py` verifies the frozen checksums.

Source | Release/object | License | Upstream
--- | --- | --- | ---
ncbi_gene2go | NCBI Gene daily export, 2026-07-13 | [NCBI molecular data; NCBI places no restrictions on use or distribution](https://www.ncbi.nlm.nih.gov/home/about/policies/) | [download](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz)
hgnc_complete_set | HGNC complete set, upstream Last-Modified 2026-07-10 | [HGNC data release policy: no restrictions on access or use](https://www.genenames.org/download/) | [download](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt)
hgnc_withdrawn | HGNC withdrawn reports, upstream Last-Modified 2026-07-10 | [HGNC data release policy: no restrictions on access or use](https://www.genenames.org/download/) | [download](https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/withdrawn.txt)
go_basic | Gene Ontology 2026-06-15; http://purl.obolibrary.org/obo/go/releases/2026-06-15/go.owl | [CC BY 4.0](https://geneontology.org/docs/go-citation-policy/) | [download](https://current.geneontology.org/ontology/go-basic.obo)
reactome_ncbi | Reactome version 96 (2026-04-01; files refreshed 2026-06-19); doi:10.5281/zenodo.19581589 | [CC0 1.0 (Reactome annotation files)](https://reactome.org/license) | [download](https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt)
reactome_pathways | Reactome version 96 (2026-04-01; files refreshed 2026-06-19); doi:10.5281/zenodo.19581589 | [CC0 1.0 (Reactome annotation files)](https://reactome.org/license) | [download](https://reactome.org/download/current/ReactomePathways.txt)
reactome_relations | Reactome version 96 (2026-04-01; files refreshed 2026-06-19); doi:10.5281/zenodo.19581589 | [CC0 1.0 (Reactome annotation files)](https://reactome.org/license) | [download](https://reactome.org/download/current/ReactomePathwaysRelation.txt)

The local GFF inputs, source releases, repository object IDs, SHA-256 values, licensing notes, and derivation provenance are frozen in `outputs/INPUT_MANIFEST.tsv`. The universe GFF version is also encoded in the audited filename (`chm13v2.0_RefSeq_Liftoff_v5.2`).

## Rebuild and validation

From repository root:

```bash
python3 paper_prep/_brainstorming/chm13_copy_enrichment/fetch_sources.py
python3 paper_prep/_brainstorming/chm13_copy_enrichment/build_term_maps.py
python3 -m unittest discover -s paper_prep/_brainstorming/chm13_copy_enrichment/tests -v
```

The build is standard-library-only and deterministic: rerunning against unchanged snapshots reproduces byte-identical gzip outputs (gzip `mtime=0`) and the same `validation_summary.json`.
