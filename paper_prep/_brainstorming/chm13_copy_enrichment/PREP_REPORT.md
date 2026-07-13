# CHM13 PHR physical gene-copy input audit

## Decision

The analysis universe is the **61,312 physical `gene` feature rows** in
`data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz`. Each row is retained once and is
identified by its GFF3 `ID`; repeated `gene_name` values are deliberately not
collapsed. The primary PHR membership rule is gene midpoint in one of the 37
intervals in `data/chm13.phrs.bed`. Any positive-base overlap is retained as a
sensitivity assignment.

The rebuilt result agrees with the draft audit totals: **61,312 universe loci,
402 midpoint PHR loci, and 412 any-overlap PHR loci**. There is therefore no
departure to reconcile. The 10 additional overlap-only rows are real GFF3
features spanning a PHR boundary, listed below; they are not inferred or
propagated copies.

This audit is CHM13-only. It does not read HPRCv2 assemblies or projected
annotations, and it does not change the manuscript.

## Audited inputs and provenance

All analysis inputs are repository files. `analysis_ready/PROVENANCE.tsv`
records the same byte sizes and SHA-256 digests in machine-readable form.

| Role | Repository path | Bytes | SHA-256 |
|---|---|---:|---|
| 37 target PHRs | `data/chm13.phrs.bed` | 3,612 | `03cc73f049e9625d131137d8ab7fbc5f52833c2aade52b9b6635d5a874b55cb9` |
| physical gene annotation | `data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz` | 56,368,414 | `a1c8e61cb4e60a3af3a18599b7d5551a72a1b0317bdffad42ae7fa36e73da968` |
| p/q boundaries and chromosome extents | `data/chm13v2.0_cytobands_allchrs.bed` | 30,753 | `09c7fa3ca3e222ebaabbecedc3086aa2543999824c0966f96b0c3f31c1126ea2` |

The PHR BED and GFF3 entered the repository in git commit
`5860210217aa01e1226a881d343aedf9fe85e5c1` (2026-05-31); the cytoband source
entered in `882c6cafc434332327078a13da6422b2ebd5cad8` (2026-07-07). Their git blob
IDs are, respectively, `f7c8769ff12e756f267e413812124f7598d1baa3`,
`cb3ca8cd582b5e13dd89b8cbc17bb34893bd2cfe`, and
`bec810221ab352144229d07d4bc80e00f7142814`.

The cytoband BED is an explicit auxiliary coordinate input, not an annotation
source. It is used only to define chromosome extents and the first q-arm base.
The gene universe comes exclusively from the specified GFF3.

## Coordinate and contig audit

- `chm13.phrs.bed` is interpreted as BED: zero-based, half-open
  `[start0, end0)`. It has exactly 37 nonempty, pairwise-disjoint intervals on
  37 distinct arms (18 p, 19 q). Their union is **6,014,981 bp**. No endpoint
  correction or manual interval edit is applied.
- GFF3 genes are interpreted as one-based, closed `[start1, end1]` and
  converted to `[start1 - 1, end1)` before intersection. Both coordinate forms
  are emitted in the locus tables. Feature length is `end0 - start0`.
- The primary zero-based midpoint is `floor((start0 + end0) / 2)`. For an
  even-length feature this is the right central base. The alternative left
  central base gives the same real-data total (402) and changes no real-data
  assignment, but the chosen convention is fixed and unit-tested.
- A midpoint hit requires `phr_start0 <= midpoint0 < phr_end0`. Any-overlap
  sensitivity requires `gene_start0 < phr_end0` and
  `gene_end0 > phr_start0`, so endpoint-only contact is not overlap.
- The GFF3 contains exactly the canonical CHM13 chromosome names `chr1` through
  `chr22`, `chrX`, and `chrY`. The PHR BED contains the canonical subset
  `chr1` through `chr22` plus `chrX`; every target contig exists in the GFF3.
  The 693 chrY genes remain in the genome-wide universe even though chrY has
  no target PHR in this BED.
- Cytobands are required to be contiguous from zero to the chromosome end,
  with all p bands before all q bands. Each gene is assigned to p or q by its
  midpoint relative to the first q-band base. The same rule assigns each
  terminal PHR. `chm13_arm_summary.tsv` exposes all 48 boundaries and the
  resulting universe/target counts for review.

## Physical-copy universe audit

The builder streams the GFF3 and emits one locus for every row whose third
column is exactly `gene`. It requires `ID`, `gene_name`, and `gene_biotype` but
does not filter on biotype, strand, coding status, copy-number metadata, or PHR
membership.

The source has 61,312 gene rows and 61,312 unique gene IDs. It has only 58,231
unique gene names, or 3,081 repeat-name rows beyond a symbol-deduplicated
universe. Those repeat-name rows are the physical copy-number information and
are retained. There are also 12 rows beyond the 61,300 unique
chromosome/start/end/strand tuples; these coincident but separately identified
GFF3 gene rows are likewise retained. A coordinate coincidence is not grounds
for synthesis or deduplication.

All observed gene biotypes remain present, including 20,008 protein-coding,
18,389 lncRNA, 16,018 pseudogene, 2,046 miRNA, and 1,262 transcribed-pseudogene
rows, plus the less abundant ncRNA and immune-segment classes. The midpoint
target contains 204 pseudogenes, 104 lncRNAs, 51 miRNAs, 24 protein-coding
genes, 18 transcribed pseudogenes, and one misc_RNA. Thus the target totals are
not coding-gene-only counts.

Every universe output row carries:

- its unique source `locus_id` and original `gff_line`;
- its source chromosome, converted coordinates, strand, biotype, and Liftoff
  copy metadata when present;
- `record_origin=gff3_gene_row`; and
- independently computed midpoint and overlap flags.

The integration test compares the ordered tuple `(gff_line, ID, gene_name,
chromosome, start0, end0, strand, gene_biotype)` for all 61,312 output rows
against the raw GFF3. This proves a one-to-one source mapping: no output row is
synthetic, no symbol is propagated to another chromosome, and no physical GFF3
gene row is lost.

## Midpoint versus overlap sensitivity

The 402 midpoint loci are all among the 412 any-overlap loci. Target PHRs do
not overlap one another, so each assigned gene has at most one primary PHR and
one overlap PHR. The 10 overlap-only boundary-spanning features are:

| Locus ID | Chromosome arm | Gene interval `[start0,end0)` | PHR overlap bp | Biotype |
|---|---|---:|---:|---|
| FAM41C | chr1_p | 297044–306720 | 662 | lncRNA |
| LINC01237 | chr2_q | 242383552–242589193 | 15,066 | lncRNA |
| ZNF595 | chr4_p | 52752–87637 | 15,335 | protein_coding |
| LOC101929756 | chr7_p | 177343–183674 | 1,032 | lncRNA |
| RPL23AP53 | chr8_p | 5120–29123 | 6,761 | transcribed_pseudogene |
| IQSEC3 | chr12_p | 58887–174477 | 13,896 | protein_coding |
| LOC101929650 | chr17_q | 84166185–84240296 | 11,192 | lncRNA |
| ROCK1P1 | chr18_p | 262854–276133 | 4,163 | transcribed_pseudogene |
| RPL23AP82 | chr22_q | 51270583–51313167 | 15,970 | transcribed_pseudogene |
| SPRY3 | chrX_q | 153837808–154007734 | 80,816 | protein_coding |

The 10-row difference is therefore explained entirely by the predeclared
membership rules: each feature intersects a PHR, but its midpoint lies outside
that PHR. The midpoint table is the primary analysis input; the overlap table
is a sensitivity input and must not be silently substituted.

## Analysis-ready artifacts

`analysis_ready/` is generated, committed, and treated as immutable. Do not
edit these tables manually.

| File | Data rows | Purpose |
|---|---:|---|
| `chm13_phr_intervals.tsv` | 37 | normalized PHR targets, arm calls, source BED lines |
| `chm13_gene_loci.tsv.gz` | 61,312 | complete physical-copy universe |
| `chm13_phr_gene_midpoint.tsv` | 402 | primary target loci |
| `chm13_phr_gene_any_overlap.tsv` | 412 | overlap sensitivity target loci |
| `chm13_arm_summary.tsv` | 48 | p/q boundaries and per-arm universe/target totals |
| `PROVENANCE.tsv` | 4 | input hashes, byte sizes, roles, builder version |
| `MANIFEST.sha256` | 6 | checksums, byte sizes, and row counts for every table above |

`MANIFEST.sha256` is the authoritative output digest list. The compressed
universe is written with an empty gzip filename and timestamp zero, making its
bytes deterministic under the audited builder.

## Rebuild and validation

Run from the repository root:

```bash
python3 paper_prep/_brainstorming/chm13_copy_enrichment/build_inputs.py
python3 paper_prep/_brainstorming/chm13_copy_enrichment/build_inputs.py --check
python3 -m unittest discover \
  -s paper_prep/_brainstorming/chm13_copy_enrichment \
  -p 'test_*.py' -v
```

The first command rebuilds every table from raw repository inputs into a
temporary directory and replaces `analysis_ready/`. The second
rebuilds into a fresh temporary directory and byte-compares every file with
the committed artifacts. The test suite covers:

- BED/GFF coordinate conversion and both PHR boundary contacts;
- positive and negative strand invariance;
- repeated-symbol physical copies remaining separate;
- pseudogene and ncRNA retention;
- p/q boundary assignment;
- exact one-to-one real GFF3 provenance and absence of propagated loci;
- required real-data row counts; and
- byte-for-byte rebuild reproducibility.

Validated audit invariants are also executable preconditions in
`build_inputs.py`: 37 PHRs, 61,312 unique-ID gene rows, 402 midpoint loci, 412
overlap loci, canonical contigs, nonoverlapping targets, valid coordinates,
and unambiguous target assignments. A source change that alters any audited
total fails loudly rather than silently changing the analysis universe.
