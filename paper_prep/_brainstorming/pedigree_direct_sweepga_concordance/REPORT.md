# Direct sweepGA parental haplotype concordance report

Generated: 2026-06-20.

## Scope

This package tests direct `sweepGA`/FastGA alignments from child haplotype windows to each transmitting-parent haplotype, with graph/`odgi untangle` outputs treated as the comparison target. The analysis uses native WashU assembly/window coordinates from FASTA record names and the graph-derived candidate tables; no CHM13/reference-projected coordinates are introduced.

## Execution

Heavy alignments were submitted through Slurm only. The completed raw run is:

| comparison | job | state | output |
|---|---:|---|---|
| PAN027 paternal hap2 vs PAN011 hap1 | 1704265 | COMPLETED | `raw_paf/PAN027pat_vs_PAN011_hap1.sweepga_many_many_j0.paf.gz` |
| PAN027 paternal hap2 vs PAN011 hap2 | 1704266 | COMPLETED | `raw_paf/PAN027pat_vs_PAN011_hap2.sweepga_many_many_j0.paf.gz` |
| PAN027 maternal hap1 vs PAN010 hap1 | 1704267 | COMPLETED | `raw_paf/PAN027mat_vs_PAN010_hap1.sweepga_many_many_j0.paf.gz` |
| PAN027 maternal hap1 vs PAN010 hap2 | 1704268 | COMPLETED | `raw_paf/PAN027mat_vs_PAN010_hap2.sweepga_many_many_j0.paf.gz` |
| PAN028 maternal hap1 vs PAN027 hap1 | 1704269 | COMPLETED | `raw_paf/PAN028mat_vs_PAN027_hap1.sweepga_many_many_j0.paf.gz` |
| PAN028 maternal hap1 vs PAN027 hap2 | 1704270 | COMPLETED | `raw_paf/PAN028mat_vs_PAN027_hap2.sweepga_many_many_j0.paf.gz` |

Command shape:

```bash
/home/erikg/.cargo/bin/sweepga --fastga --num-mappings many:many --scaffold-jump 0 --temp-dir "$TMPDIR" --output-file "$TMPDIR/<comparison>.paf" <query.fa> <target.fa>
```

The installed binary is `/home/erikg/.cargo/bin/sweepga` (`/export/local/home/erikg/.cargo/bin/sweepga`), reporting `sweepga 0.1.1`. The current binary supports `--num-mappings many:many --scaffold-jump 0`. It does not stream output with `--output-file -`, so the final wrapper writes a scratch PAF and compresses it explicitly.

## Output tables and plots

- Raw first-class PAFs: `raw_paf/*.sweepga_many_many_j0.paf.gz`
- Derived PAF filters: `filtered_paf/*.paf.gz`
- Per-PAF summary: `summaries/paf_file_summary.tsv`
- Graph-vs-direct segment concordance: `summaries/direct_vs_graph_concordance.tsv`
- Focused review plot: `plots/direct_sweepga_concordance_focused.{svg,pdf}`
- Full-window overview plot: `plots/direct_sweepga_full_genome_overview.{svg,pdf}`

## Concordance summary

`summaries/direct_vs_graph_concordance.tsv` compares every row of `fig5_synteny_recombination_schematic/selected_segments.tsv` against the direct PAFs for the corresponding parent-haplotype comparison. A row is called `agree` when the direct alignment overlaps the graph-derived query interval and lands on the same parent haplotype and target arm. It is `discordant_target` when the best direct overlapping signal goes to another haplotype/arm, and `inconclusive_partial_overlap` when the expected target is seen but the overlap is too small for a clean segment-level replacement.

Counts across selected graph-derived segments:

| event | agree | discordant | inconclusive | interpretation |
|---|---:|---:|---:|---|
| `PAR1_XY_positive_control` | 6 | 0 | 0 | Direct sweepGA cleanly recovers the expected PAR X/Y inheritance sanity-check structure. Direct PAF can be primary for these segment calls. |
| `PAN027_chr9q_chr3q_PHR_candidate` | 8 | 2 | 0 | Direct sweepGA strongly supports the graph-derived PAN027 paternal chr9q/chr3q event, including same-chromosome context and primary chr3q donor calls. The discordant rows are a side-fragment/low-confidence-tail class, so graph context should remain primary for those small fragments pending manual review. |
| `PAN028_chr9q_chr3q_PHR_candidate` | 14 | 7 | 1 | Direct sweepGA supports the broad PAN028 maternal chr9q/chr3q structure and multiple same-chromosome/primary-donor rows. Several side-fragment and small donor rows are discordant or partial, so use direct PAFs as primary only for the clean agreeing segments and keep graph/untangle primary for ambiguous fragments. |

Overall row counts: 28 `agree`, 9 `discordant_target`, 1 `inconclusive_partial_overlap`.

## Evidence-source recommendation

The table includes `evidence_source_recommendation`:

- `direct_sweepga_can_be_primary_for_segment`: 28 rows. These are clean direct overlaps to the expected parent haplotype and target arm, usually with identity above 0.95.
- `keep_graph_primary_pending_manual_review`: 9 rows. The direct signal is present but the strongest overlap points to a different target arm/haplotype than the graph-derived small fragment.
- `inconclusive_keep_graph_primary`: 1 row. The expected target is visible but only partially overlaps the graph interval.

Recommended source policy for the downstream Fig. 5 evidence decision:

- PAR1 positive control: direct sweepGA can be primary.
- PAN027 paternal chr9q/chr3q candidate: direct sweepGA can be primary for same-chromosome context and primary chr3q donor segments; retain graph/untangle as primary for the small side/low-confidence fragments.
- PAN028 maternal chr9q/chr3q candidate: direct sweepGA can be primary for the broad same-chromosome context and clean chr3q donor segments; retain graph/untangle as primary for discordant side fragments and partial donor/tail calls.
- PAN027 maternal vs PAN010: raw and filtered direct PAFs were produced as required, but the current graph selected-segment table has no corresponding candidate event row, so this comparison has file-level summaries but no event-level concordance row.

## Validation notes

- Heavy alignment was not run on the head node; all six raw direct comparisons completed through Slurm.
- Raw many:many/no-scaffold PAFs are preserved and compressed.
- Filters were derived from raw PAFs without rerunning alignment.
- Coordinates remain native assembly/window coordinates from source FASTA names and graph-derived selected segment intervals.
- `submission/` and manuscript figure directories were not modified.
