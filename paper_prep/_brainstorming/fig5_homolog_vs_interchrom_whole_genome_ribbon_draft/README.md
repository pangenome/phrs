# Fig5 whole-genome 10:10 ribbon draft

Draft whole-genome ribbon view for the corrected `PAN027pat_vs_PAN011_joint`
SweepGA/F32 10:10 IMPG class-winner scan.

Inputs are vendored in `data/` so the figure is self-contained (no moosefs
needed); the script falls back to the moosefs sources, or override on the CLI:

- IMPG 10:10 class-winner scan (`--class-winners`):
  `data/fig5_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz`
  (14 MB; the `PAN027pat_vs_PAN011_joint.sweepga_f32.10to10.query_2000bp.predepth_class_winners.impg_similarity.tsv.gz` scan)
- child query chromosome lengths (`--query-fai`):
  `data/fig5_PAN027pat_vs_PAN011_joint.query.fa.fai`
- father joint-target chromosome lengths (`--target-fai`):
  `data/fig5_PAN027pat_vs_PAN011_joint.target.fa.fai`

Lane layout:

- top genome: PAN011 h1
- middle genome: PAN027 paternal haplotype, the meiotic product inherited by
  the child/query haplotype
- lower genome: PAN011 h2
- chromosomes are concatenated in chromosome order with actual chromosome-length
  scaling within each genome track
- colored ribbons show adjacent 2 kb query windows where the best interchromosomal IMPG
  match beats the best same-chromosome/homologous match
- drawn ribbons are filtered to runs of at least 10 kb and mean interchromosomal
  identity of at least 0.95

Regenerate:

```bash
# python3 (stdlib) reads the vendored data/ inputs; rsvg-convert (librsvg, via
# guix) is needed for PDF/PNG, otherwise it writes SVG only. The wrapper passes
# the published short track labels (PAN011 h1 / PAN027 paternal / PAN011 h2).
bash paper_prep/_brainstorming/fig5_homolog_vs_interchrom_whole_genome_ribbon_draft/scripts/make_whole_genome_ribbon_draft.sh
```

Label provenance: the source FASTA manifest reliably maps the child/query
records for this panel (`PAN027#2`) to the paternal haplotype, but the displayed
PAN011 parent target records are only represented as recovered source haplotype
numbers after joint-target collapsing (`PAN011#joint#h1_chr*` and
`PAN011#joint#h2_chr*`). Because the available manifest does not establish a
parental-origin mapping for PAN011 h1/h2, the parent tracks intentionally remain
`PAN011 h1` and `PAN011 h2` rather than paternal/maternal labels.

Note: the committed SVG/PDF/PNG were rendered by an earlier version of the label
logic; re-running the current script reproduces the same ribbons/coordinates but
differs by a few cosmetic label-text lines (the right-edge layer labels and a
comparison-id subtitle prefix).

Outputs:

- `fig5_homolog_vs_interchrom_whole_genome_ribbon_draft.svg`
- `fig5_homolog_vs_interchrom_whole_genome_ribbon_draft.pdf`
- `fig5_homolog_vs_interchrom_whole_genome_ribbon_draft.png`
- `fig5_homologous_recombination_context_ribbon_draft.svg`
- `fig5_homologous_recombination_context_ribbon_draft.pdf`
- `fig5_homologous_recombination_context_ribbon_draft.png`
- `whole_genome_ribbon_runs.tsv`
- `whole_genome_ribbon_summary.tsv`
- `whole_genome_homologous_context_runs.tsv`
- `whole_genome_homologous_context_summary.tsv`
- `whole_genome_ribbon_merge_audit.tsv`

The homologous-context variant uses the same geometry but adds full
same-chromosome father-child homologous chains as a light-gray ribbon layer.
Homologous chains are grouped from all `same_chrom` class-winner rows with
identity at least 0.95 and drawn when the grouped run is at least 10 kb. Each
light-gray homologous ribbon uses the exact native grouped donor interval and
exact native grouped child interval in the TSV outputs. For visibility on the
whole-genome track, the plotted homologous ribbon and endpoint marks apply a
display-only minimum width that scales with chain length; this keeps long
homology visually legible without changing the exact end-to-end merge or
reported coordinates. The colored interchromosomal/non-homologous winners are
preserved on top.

Raw 2 kb windows are grouped without gap-tolerant coalescing. The run builder
collapses windows with the same child sequence and donor sequence only when the
child endpoint and donor endpoint both touch exactly in a consistent donor
direction (`end_to_end_merge_gap_bp = 0`). This makes adjacent same donor-child
chains render as one wider ribbon while preserving all fragments in the run's
interval, base-pair total, and window count.
`whole_genome_ribbon_merge_audit.tsv` records the raw, end-to-end merged, and
drawn counts and the maximum absorbed query/donor endpoint gaps.
