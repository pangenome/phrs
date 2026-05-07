# Recommended v7 Slide Patch

Do not edit the deck in this task. This patch note is for `review-zoom-v7-fanin-render`.

## Target Area

The current deck has a gene-enrichment section around `slides/v2-review-zoom/_typst/zoom_review_deck.typ` slides 14m-14c:

- `14m`: method transition
- `14a`: copy-aware method boundary
- `14b`: ranked copy-aware candidate signals
- `14c`: community/family map

Erik's feedback is that this section underuses the methods/results report. For v7, replace or tighten that section with report-backed wording from `subtelomeric_analysis_report.md` section 9 and `end-to-end-report/report/03_gene_enrichment.md`.

## Recommended Slide 1

Slide number: keep near `14m` or use `14g.1` if preserving existing numbering.

Title:

> Gene enrichment: report-backed biology, conservative statistics

Body:

- Canonical HPRCv2 community-family Fisher screen: 116 tested rows, 0 BH-significant rows.
- C3 OR and C7 MTCO are near-miss/candidate presence patterns with BH q = 0.07118, not definitive enriched classes.
- The biology is still slide-worthy as recurrent gene-family architecture: OR4F/olfactory receptors, C1 D4Z4/DUX4L, C5 DDX11L/WASH/FAM138/IQSEC3, C15 PAR1 coding genes, C7 MTCO pseudogenes, and a pseudogene/ncRNA-rich backbone.

Suggested visual:

- Use `gene_enrichment_report_backed_summary.svg`, or render `gene_family_function_slide_table.md` as a compact table.

Source line:

> Sources: `subtelomeric_analysis_report.md:477-600`; `end-to-end-report/report/03_gene_enrichment.md`; HPRCv2 enrichment TSV rows listed in `v7/gene_enrichment_report_backed/README.md`.

Speaker note:

> The main correction is statistical language. The canonical family Fisher screen does not give BH-significant results, so this is not a GO-proof slide. But section 9 gives strong report-backed gene architecture. We can say recurrent OR4F architecture, C1 D4Z4/DUX4L, C5 DDX11L/WASH/FAM138 plus IQSEC3, PAR1 coding genes, and C7 acrocentric MTCO candidate signal.

## Recommended Slide 2

Slide number: immediately after slide 1 if the section keeps two slides.

Title:

> GO context is useful, but not the proof layer

Body:

- Standard GO/ORA tables are historical context; they collapse copy number and should be labeled standard/deduplicated.
- Coding-only ORA recovers olfactory and GTP-related terms from a 9-gene query; useful context, but not the canonical HPRCv2 community-family analysis.
- Exploratory copy-weighted ORA gives striking fold enrichments, but validation notes anti-conservative p-values and failed FDR control unless revalidated with calibrated/permutation methods.
- Therefore: lead with report-backed HPRCv2 community-arm support; show GO/ORA only as method contrast.

Suggested visual:

- Use `go_function_context_table.md` as a small table.

Source line:

> Sources: parked `paper_prep/_brainstorming/phr_GO_*`, `phr_coding_only_GO_*`, `phr_copy_weighted_enrichment.csv`, and `improved_copy_weighted_enrichment.csv`; caveat from `copy_number_weighted_ora_methodology_synthesis.md`.

Speaker note:

> The older GO outputs are useful because they explain why OR and GTP terms came up during brainstorming. They are not the final statistical evidence. Standard and coding-only ORA collapse copy number; the improved copy-weighted table is explicitly exploratory and includes a known smell/olfactory mapping caveat where the listed genes are IL9R/IL9RP, not OR4F.

## Optional Compression

If only one slide can be added, use the title and body from Recommended Slide 1, then add one small note at the bottom:

> GO/ORA context: standard/coding-only outputs are historical copy-collapsed tables; exploratory copy-weighted outputs are not final statistics without permutation/calibration.

## Exact Phrases To Use

- "recurrent gene-family architecture"
- "copy-aware candidate signal"
- "report-backed community-arm support"
- "presence pattern that does not survive BH correction"

## Exact Phrases To Avoid

- "BH-significant GO enrichment proves..."
- "definitive enriched class"
- "validated copy-weighted ORA p-value"
- "OR4F GO enrichment from the current canonical HPRCv2 family tests"

## Integration Notes

- Keep final deck integration in the fan-in render task.
- It is safe to keep the existing D4Z4 and OR4F browser backups, but speaker notes should match the v7 caveat language.
- If replacing the v5 figures directly, use `gene_enrichment_report_backed_summary.svg` as the fastest report-backed visual and cite the README.
- If preserving v5 figure 14b, update its title or source note so bar length is explicitly support count, not q-value or validated enrichment effect.
