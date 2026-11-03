# Roadmap: integrating the heterolog/coalescence work into the manuscript

Assembled 2026-11-03. This file is the durable record of what is done, what is open, and in what order
- so that a compacted context cannot lose it again. Nothing here is applied to `paper.tex`; manuscript
text is produced as patch files for review.

## A. Artifacts that are done and current

| artifact | content |
|---|---|
| `response/notes/hprc_heterolog_divergence.md` | cohort-scale allelic vs paralog divergence; 2x2 classes; intra/inter community; best-partner per end; concert ladder; FINAL NUMBERS section; corrections log |
| `response/notes/community_tier_table.md` | the three tiers (allelic-like C7/C15/C14; near-allelic cohesive C1/C11/C3; weak; residual C4) |
| `response/notes/community_coalescence_times.md` | T = d/2mu bound (5-12 Myr upper, 2-5 Myr at mu=3e-8); distributions; skew |
| `response/notes/telomere_distance_profile.md` | ratio vs distance: step structure, breakpoints ~10 and ~20 kb, paralog maximum at 12-14 kb |
| `response/notes/telomere_gene_vs_intergenic.md` | gene vs intergenic split; the terminal ~6 kb is unmeasurable |
| `response/notes/9q3q_locus_full_report.md` | the locus record: alternatives A-G, ONT reads, dice switch, flagger, corrections |
| `response/notes/M2b_9q3q_read_support.md` | companion read/multi-alignment note |
| `response/notes/lit_Tan2026_methods_extract.md` | verified full-text extract of the companion paper; both deference errors retracted |
| `response/notes/heterolog_scope.pdf/.png`, `telomere_profile.pdf/.png`, `coalescence_distribution.pdf/.png` | figures |
| scripts | `scripts/ci/hprc_heterolog_divergence.py`, `telomere_distance_profile.py`, `telomere_gene_profile.py`, `telomere_end_ratio_curves.py`, `make_heterolog_figure.py`, `make_telomere_profile_figure.py`, `scripts/cladistics/community_coalescence.py` |

## B. Open items, ordered

### Batch 1 - deliverables that need no new compute

1. **Extended Data table from the community tier table.** Three tiers, ratio, n arms, n pairs, out-ratio,
   with C4 flagged as residual. Source of truth: `community_tier_table.md`.
2. **Index the orphan notes.** Six of seven notes are referenced by nothing in `README.md`,
   `POINT_BY_POINT.md` or `PROVENANCE.md`.
3. **Propagate the C4 demotion** to anything that lists communities (ED table caption, any figure that
   colours by community).

### Batch 2 - the new figure (compute running: `proc end-ratio-curves`)

4. **Spaghetti ratio-vs-distance figure.** Design:
   - **Per-end curves**: one curve per individual chromosome end (individual x arm), x = distance from
     the telomere-subtelomere junction (log), y = paralog/allelic divergence ratio. From
     `telomere_end_ratio_curves.py` (END level, weighted ratios). Draw a light-grey subsample (~300) for
     texture, not all ~19k.
   - **Overlaid summaries**: per-arm median curves, then per-community smooths (C3, C11, C1, C7, C15,
     C14) on top, plus the cohort median.
   - **Orientation**: distance outward from the telomere, so every end is normalised to the same zero;
     p and q arms are directly comparable because both are measured from their own telomere.
   - **MANDATORY masking**: the first ~6-10 kb must be hatched/greyed as unmeasurable (allelic
     divergence 0.3-2.5% there versus 0.09-0.16% interior, which cannot be orthologous). The JSON
     carries per-bin allelic values so the plotting script masks automatically. Plotting those points
     as measurements would be the single most misleading thing we could do with this figure.

### Batch 3 - correctness and bounded scientific loose ends

5. **Tip control** (long >= 5 kb and identity >= 95% records only): decides whether anything inside the
   terminal few kb survives. Cheap flag plus one run.
6. **chr18_q candidate false negative**: it has a 0.40% partner (chrX_p, 3.3 Myr) yet carries no PHR
   call. Check whether the caller's window rule excludes it and report either way.
7. **Copy-resolved test for WASH/DDX11L** - not feasible with multi-mapping assemblies; record as a
   stated limitation rather than an open promise (already written into
   `telomere_gene_vs_intergenic.md`).

### Batch 4 - manuscript text, as patches under `response/patches/` (never applied directly)

8. **Methods** for everything above: aggregator (binning, lexicographic dedupe, CIGAR d_sub, self-pair
   exclusion, inter-chromosomal-only), class definitions, distance bins and the both-sides rule, gene
   mask projection, mu and the upper-bound argument, and the discarding rule for invalid allelic bins.
9. **Results**: allelic baseline 0.15%; paralog/allelic 5.93x; intra vs inter chromosomal 1.07x;
   the three tiers with C4 demoted; best partner per end (0.52% median, 11/45 within 2x); the
   coalescence bound; the distance-profile step structure.
10. **Discussion**: combine with Tan et al. - they establish population frequency and qualitative
    antiquity; we establish inheritance at haplotype/molecule resolution, a molecular-clock bound, and
    the measurement limit near the telomere. State explicitly that neither subsumes the other.

### Batch 5 - previously-flagged manuscript defects still outstanding

11. **Fig 5B axis mislabel**: the panel plots dice (node overlap) but the caption says "best-match
    sequence identity". Real defect, not cosmetic.
12. **Four earlier pending prose items**: concerted-evolution passage; specificity paragraph (telomeric
    0.504 vs interstitial 0.129, delta +0.393, p = 0.006); censoring correction with the agreed framing
    (3.51 Mb excludes acrocentric p and X arms; no constituent PHR > 295 kb; lengths on the six
    saturated arms are lower bounds); Methods consolidation (dyadic clustering, average-rank ties, tool
    pins, PCLAI scoping, data-integrity statement).
13. **`data/chm13.phrs.bed` is the wrong call** (15,088 regions / 38 arms, `chr13_p` absent) and it
    feeds Fig 2b/2c-adjacent side tables and the permutation controls - regenerate before reuse.
14. **Open editorial question**: which Fig 5 panel to drop (A pedigree, B identity zoom, C ribbon), and
    whether the chr1p/chr5q dot plot becomes an Extended Data figure framed as "a second candidate of
    the same class, resolved as inherited by the parental-haplotype test".

### Batch 6 - bookkeeping

15. Extend `scripts/ci/audit_response_numbers.py` to machine-check the new notes (it has already caught
    one real error).
16. Flip the stale `to write` entries in `POINT_BY_POINT.md`, add the unsynced lanes, then `compile.sh`
    and the audit together.
17. Tan Methods + reference list once the PDF can be obtained (browser download); only then can the
    citation check and the rate-comparison sentence be written.

## C. What the spaghetti figure can and cannot show

**Can:** whether the paralog/allelic ratio is systematically lower near the junction for every end or
only for some; whether the tier structure (C7/C15 near 1, C3/C11 near 3.5, weak communities higher) is
visible end by end; the spread between ends at a given distance, which no pooled average currently
shows.

**Cannot:** say anything about the first few kb, because the denominator fails there. And it cannot
substitute for copy-resolved data at WASH/DDX11L, because per-end curves still compare window
sequences, not individual repeat copies.
