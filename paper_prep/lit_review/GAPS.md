# GAPS — open questions and missing-substrate notes

**Author:** lit-review-synthesis (agent-877), 2026-05-06.
**Companion files:** `SYNTHESIS.md`, `CHRONOLOGY.md`, `REFERENCES_v2.bib`, `paper_prep/synthesis/CROSSWALK.md`.

This is a one-page note flagging anything the synthesis pass surfaced as missing, mis-cited, or open. It is not a follow-up plan; only a list of items that future passes (or revisions to the manuscript) should resolve.

## Substrate provenance

All 14 topical reviews are now committed under `paper_prep/lit_review/`:
`topic_01` cytogenetic foundations (21 bib entries), `topic_02` TAR1 / TTAGGG / ITS structure (16), `topic_03` pseudohomologous regions concept (14), `topic_04` sex-chromosome pseudoautosomal regions (28), `topic_05` acrocentric / rDNA / Robertsonian (16), `topic_06` D4Z4 / DUX4 / FSHD (13), `topic_07` concerted evolution / NAHR / gene conversion (18), `topic_08` meiotic bouquet / nuclear-envelope tethering (19), `topic_09` Hi-C / Pore-C / single-cell / meiotic 3D (16), `topic_10` pangenome graphs / IMPG (15), `topic_11` pedigree-based recombination detection (20), `topic_12` HPRC v1 / v2 population pangenomes (17), `topic_13` subtelomere population genetics / Fst / out-of-Africa (12), `topic_14` OR / OR4F gradient (17).

`REFERENCES_v2.bib` merges `REFERENCES.bib` (24) plus all 14 topic bibs (266 raw, 242 unique after intra-topic dedup) into 270 deduplicated entries. Dedup is by DOI (when present) and by (first-author surname, year, normalized-title prefix) with the canonical / prefix-free key preferred on collision so existing `SYNTHESIS.md` citations resolve unchanged. 141 of the 270 entries use canonical keys (carried over from `REFERENCES.bib` + the prior synthesis pass that filled in canonical keys for milestone refs); the remaining 129 retain their topic-prefix keys (e.g. `acrocentric_*`, `bouquet_*`, `dux4_d4z4_fshd_*`) and represent unique adds from the topic reviews.

Earlier passes of this synthesis (agent-898, agent-913) operated when only 2 or 3 of the 14 topic reviews were committed, so much of the per-paper substrate had to be reconstructed from `REFERENCES.bib`, the `CROSSWALK.md` claim ↔ substrate map and the end-to-end-report literature ledger (`12_literature.md`). The current pass (agent-928) restores grounding by merging in the now-committed topical reviews; the narrative depth of Parts I–IV is competent at the field-history level, and the per-paper depth for any sub-topic lives in the corresponding `topic_NN_*.md` file. A future revision pass that wants to thread topic-specific citations into the prose at `topic_NN_*.md` density would primarily edit Parts I–IV in place, replacing canonical paper-anchored citations with the matching `topic_NN_<key>` substrate where richer and not yet integrated.

## Bibliography-level gaps

- **`Cechova2025`** is cited as "(in press)" with no DOI; the WashU pedigree paper is not yet on bioRxiv per the lead author's clarification. Update DOI once available.
- **`Porubsky2025`** is cited as "(in press)" with no DOI; the CEPH1463 paper was published in *Nature* (April 2025) per the lead author's recent comment. Update with the published DOI before submission.
- **`deLima2025`** is "(in submission)"; verify status before submission.
- **`Xu2025`** lists "Hanbo Xu and others" — replace with full author list once the manuscript is in its published form.
- **`Francis2025`** has full DOI (`10.1038/s41588-025-02358-1`) but author list is currently truncated to "Brittany A. Francis and others" — replace with full author list before submission.
- **Stergachis Fiber-seq telomere companion paper.** `REFERENCES_v2.bib` (and `REFERENCES.bib`) cite `Stergachis2020` with a note that the specific Fiber-seq-at-39-of-46-telomeres paper is the relevant downstream work — track that paper down before submission. Likely candidates: a 2022–2024 *Cell* / *Cell Genomics* / *Nature* Stergachis-lab paper extending Fiber-seq to centromeres and telomeres.
- **`StankiewiczLupski2002` and `StankiewiczLupski2010`** are both cited; verify the 2002 *Trends in Genetics* paper and the 2010 *Annual Review of Medicine* paper resolve to the intended NAHR vocabulary in context.
- **`Computational2018` (consortium-author key)** is the canonical citation for the *Computational pan-genomics* review (Marschall, Marz, Abeel et al., *Briefings in Bioinformatics* 2018, doi:10.1093/bib/bbw089). The duplicate author-led `Marschall2018` key was removed from `REFERENCES_v2.bib` to satisfy the deduplication rule.
- **Topic-prefix vs canonical-key deduplication.** Several topical reviews use prefixed keys (e.g. `sexchrompars_charchar2003`, `acrocentric_Henderson1972`, `dux4_d4z4_fshd_lemmers2010`) for papers that the existing `REFERENCES.bib` / prior synthesis pass cite as canonical PascalCase (`Charchar2003`, `Henderson1972`, `Lemmers2010`). The merge into `REFERENCES_v2.bib` keeps the canonical key on each duplicate (per the dedup rule: prefer prefix-free) so `SYNTHESIS.md` and `CHRONOLOGY.md` citations continue to resolve. The `topic_NN_*.md` files (read-only) still use their prefix keys and must therefore be compiled against `topic_NN_*.bib` directly, not against `REFERENCES_v2.bib` alone. The unique-to-topic prefix keys (~129 entries) survive in `REFERENCES_v2.bib` unchanged.

## Substantive gaps

- **C7 — direct LAD / Lamin B1 ChIP overlay.** The abstract's "facilitated by the physical proximity of subtelomeres at the nuclear envelope" wording rests on Dip-C radial position as a proxy and on the Masny 2004 / Ottaviani 2009 mechanism for D4Z4 — chapter 05 §"Nuclear lamina cross-reference" does not run a genome-wide LAD intersection. `REWRITE_PLAN.md` TASK-22 is the correct follow-up. Fallback per `CROSSWALK §6`: soften to "consistent with peripheral telomere positioning."
- **C5 — neighbour-joining tree** is named in the abstract but Andrea uses Leiden + UPGMA (12 of 15 communities agree exactly). One-line fix is `ape::nj()` on the existing arm-distance matrix; alternatively, relax the abstract wording to "cladistic analysis (Leiden / UPGMA)." Per `REWRITE_PLAN.md` TASK-01 / TASK-13.
- **C2 — "~12% pairwise sampling" computation.** Methods must compute the wfmash k-mer-evaluation rate from the on-disk PAFs (the realised value, not asserted) and write the Erdős-Rényi argument explicitly (n=18,827; p* = log(n)/n ≈ 5.21e-4; 12% is ~230× above). Per `REWRITE_PLAN.md` TASK-10 / TASK-11. The argument is in `CROSSWALK §7b`.
- **C2 / C3 dataset count.** Andrea reports 465 throughout; the abstract reports 466. Resolve canonically as 466 = 233 individuals × 2 haplotypes + CHM13 reference (per `CROSSWALK §7c`).
- **C8 — "concerted evolution" loose-sense framing.** The Discussion paragraph should explicitly state the loose-sense usage and tie it to the pedigree (chapter 14) as the empirical anchor. Per `CROSSWALK §3` and `REWRITE_PLAN.md` TASK-19.

## Open empirical questions surfaced by the synthesis

- **Crossover-rate ↔ cross-arm-affinity correlation honest negative.** End-to-end-report chapter 12 testable prediction #7 reports `rho = −0.43, p = 0.006` across all 39 arms but `rho = 0.00` after excluding 7 short-read-confounded arms. Open question: does a higher-quality recombination map (deeper sequencing, longer reads at acrocentric / PAR arms) recover the predicted negative correlation? Probably yes, but cannot be answered without new data.
- **LINC-complex / SUN1-mutant test of meiotic alignment at PHR scale.** Zuo 2021 wild-type vs SUN1-W151R zygotene Hi-C is available (GEO: GSE155142, GSE155638, GSE155967). The PHR-median 105 kb scale is well below the ~5% mutant tip zone; tip contacts may be maintained without LINC-mediated force transmission, so the predicted effect is uncertain. Per chapter 12 testable prediction #1.
- **Somatic vs germline exchange.** HPRC v2 sample-source metadata (LCL vs blood) is not currently available; the manuscript cannot directly distinguish meiotic from somatic-LCL exchange in chapter 14 patches. The pedigree's *inheritance* across three generations does rule out purely somatic origin for the transmitted patches but does not directly establish meiotic timing in cells.
- **Subtelomeric Fiber-seq / nucleosome footprint at PHR boundaries.** The Stergachis 2020 framework supports single-molecule CTCF / nucleosome footprinting; a per-PHR-boundary analysis at 39 / 46 telomeres would either confirm or reject the CTCF-cohesin boundary prediction (chapter 12 testable prediction #4). Not currently in scope but a clean follow-up.

## What this synthesis does *not* claim

This document does not extend any of the abstract's claims beyond what `paper_prep/synthesis/CROSSWALK.md`, `paper_prep/synthesis/REFERENCES.bib`, the 14 committed topical reviews (`topic_01_*` through `topic_14_*`) and `end-to-end-report/report/12_literature.md` already establish. The prose in `SYNTHESIS.md` operates at the field-history-and-framework level; where deeper per-paper analysis is needed, the corresponding `topic_NN_*.md` provides it. No new analytical claim is made; this is a literature synthesis only.

*End of GAPS.md.*
