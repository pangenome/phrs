# GAPS — open questions and missing-substrate notes

**Author:** lit-review-synthesis (agent-877), 2026-05-06.
**Companion files:** `SYNTHESIS.md`, `CHRONOLOGY.md`, `REFERENCES_v2.bib`, `paper_prep/synthesis/CROSSWALK.md`.

This is a one-page note flagging anything the synthesis pass surfaced as missing, mis-cited, or open. It is not a follow-up plan; only a list of items that future passes (or revisions to the manuscript) should resolve.

## Process gap: 12 of 14 topical reviews never committed

Of the 14 topical lit reviews originally fanned out, only `topic_03` (pseudohomologous regions concept) and `topic_11` (pedigree-based recombination detection) are committed to the repository. The other 12 — `topic_01` (cytogenetic foundations), `topic_02` (TAR1 / TTAGGG / ITS), `topic_04` / `topic_04_2` (PARs), `topic_05` (acrocentric / rDNA / Robertsonian), `topic_06` (D4Z4 / DUX4 / FSHD), `topic_07` (concerted evolution / NAHR / gene conversion), `topic_08` (meiotic bouquet / nuclear-envelope tethering), `topic_09` (Hi-C / Pore-C / single-cell / meiotic 3D), `topic_10` (pangenome graphs / IMPG / transitive closure), `topic_12` (HPRC v1 / v2 population pangenomes), `topic_13` (subtelomere population genetics / Fst / out-of-Africa), `topic_14` (OR / OR4F gradient) — were marked done by their dispatched agents in 30–60 s without producing any committed output. Their evaluator scores ranged 0.02–0.17.

Consequence for this synthesis: Parts I, III and IV draw their substrate from `paper_prep/synthesis/REFERENCES_v2.bib` (123 entries, augmented from the existing 24-entry `REFERENCES.bib`), the end-to-end-report's literature ledger (`end-to-end-report/report/12_literature.md`), Andrea's appendix references (`13_appendix.md`), the `CROSSWALK.md` claim ↔ substrate map, and the lead-author clarifications recorded therein. The depth is competent at the chronological-narrative level but does not go to the per-paper case-by-case depth that `topic_03` and `topic_11` do. A future pass redo of those 12 topics would substantially deepen Parts I–IV with primary-paper-by-primary-paper analysis.

## Bibliography-level gaps

- **`Cechova2025`** is cited as "(in press)" with no DOI; the WashU pedigree paper is not yet on bioRxiv per the lead author's clarification. Update DOI once available.
- **`Porubsky2025`** is cited as "(in press)" with no DOI; the CEPH1463 paper was published in *Nature* (April 2025) per the lead author's recent comment. Update with the published DOI before submission.
- **`deLima2025`** is "(in submission)"; verify status before submission.
- **`Xu2025`** lists "Hanbo Xu and others" — replace with full author list once the manuscript is in its published form.
- **`Francis2025`** has full DOI (`10.1038/s41588-025-02358-1`) but author list is currently truncated to "Brittany A. Francis and others" — replace with full author list before submission.
- **Stergachis Fiber-seq telomere companion paper.** `REFERENCES_v2.bib` (and `REFERENCES.bib`) cite `Stergachis2020` with a note that the specific Fiber-seq-at-39-of-46-telomeres paper is the relevant downstream work — track that paper down before submission. Likely candidates: a 2022–2024 *Cell* / *Cell Genomics* / *Nature* Stergachis-lab paper extending Fiber-seq to centromeres and telomeres.
- **`StankiewiczLupski2002` and `StankiewiczLupski2010`** are both cited; verify the 2002 *Trends in Genetics* paper and the 2010 *Annual Review of Medicine* paper resolve to the intended NAHR vocabulary in context.
- **`Computational2018` (consortium-author key)** is the canonical citation for the *Computational pan-genomics* review (Marschall, Marz, Abeel et al., *Briefings in Bioinformatics* 2018, doi:10.1093/bib/bbw089). The duplicate author-led `Marschall2018` key was removed from `REFERENCES_v2.bib` to satisfy the deduplication rule.

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

This document does not extend any of the abstract's claims beyond what `paper_prep/synthesis/CROSSWALK.md`, `paper_prep/synthesis/REFERENCES.bib`, the committed topical reviews (`topic_03`, `topic_11`) and `end-to-end-report/report/12_literature.md` already establish. Where the substrate is thin (Parts I, III and IV at the per-paper level), the prose in `SYNTHESIS.md` operates at the field-history-and-framework level and signposts where the missing topical reviews would have provided per-paper depth. No new analytical claim is made; this is a literature synthesis only.

*End of GAPS.md.*
