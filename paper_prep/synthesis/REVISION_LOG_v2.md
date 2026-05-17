# Revision log: draft-v2 pass B content additions + REFERENCES_v4 integration

Applied by agent-112 on 2026-05-17. Inputs: `CITATION_UPGRADE_PLAN.md`, `LITERATURE_REFRESH_v1.md`, `CONSISTENCY_AUDIT_v1.md`. Output file: `paper_prep/synthesis/NATURE_DRAFT_v2.md` (renamed from `NATURE_DRAFT_v1.md` so v1 freeze remains available for diff).

Word count: main text 3406 -> 3937 (delta +531; within the [3000, 4800] target). Abstract sha unchanged (214 words). Bibkey universe: 130 unique inline keys, 130 References-section entries, 0 leaks vs `REFERENCES_v4.bib` (364 entries).

## A. Citation upgrades from CITATION_UPGRADE_PLAN.md (per paragraph)

- P1 (history) -- ADD `@Salsi2026fshd` on the D4Z4 clause with a hedging clause ("the canonical D4Z4 macrosatellite ... with degenerate D4Z4-like copies on at least ten additional chromosomes revealed by T2T-CHM13"); ADD `@acrocentric_rdna_robertsonian_hartley2026biobank` on the acrocentric clause; ADD `@logsdon2025hgsvc` after the HPRC v2 citation block. `@yang2025chr2fusion` skipped (P1 does not currently invoke the chr2 fusion clause).
- P2 (methods substrate) -- ADD `@andreace2023pangenome`, `@heumos2024nfcore` on the PGGB clause.
- P3 (heatmaps, scale) -- no changes per plan.
- P4 (NJ tree) -- ADD `@hebbar2026marmoset`, `@degennaro2026ape` on the acrocentric clade clause; ADD `@Salsi2026fshd` on the 4q/10q clade clause (now reads "canonical D4Z4 macrosatellite").
- P5 -- no changes per plan.
- P6 (F_ST) -- ADD `@jeong2025segdup`, `@porubsky2026chr22q11`, `@hprc_siren2025`, `@bird2023africa` on the out-of-Africa F_ST citation block. P6 non-AFR F_ST range "-0.05 to 0.01" verified against `04_heterogeneity.md` L103-111 (corroborates Pass A edit 6); no further change.
- P7 (3D Hi-C) -- ADD `@hic3d_cheng2024`, `@bouquet_KaiserCTCF2025`, `@hic3d_kitamura2025` to the reproducibility-heatmap sentence.
- P8 (flanking paradox, bouquet) -- ADD `@bouquet_GarnerKASH52023`, `@bouquet_LiuSPDYA2025`, `@bouquet_MengSUN1NOA2023`, `@bouquet_JimenezCentromere2025` to the bouquet citation block; ADD `@subtelstruct_Smeds2025nonBDNA` to the D4Z4-CTCF-lamin clause.
- P9 (pedigree exchange) -- ADD `@pedigree_Schweiger2024spermNCO`, `@noyes2026sd`, `@chen2025paraphase` on the 133 gene-conversion-like clause; ADD `@Tardy2026fshd` on the C1 D4Z4 crossover-like sentence; ADD `@Porubsky2025` no-crossover-SV defusing sentence.
- P10 (RPE-1, mouse) -- ADD `@hic3d_Volpe2025RPE1` to the RPE-1 reference assembly clause; ADD `@t2t_Zhang2025macaque` to the mouse pipeline clause.
- P12 (causal loop, limits; was P11 in v1 before the gene-enrichment insertion) -- ADD `@palsson2025recomb` as a SUPPLEMENT to `@Sasani2019; @Smolka2024` on the Lalli-collapse clause and again in the outlook; ADD `@hic3d_Chen2026HiChew`, `@pangenome_Loegler2025review` to the methods outlook.
- Methods -- ADD `@logsdon2025hgsvc` (M1 samples); `@kaushan2026tracepoints` (M4 IMPG); `@andreace2023pangenome`, `@heumos2024nfcore` (M6 PGGB); `@hprc_siren2025` (M11 F_ST); `@hic3d_kitamura2025` (M14 single-cell 3D); `@Tan2018` left intact.

Net inline citation count: 101 unique keys (v1) -> 130 unique keys (v2). Net adds: 29 distinct new bibkeys; no REMOVE applied. `@Porubsky2025` is the in-place updated v3 entry now resolved to the published Nature 643:427-436 (handled inside `REFERENCES_v4.bib`).

## B. Gene-enrichment paragraph inserted between P10 (RPE-1/mouse) and P11 (causal loop)

New ~150-word paragraph. Names all six required entities:
- OR4F olfactory receptor family (Mefford canonical prediction, recovered in 7 of 15 communities; 10 OR4F members; OR4F5 and OR4F8P on 14 arms each). Source: `03_gene_enrichment.md` L33-39.
- Hub pseudogenes: RPL23AP45 (10 communities, 21 arms), SEPTIN14P22 (9 communities, 22 arms), DDX11L16 (9 communities, 20 arms). Source: `03_gene_enrichment.md` L94-100.
- MTCO pseudogenes (MTCO1P34, MTCO3P26/33/34) specifically enriched in acrocentric short-arm community C7. Source: `03_gene_enrichment.md` L70.
- SHOX as a PAR1 (C15) protein-coding gene. Source: `03_gene_enrichment.md` L76.
- Biotype composition: 32.1% protein-coding in C15 (PAR1) vs <=9% in the other 14 communities. Source: `03_gene_enrichment.md` L8-31.
- 11 Ambrosini subtelomere-specific duplicon blocks mapping into 15 Leiden arm-level communities. Source: `03_gene_enrichment.md` L48-60.
Cites `@Ambrosini2007`, `@Trask1998`, `@Mefford2001`, `@MeffordTrask2002` (all present in `REFERENCES_v4.bib`).

## C. C4 minimal-PHR positive-control sentence

Inserted in P7 (3D Hi-C paragraph), at the end of the reproducibility-heatmap clause, before the boundary to the multi-mapping paragraph. Verbatim: "The C4 community (chr7_q paired with chr12_q) is a minimal-PHR positive control: only 5 to 25 kb of tip-only shared sequence with zero gene annotations, yet significant 3D enrichment in 4 of 5 diploid Hi-C samples falsifies the alternative that the within-community signal is gene-content-driven."

Sources: `11_summary.md` L46-48 (finding 10) and `05_hic_validation.md` per-community enrichment table.

## D. Limitations clause: three new clauses appended

Original 4 limitations -> 7 limitations (one continuous sentence-clause structure preserved; no bullet list). New clauses (Pass B inserts):

1. (Second clause) Somatic exchange during cell-line propagation could inflate cross-arm affinity at the C1 4q/10q D4Z4 locus in LCL-derived HPRC v2 assemblies; matched blood-derived controls not yet available. Source: `10_limitations.md` limitation 8 (L23).
2. (Fourth clause) GM12878 Dip-C uses 16 cells after explicit exclusion of cell 12 as a duplicate of cell 10 (shared SRR7226706 long-insert library). Source: `10_limitations.md` limitation 17 (L43). NOTE: the task description stated "n=15 cells, not 16" but the report consistently says 16 cells used after exclusion (`06_dipc_validation.md` L7 confirms); v2 follows the report value.
3. (Fifth clause) PBMC Dip-C negative control inherits hg19-to-T2T coordinate projection noise (PHR boundaries projected via impg) which contributes to the non-significant W/B = 0.983 alongside small N = 18 and mixed cell-type composition. Source: `10_limitations.md` limitation 12 (L33).

## E. Mouse window-size sentence

Added in P10 (mouse paragraph), immediately before the rho = 0.715 Spearman correlation. Verbatim: "A window-size sweep showed 30 of 49 mouse PHRs saturate the 1 Mb extraction window, 25 of 49 at 2 Mb and 19 of 49 at 4 Mb, supporting the 1 Mb choice for the published rho = 0.715 correlation." Source: `08_mouse.md` L122-128 (PHR size comparison table, "Saturated (>=90% window)" column).

## F. RPE-1 self-vs-HPRC sentence

Added in P10 (RPE-1 paragraph), immediately after the t(X;10) sentence. Verbatim: "Mantel correlation between RPE-1 async CiFi Hi-C and the self-discovered partition is 0.548, vs 0.457 when the HPRC v2 population partition is transferred to RPE-1 without re-derivation: the population partition transfers to a single individual at a similar rho." Source: `09_rpe1_self.md` L56-71 (self vs HPRC table). NOTE: the task description listed the values 0.457 / 0.548 swapped; the report has HPRC-partition transferred async CiFi rho = 0.457 and self-discovered async CiFi rho = 0.548. v2 follows the report.

## G. PBMC Dip-C negative result

Added in P8 (3D paragraph), immediately after the S_all 11%/40% sentence. Verbatim: "A non-GM12878 Dip-C control on 18 PBMC cells [@Tan2018] gives W/B = 0.983, p = 0.305 at PHR-specific hg19-projected coordinates, confirming that the within-community signal is not a generic Dip-C artefact." Source: `05_hic_validation.md` L455-469.

## H. Internal report inconsistency reconciliation

Per instruction, the draft does not cite a specific count of findings (avoids the upstream `11_summary.md` "10 vs 12 findings" mismatch) and does not cite a specific count of novel contributions (avoids the upstream `12_literature.md` "24 vs 27" mismatch). No edit needed in the draft text; the existing v1 phrasing already avoids both counts.

## I. Verification

- All `[@bibkey]` clusters resolve to entries in `REFERENCES_v4.bib`: 130 unique inline keys vs the 364-entry bib; `comm -23 inline.txt bib.txt` returns empty.
- References section in the draft updated to list all 130 keys.
- Word count main text 3937 (target [3000, 4800]).
- Abstract sha unchanged (`diff` returns identical content; 214 words).
- Frontmatter updated: `version: v2`, `main_text_words: 3937`, `references: paper_prep/synthesis/REFERENCES_v4.bib`, `citation_plan: paper_prep/synthesis/CITATION_UPGRADE_PLAN.md`, `lit_refresh: paper_prep/synthesis/LITERATURE_REFRESH_v1.md`, `pass_history` entries added.
- Zero em-dashes in the draft.
- Zero `---` outside the YAML fences (lines 1 and 20 only).
- Pass A edits persist: `(chr5p, chr6q, chr7q, chr12q, chr14q, chr20p, chr20q)` does not appear; `12 of 15` (UPGMA agreement) and other Pass A fixes intact (Pass A's `REVISION_LOG_v1.5.md` unchanged).
