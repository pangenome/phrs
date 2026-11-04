# Materials index: every analysis in the heterolog / coalescence / telomere cycle

Purpose: a complete, auditable map from question -> script -> input -> output -> note -> figure -> status,
so that each analysis contributing to the response can be located, re-run and cited. Companion to
`ROADMAP_analysis_integration.md` (what remains) and `PROVENANCE.md` (data lineage, P1-P18).

All scripts are in the repo. Large intermediate outputs are on `/moosefs/erikg/phrs_tmp/` and are
reproducible from the scripts; none are committed.

## A. Divergence and the allelic/paralogous framework

| analysis | script | input | output | note |
|---|---|---|---|---|
| cohort allelic vs paralogous vs community divergence (1e-4 bins, self-pairs excluded, both d_sub and gap-compressed) | `scripts/ci/hprc_heterolog_divergence.py` | `/moosefs/.../PHR_III/all-vs-all.1Mb.p95/*.paf.gz` (18,827 PAFs) | `hprc_het_final.json/.tsv`; committed pair table `scripts/ci/hprc_heterolog_divergence_pairs.tsv` | `hprc_heterolog_divergence.md` |
| community tiers | (derived from the above + k=15 assignments) | as above | `submission/supporting_material/ED_community_tier_table.tsv` | `community_tier_table.md`; figure patch `patches/ED_table_community_tiers.tex` |
| ED Fig 5 (scope figure, 5 panels) | `paper_prep/figures/ed5_paralog_allelic/figure_ed5_paralog_allelic.py` | `hprc_het_final.json` | `figure_ed5_paralog_allelic.pdf/.png` -> `submission/fig/ExtendedDataFigures/ED_Fig5_*` | caption in the figure dir; patch `patches/paper_tex_insert_ed5.tex` |
| exchange rate / transfer size / R-to-rate mapping | (arithmetic in the note; no separate script) | derived | - | `exchange_rate_estimates.md` |
| coalescence times per community | `scripts/cladistics/community_coalescence.py` | 500 kb windows of one genome, wfmash all-vs-all | `community_coalescence_communities.tsv` (committed), `_pairs.tsv` | `community_coalescence_times.md` |

## B. Positional structure toward the telomere

| analysis | script | input | output | note |
|---|---|---|---|---|
| ratio vs distance from the telomere (pooled) | `scripts/ci/telomere_distance_profile.py` | same PAF set | `/moosefs/erikg/phrs_tmp/telomere_profile.json` | `telomere_distance_profile.md` |
| gene vs intergenic split of the same profile | `scripts/ci/telomere_gene_profile.py` | PAF set + `data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz` | `telomere_gene_profile.json` | `telomere_gene_vs_intergenic.md` |
| per-arm and per-end ratio curves | `scripts/ci/telomere_end_ratio_curves.py` | PAF set | `telomere_end_curves.json` (48 arms, 14,814 ends) | `telomere_distance_profile.md` |
| ED Fig 6 (profile: per-arm, per-end, tiers) | `paper_prep/figures/ed6_telomere_profile/figure_ed6_telomere_profile.py` | `telomere_end_curves.json` + k=15 assignments | `figure_ed6_telomere_profile.pdf/.png` -> `submission/fig/ExtendedDataFigures/ED_Fig6_*` | caption in the figure dir; patch `patches/paper_tex_insert_ed6.tex` |
| terminal-zone diagnostic (telomeric repeat content, composition, CpG O/E, allelic alignability) | inline in this cycle (documented in the note with exact parameters) | window FASTA + PAFs | tables in the note | `telomere_distance_profile.md`, `telomere_gene_vs_intergenic.md` |
| gene content at the tip | inline (GFF3 + window coordinates) | `data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz` | tables in the note | `telomere_gene_vs_intergenic.md` |

## C. The pedigree locus

| analysis | script | input | output | note |
|---|---|---|---|---|
| ONT read support at all six involved ends | `scripts/pedigree/ont_read_support_9q3q.py` | ONT BAMs, `washu_v131/ont/` | `response/notes/read_support.json` (committed) | `M2b_9q3q_read_support.md`, `9q3q_locus_full_report.md` |
| multi-alignment class winners (dice) at 9q/3q and 5q/1p | (pipeline outputs) | `data/fig5_*.class_winners.impg_similarity.tsv.gz` | tables in the notes | `9q3q_locus_full_report.md` §6 |
| flagger/nucflag overlap QC over the 538 WashU patch calls | `patch_flagger_qc.py` (working copy) | flagger/nucflag BEDs (UCSC public) | `washu_patches_flagger_overlap.tsv` | `01_evidence_log_pedigree.txt` |
| pedigree patch untangling and copy correction | `scripts/pedigree/untangle_multimap_tracts.py --kid-correct-copies` | assemblies v1.1 | tract tables | `01_evidence_log_pedigree.txt` rev 10-12 |

## D. Comparanda and decisions

| item | artefact | note |
|---|---|---|
| companion paper (Tan et al. 2026) full-text extract | `lit_Tan2026_methods_extract.md` | retrieval route and limits recorded in PROVENANCE P17 |
| 1p/5q resolution | `1p5q_resolution.md` | |
| community tiers response paragraph (acrocentric/PAR objection) | `patches/response_acrocentric_PAR.tex` | |
| roadmap and outstanding work | `ROADMAP_analysis_integration.md` | |

## E. Reproduction rules that matter

1. `d_sub = X/(=+X)` from CIGARs only. The `md:f` tag is **not** substitutions/alignment length
   (52% disagreement with the CIGAR) and must not be used; `gi:f` is gap-compressed identity.
2. Bin at 1e-4. At 1e-3 every per-record median is quantised and the allelic baseline is wrong
   (0.10% instead of 0.15%).
3. Exclude allelic self-pairs (`arm_a == arm_b`) before any "best partner" search, or each arm matches
   its own alleles.
4. Telomere orientation: the window interval endpoint coinciding with a contig end is the telomere
   (p arms position 0, q arms position L); windows are stored in original contig orientation.
5. Mask the terminal ~6 kb: allelic alignments there are sparse and fragmentary (2 orders of magnitude
   fewer individuals, blocks ~18x shorter), so the ratio is undefined.
6. Report per-arm equal-weight medians for positional claims, not pooled values - pooling lets a few
   high-record arms drive the curve.
