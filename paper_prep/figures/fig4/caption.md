# Figure 4 — Pedigree-resolved exchanges and cross-species generalisation

Direct observation of the events that build the communities, in human
pedigrees and across mammals.

**(a)** WashU 3-generation T2T pedigree (Cechova et al. 2025): per-flank
`odgi untangle` inter-chromosomal HQ patches for PAN027 mat ← PAN010, PAN027
pat ← PAN011, PAN028 mat ← PAN027. Coloured by Leiden community of the query
arm; **494/538 = 92 %** sit in a Leiden community
(`all_pedigrees_patches.tsv`).

**(b)** CEPH1463 4-generation pedigree (Porubsky et al. 2025): **11**
parent×chr-pair features detected by hifiasm AND verkko in ≥1 child each,
same Leiden community. chr10/chr18 (C2; Linardopoulou) in NA12877 paternal
AND NA12878 maternal independently; chr12/chr9 (C5) in NA12889 + NA12890 G1
grandparents (SURVEY_14 §1.6).

**(c)** RPE-1 (single diploid): 46-arm Jaccard (`rpe1.dist_matrix.tsv`) and
async CiFi Hi-C contact at 50 kb
(`rpe1_self_async_cifi_contact_matrix.tsv`). Red boxes mark **chrX_q × chr10_q**:
Leiden C2 = {chr10_q, chrX_q} and elevated 3D contact — **t(X;10)
rediscovered from sequence alone in one individual** (SURVEY_09 §1.3).

**(d)** Mouse (B6 + CAST T2T, Francis et al. 2025) zygotene Hi-C (Zuo et al.
2021, 50 kb). Per-PHR-pair Jaccard vs contact: **Spearman ρ = 0.715,
p = 4.4 × 10⁻⁵⁵, n = 344 inter-chr pairs**
(`zuo2021_zygotene_phr_pair_correlation.tsv`). Same correlation across all
4 meiotic stages (ρ 0.574–0.715; SURVEY_08 §1.7).

See `sources.tsv` for input paths.
