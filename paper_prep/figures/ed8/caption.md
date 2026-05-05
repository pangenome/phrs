# Extended Data Figure 8 — Discussion synthesis: feedback loop, D4Z4 model, recombination null, compartment diagnostic

**a, Causal feedback loop (`SURVEY_07 §1.6`).** Four-link cycle: sequence sharing → 3D proximity → ectopic exchange → new shared segments → propagation. Edges encode support: solid blue = direct measurement (CHM13 Hi-C ρ = 0.674; HG002 Pore-C ρ = 0.485); solid green = established literature (FSHD D4Z4 translocations, Lemmers 2010); dashed olive = inferred.

**b, D4Z4-CTCF-lamin tethering for C1 (chr4q ↔ chr10q).** D4Z4 macrosatellite at the tip → CTCF binding (Ottaviani 2009) → lamin A/C tether (Masny 2004) → co-peripheralization (Dip-C radial 0.732). C1: median 22 DUX4L vs 0–2 in 7 outliers (Mann-Whitney p = 5.3e-6); 43.4 % chr4q discordance. Inset: inter-arm sharing peaks at 0–15 kb (D4Z4 location).

**c, Recombination vs cross-arm affinity — honest null.** Lalli 2025 cM/Mb (`subtelomeric_recomb_rates.tsv`) vs per-arm fraction of sequences in cross-arm communities (from `cross_arm_affinity_sequences.tsv` + assignments). Left: all 46 arms, ρ = −0.35, p = 0.017. Right: well-callable arms only (>12 callable variants; n = 40), ρ = −0.01, p = 0.97. Anti-correlation is driven by seven low-callability arms (chr14p/22p/Xp = 0; chr15p = 1; chr13p = 9; chr21p = 11). Survey-reported published values: ρ = −0.43 (39 arms) → ρ = 0.00 (N = 32). Long-read maps required.

**d, Compartment identity at tips (HG002, n = 92 arm × hap).** 63 / 92 = 68 % A-compartment (e1 > 0); mean e1 = +0.0073 — weak signature. Source: `compartment_analysis.tsv` (HG002 100 kb Hi-C eigenvector, GC-oriented). Dip-C radials are bimodal (C1 = 0.732 peripheral; C6 / C10 = 0.505 / 0.474 interior): tips are A-leaning by GC but internally positioned, dissociating compartment identity from envelope tethering.

(Word count ≈ 200.)
