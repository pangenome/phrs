# Extended Data Figure 5 — Multi-resolution and confound-exclusion robustness for Hi-C

**a, W/B contact across five mcool resolutions, eight datasets.** Mean inter-chromosomal within-community contact ÷ mean between, from `community_based/{5,10,20,50,100}000bp/<sample>_global_test.tsv`. All 40 (sample × resolution) cells have W/B ≥ 11.6 (median 18.5, max 48.4); rank ordering of samples is essentially preserved across the 5 kb–100 kb axis — no single resolution drives the signal.

**b, Mantel ρ before vs after acrocentric+sex exclusion (50 kb).** Each point is one HPRC sample (n = 7; HG002 CiFi was not run with this control). Y > X for 7/7 — excluding chr13–22 p-arms and chrX/Y *strengthens* the arm-similarity ↔ Hi-C correlation (CHM13 0.66 → 0.80; HG02148 0.15 → 0.21; NA19036 0.27 → 0.49). The signal is not driven by acrocentric/nucleolar/PAR contacts. Source: `no_acrocentric/50000bp/<sample>_global_test.tsv`.

**c, O/E-normalised contact, within vs between sequence community.** Mean O/E (chromosome-size and marginal normalised) over all *trans* arm pairs. Within > between by 5.9× (HG02559) to 18.4× (HG002 CiFi); robust across Hi-C / CiFi / Pore-C platforms (8/8). Source: `community_based/50000bp/<sample>_oe_matrix.tsv`.

**d, Per-community reproducibility heatmap, 15 sequence communities × 11 datasets, 50 kb.** Cell colour = log₂(observed / null mean within-community contact) from 10 000 random-label permutations; * q < 0.05, ** q < 0.001 (Benjamini–Hochberg). Grey = community absent. Source: `community_based/50000bp/<sample>_community_bootstrap_tests.tsv` (+ RPE-1 equivalent).
