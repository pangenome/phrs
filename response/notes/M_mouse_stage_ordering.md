# M — Mouse stage-profile disambiguation (statistic vs axis vs satellite)

Status: IN PROGRESS. Notes appended after every step.

## Skeleton / plan
- [ ] STEP A — axis audit (odgi graph-node Jaccard vs re-derived PAF Jaccard)
- [ ] STEP B — same axis, two statistics (per-PHR-pair Spearman vs arm-level Mantel)
- [ ] STEP C — referee's test R1.7 (stratified pair-label permutation + arm bootstrap CI)
- [ ] STEP D — satellite control at pair resolution
- [ ] STEP E — verdict, manuscript sentence, what is unmeasurable

Outputs: `scripts/hic/mouse_stage_profile_disambiguation.py`, `scripts/hic/mouse_stage_profile_disambiguation.json`

## Log

## Log

### L0 — inputs read (before any computation)
- `PEER_REVIEW.md` (authority; R1 major points 6 and 7 read verbatim).
  `paper_prep/synthesis/SIMULATED_PEER_REVIEW_v1.md` **never opened**.
- `response/notes/M6_mouse_satellite.md` (authority for the satellite lane).
- `scripts/hic/mantel_mouse_d_m5.py` (105 lines) + `scripts/hic/mantel_mouse_d_m5.json`.
- `scripts/hic/mouse_satellite_sensitivity.py` (1038 lines) + its JSON.
- `submission/scripts/figures/make_fig4c_mouse_zygotene.R` (the submitted analysis).
- `data/mouse_{leptotene,zygotene,pachytene,diplotene}_phr_{20000,50000}bp_seqlevel.tsv`
  (1135 rows at 20 kb; 1088 / 1135 / 1135 / 999 at 50 kb).
- `data/mouse_meiosis_sweep/seqlevel/{2Mb,4Mb}/*_seqlevel_{20000,50000}bp_seqlevel.tsv`
  (robustness windows).
- `/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/similarity_1Mb/mouse.dist_matrix.tsv`
  (27 arms) and 4 PAF files under `subtelo_1Mb/all-vs-all.p95/`.

### L1 — first-pass audit findings (before writing the script)

Two of these change the framing of the whole question.

**(a) `mantel_mouse_d_m5.py` does NOT use a re-derived PAF Jaccard.** Line 26:
`DIST = Path(".../similarity_1Mb/mouse.dist_matrix.tsv")`. That is the *published odgi
graph-node distance matrix* — the same quantity analysis (1) uses. The "re-derived
PAF-based Jaccard" exists only inside `mouse_satellite_sensitivity.py` and is **void**
(M6 §5b). So at the arm level there is no competing axis at all.

**(b) The published arm-level axis IS the arm-pair mean of the pair-level `jaccard`
column.** Max |difference| over the 344 arm pairs = **1.1e-16** (float noise). So
analysis (1) and analysis (2) use *literally the same numbers*, aggregated differently.
There is one axis, not two.

**(c) The published Mantel routine ranks with `np.argsort(np.argsort(x))`,** which breaks
ties by row order. 129 of the 344 arm pairs are tied at similarity 0, so the published
rho is not the Spearman coefficient; it is row-order dependent. Tie-corrected (average
ranks) the eight published values move by up to 0.049 with the same stage order.

**(d) The satellite-poor subset is degenerate in the sequence axis.** At 20 kb the
42 both-satellite-poor (q×q) PHR pairs have only **8 distinct Jaccard values, and 35 of
42 are exactly 0**. The 21 arm-level q×q pairs M6 reports (17 of them tied at
similarity 0) have the same problem. M6's "low-satellite flat 0.536/0.539/0.486/0.486"
is therefore a rank correlation resting on 4 non-tied arm pairs.

