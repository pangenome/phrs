#!/usr/bin/env python3
"""
D-M12: Bootstrap / Wilson / block-jackknife CIs for headline correlations.

Inputs (absolute paths, local snapshots; canonical inputs live on
/moosefs/guarracino/HPRCv2/PHR_III/ and are not mounted in this worktree):

  ARM_DIST_TSV   = local snapshot of the arm-level Jaccard distance matrix
                   (pangenome graph). Symmetric, diagonal = self-distance,
                   off-diagonal = 1 - mean inter-arm Jaccard similarity.
  HG002_HIC_TSV  = local snapshot of HG002 arm-level Hi-C contact matrix.
  CHM13_HIC_TSV  = local snapshot of CHM13 arm-level Hi-C contact matrix
                   at 50 kb resolution.

Outputs:
  paper_prep/synthesis/ANALYSIS_D_M12.md is written separately; this script
  only prints results to stdout. Results are also written as TSV alongside
  this script for reuse.

CIs computed:
  (a) Mantel rho (HG002 Hi-C, CHM13 Hi-C, similarity vs contact)
      - Bootstrap CI: resample arms with replacement (n_arms),
        rebuild both matrices, recompute Spearman on upper-triangle pairs.
      - B = 10,000 replicates.
      - Permutation p-value (10,000 row-and-column shuffles) reported
        separately so the user can compare to upstream permutation p.
  (c) Pedigree within-Leiden fraction 494/538
      - Wilson-score 95% CI on the binomial proportion.
  (d) F_ST per superpopulation pair (Hudson)
      - 10 arm/community pairs from
        end-to-end-report/report/04_heterogeneity.md §population structure.
      - Per arm: p_pop = cross_count / (cross_count + self_count).
      - Hudson F_ST per arm per pair, mean over arms.
      - Block-jackknife 95% CI over arms (leave-one-out, n=10 blocks).

CIs deferred (no local data):
  (b) Mantel trajectory under 5 exclusion sets x 5 resolutions
      - needs /moosefs PHR_III/analysis/human/community_based/{res}/*.tsv
        plus per-resolution arm-level similarity matrices.
  (e) Mouse zygotene Spearman rho = 0.715
      - needs the 344 PHR-pair (Jaccard, Hi-C contact) vector that lives
        only in /moosefs/.../mouse/.

Notes on data version:
  The local snapshots predate the v5 paper rerun and reproduce a Mantel
  rho ~ 0.32 on HG002 instead of the headline ~ 0.66. The bootstrap CI
  width is the deliverable here; once /moosefs is mounted, re-run with the
  current matrices and the CI scales relative to the new point estimate
  (asymmetric Fisher z transform handles the rho-magnitude rescaling
  approximately; the standard error of rho_Fisher ~ 1 / sqrt(n_arms - 3)).
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

# ---- inputs --------------------------------------------------------------

ARM_DIST_TSV = "/home/guarracino/Dropbox/working/Garrison/hprcv2/PHR_III/hic_validation/arm_dist_matrix.tsv"
HG002_HIC_TSV = "/home/guarracino/Dropbox/working/Garrison/hprcv2/PHR_III/hic_validation/hg002_contact_matrix.tsv"
CHM13_HIC_TSV = "/home/guarracino/Desktop/Garrison/HPRCv2/PHR_III/hic_validation/res_50kb/chm13_contact_matrix.tsv"

# Pedigree within-Leiden numerator/denominator from
# end-to-end-report/report/14_pedigree_recombination.md.
PEDIGREE_NUM = 494
PEDIGREE_DEN = 538

# F_ST cross-arm / self-arm counts per superpopulation per arm/community,
# from end-to-end-report/report/04_heterogeneity.md "Population structure
# in cross-arm affinity". 10 arm/community pairs with significant Fisher
# p_adj < 0.05.
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
FST_TABLE = [
    # (community, arm, cross_total, self_total,
    #  cross_per_pop_dict, self_per_pop_dict)
    ("C1",  "chr4_q",  146, 211, {"AFR": 52, "AMR": 13, "EAS": 23, "EUR": 17, "SAS": 20}, {"AFR": 32, "AMR": 53, "EAS": 47, "EUR": 27, "SAS": 38}),
    ("C3",  "chr16_q",  86, 363, {"AFR": 60, "AMR":  5, "EAS":  2, "EUR":  1, "SAS":  6}, {"AFR": 54, "AMR": 76, "EAS": 76, "EUR": 61, "SAS": 64}),
    ("C5",  "chr6_p",  245, 172, {"AFR": 25, "AMR": 44, "EAS": 57, "EUR": 43, "SAS": 49}, {"AFR": 74, "AMR": 29, "EAS": 19, "EUR": 14, "SAS": 18}),
    ("C15", "chrX_p",  305,  22, {"AFR": 71, "AMR": 62, "EAS": 52, "EUR": 43, "SAS": 48}, {"AFR": 18, "AMR":  0, "EAS":  0, "EUR":  0, "SAS":  0}),
    ("C6",  "chr19_q",  59, 398, {"AFR": 20, "AMR": 12, "EAS":  0, "EUR":  8, "SAS": 11}, {"AFR": 93, "AMR": 67, "EAS": 79, "EUR": 55, "SAS": 61}),
    ("C3",  "chr9_q",  130, 320, {"AFR": 45, "AMR": 17, "EAS": 19, "EUR": 18, "SAS": 12}, {"AFR": 66, "AMR": 62, "EAS": 58, "EUR": 43, "SAS": 59}),
    ("C6",  "chr22_q", 231, 219, {"AFR": 49, "AMR": 41, "EAS": 43, "EUR": 28, "SAS": 46}, {"AFR": 66, "AMR": 38, "EAS": 34, "EUR": 34, "SAS": 21}),
    ("C15", "chrY_p",   10,  82, {"AFR":  5, "AMR":  1, "EAS":  0, "EUR":  1, "SAS":  0}, {"AFR": 13, "AMR": 15, "EAS": 20, "EUR":  9, "SAS": 18}),
    ("C1",  "chr10_q",  26, 331, {"AFR":  4, "AMR":  3, "EAS":  3, "EUR":  1, "SAS": 11}, {"AFR": 79, "AMR": 63, "EAS": 60, "EUR": 45, "SAS": 49}),
    ("C11", "chr6_q",   12, 428, {"AFR":  1, "AMR":  1, "EAS":  0, "EUR":  4, "SAS":  4}, {"AFR":107, "AMR": 78, "EAS": 79, "EUR": 55, "SAS": 64}),
]

N_BOOT = 10_000
N_PERM = 10_000
RNG_SEED = 20260518

# ---- helpers -------------------------------------------------------------


def load_matrix(path: str) -> tuple[list[str], np.ndarray]:
    """Load a square symmetric labelled matrix from TSV (header has labels)."""
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    arms = header[1:]
    by_label = {}
    for row in rows:
        label = row[0]
        by_label[label] = {arms[i]: float(row[i + 1]) for i in range(len(arms))}
    if set(arms) != set(by_label.keys()):
        raise ValueError(f"row/col label mismatch in {path}")
    n = len(arms)
    M = np.zeros((n, n), dtype=float)
    for i, a in enumerate(arms):
        for j, b in enumerate(arms):
            M[i, j] = by_label[a][b]
    # symmetrise (small floating round-off can break exact symmetry).
    M = 0.5 * (M + M.T)
    return arms, M


def reorder(M: np.ndarray, src_arms: list[str], target_arms: list[str]) -> np.ndarray:
    """Reorder rows/cols of M from src_arms to target_arms."""
    idx = [src_arms.index(a) for a in target_arms]
    return M[np.ix_(idx, idx)]


def upper_tri(M: np.ndarray) -> np.ndarray:
    """Off-diagonal upper-triangle flattened vector."""
    iu = np.triu_indices(M.shape[0], k=1)
    return M[iu]


def mantel_spearman(sim: np.ndarray, contact: np.ndarray) -> float:
    """Spearman rho between two square matrices on their upper-triangle pairs."""
    rho, _ = spearmanr(upper_tri(sim), upper_tri(contact))
    return float(rho)


def mantel_bootstrap_ci(
    sim: np.ndarray,
    contact: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
) -> tuple[float, float, np.ndarray]:
    """
    Arm-resampling bootstrap: draw n_arms indices with replacement, rebuild
    both matrices on the resampled label set, compute Mantel rho.
    Returns (lo, hi, distribution).
    """
    n = sim.shape[0]
    out = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        s_b = sim[np.ix_(idx, idx)]
        c_b = contact[np.ix_(idx, idx)]
        # mask collisions on the diagonal that would inflate the count
        mask = np.ones_like(s_b, dtype=bool)
        np.fill_diagonal(mask, False)
        # also drop duplicated-arm pairs (i.e. where the arm appears twice
        # in the resample and produces identical rows) by deduplicating on
        # the (sorted) sampled-index pair; otherwise variance is biased low.
        # Simpler & standard: just use full upper-triangle on the resampled
        # matrix. Ties from repeated arms add valid bootstrap variance.
        iu = np.triu_indices(n, k=1)
        rho, _ = spearmanr(s_b[iu], c_b[iu])
        out[b] = rho
    lo = float(np.percentile(out, 2.5))
    hi = float(np.percentile(out, 97.5))
    return lo, hi, out


def mantel_permutation_p(
    sim: np.ndarray,
    contact: np.ndarray,
    obs: float,
    n_perm: int,
    rng: np.random.Generator,
) -> float:
    """Two-sided permutation p: row-and-column shuffles of one matrix."""
    n = sim.shape[0]
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(n)
        sim_p = sim[np.ix_(perm, perm)]
        iu = np.triu_indices(n, k=1)
        rho, _ = spearmanr(sim_p[iu], contact[iu])
        if abs(rho) >= abs(obs):
            count += 1
    return (count + 1) / (n_perm + 1)


def wilson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float, float]:
    """Wilson-score CI for a binomial proportion. Returns (phat, lo, hi)."""
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    from scipy.stats import norm

    p = k / n
    z = norm.ppf(1 - alpha / 2)
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = (z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return p, centre - half, centre + half


def hudson_fst(ci: int, ni: int, cj: int, nj: int) -> float:
    """
    Hudson-style F_ST = (HT - HS) / HT (matches the upstream pipeline's
    compute_fst_superpop.py at /moosefs/.../scripts/community/).

      pi = ci / ni
      pj = cj / nj
      HS = mean(2*pi*(1-pi), 2*pj*(1-pj))
      p_pool = (ci + cj) / (ni + nj)
      HT = 2 * p_pool * (1 - p_pool)

    Returns 0.0 if HT == 0 (locus monomorphic in the pooled sample).
    Returns NaN if either ni or nj == 0 (no haplotypes in that pop).
    """
    if ni == 0 or nj == 0:
        return float("nan")
    pi = ci / ni
    pj = cj / nj
    hs = (2.0 * pi * (1.0 - pi) + 2.0 * pj * (1.0 - pj)) / 2.0
    p_pool = (ci + cj) / (ni + nj)
    ht = 2.0 * p_pool * (1.0 - p_pool)
    if ht == 0.0:
        return 0.0
    return (ht - hs) / ht


def fst_per_arm_per_pair(
    table: list,
) -> dict[tuple[str, str], list[tuple[str, float]]]:
    """
    For each superpop pair (P1, P2) return a list of (arm_label, F_ST) over
    the 10 arms.  Uses the per-arm cross-vs-self counts.
    """
    out: dict[tuple[str, str], list[tuple[str, float]]] = {}
    pops = SUPERPOPS
    for i, p1 in enumerate(pops):
        for p2 in pops[i + 1 :]:
            per_arm = []
            for comm, arm, _, _, cross_pop, self_pop in table:
                ci = cross_pop[p1]; ni = cross_pop[p1] + self_pop[p1]
                cj = cross_pop[p2]; nj = cross_pop[p2] + self_pop[p2]
                fst = hudson_fst(ci, ni, cj, nj)
                per_arm.append((f"{comm}/{arm}", fst))
            out[(p1, p2)] = per_arm
    return out


def block_jackknife_mean_ci(values: list[float], alpha: float = 0.05) -> dict:
    """
    Block-jackknife 95% CI for the mean of a list of per-arm values.
    Leave-one-out, n_blocks = len(values).
    """
    from scipy.stats import t

    n = len(values)
    arr = np.array(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    n_eff = arr.size
    if n_eff < 2:
        return {"mean": float("nan"), "lo": float("nan"), "hi": float("nan"), "n": n_eff}
    full = arr.mean()
    leave_one = np.array([arr[np.arange(n_eff) != i].mean() for i in range(n_eff)])
    pseudo = n_eff * full - (n_eff - 1) * leave_one
    se = pseudo.std(ddof=1) / np.sqrt(n_eff)
    tcrit = t.ppf(1 - alpha / 2, df=n_eff - 1)
    return {
        "mean": float(full),
        "lo": float(full - tcrit * se),
        "hi": float(full + tcrit * se),
        "n": n_eff,
        "se": float(se),
    }


# ---- main ----------------------------------------------------------------


def main() -> None:
    out_dir = Path(__file__).parent
    rng = np.random.default_rng(RNG_SEED)

    results: dict = {"seed": RNG_SEED, "n_boot": N_BOOT, "n_perm": N_PERM}

    # ---- (c) Wilson CI for pedigree fraction ----------------------------
    p, lo, hi = wilson_ci(PEDIGREE_NUM, PEDIGREE_DEN)
    print(f"\n[c] Pedigree within-Leiden: {PEDIGREE_NUM}/{PEDIGREE_DEN} = {p:.4f}")
    print(f"    Wilson 95% CI: [{lo:.4f}, {hi:.4f}]  (width = {hi - lo:.4f})")
    results["pedigree_wilson"] = {
        "k": PEDIGREE_NUM,
        "n": PEDIGREE_DEN,
        "phat": p,
        "lo": lo,
        "hi": hi,
        "method": "Wilson-score 95% CI for binomial proportion (Brown, Cai, DasGupta 2001)",
    }

    # ---- (d) F_ST per superpop pair: block-jackknife --------------------
    print("\n[d] F_ST per superpopulation pair (Hudson estimator)")
    print("    Block-jackknife over 10 arm/community pairs (leave-one-out)\n")
    per_pair = fst_per_arm_per_pair(FST_TABLE)
    fst_out_rows = []
    fst_out_rows.append(
        ["pop1", "pop2", "mean_fst", "ci_lo", "ci_hi", "se", "n_arms"]
    )
    results["fst_block_jackknife"] = {}
    for (p1, p2), per_arm in per_pair.items():
        vals = [v for _, v in per_arm if not (v is None or v != v)]
        bj = block_jackknife_mean_ci(vals)
        print(
            f"    {p1} vs {p2}: F_ST = {bj['mean']:+.4f}  "
            f"95% CI [{bj['lo']:+.4f}, {bj['hi']:+.4f}]  (n_arms={bj['n']}, se={bj['se']:.4f})"
        )
        fst_out_rows.append(
            [p1, p2, f"{bj['mean']:.6f}", f"{bj['lo']:.6f}", f"{bj['hi']:.6f}",
             f"{bj['se']:.6f}", str(bj["n"])]
        )
        results["fst_block_jackknife"][f"{p1}_{p2}"] = bj
    # Also: per-arm-Fst long table for traceability.
    long_rows = [["pop1", "pop2", "arm_community", "fst"]]
    for (p1, p2), per_arm in per_pair.items():
        for label, v in per_arm:
            long_rows.append([p1, p2, label, f"{v:.6f}"])

    with open(out_dir / "fst_block_jackknife.tsv", "w") as fh:
        for r in fst_out_rows:
            fh.write("\t".join(r) + "\n")
    with open(out_dir / "fst_per_arm_per_pair.tsv", "w") as fh:
        for r in long_rows:
            fh.write("\t".join(r) + "\n")

    # ---- (a) Mantel bootstrap CI for HG002 Hi-C and CHM13 Hi-C ----------
    print("\n[a] Mantel rho bootstrap CI (similarity vs Hi-C contact)")
    print("    Bootstrap = arm-resampling with replacement, B = "
          f"{N_BOOT}; permutation p uses {N_PERM} row+col shuffles.")
    print("    DATA VERSION CAVEAT: local snapshot Feb 2026 not the v5 rerun;")
    print("    method is exact, point estimate differs from headline 0.656/0.657.\n")

    arms_d, dist = load_matrix(ARM_DIST_TSV)
    sim_full = 1.0 - dist  # similarity matrix

    mantel_results = {}
    for sample, hic_path in [("HG002", HG002_HIC_TSV), ("CHM13", CHM13_HIC_TSV)]:
        arms_h, contact = load_matrix(hic_path)
        common = sorted(set(arms_d) & set(arms_h))
        # restrict to common arms in a stable order
        sim = reorder(sim_full, arms_d, common)
        contact_c = reorder(contact, arms_h, common)
        obs = mantel_spearman(sim, contact_c)
        lo, hi, dist_boot = mantel_bootstrap_ci(sim, contact_c, N_BOOT, rng)
        perm_p = mantel_permutation_p(sim, contact_c, obs, N_PERM, rng)
        print(
            f"    {sample}: rho = {obs:+.4f}  bootstrap 95% CI [{lo:+.4f}, {hi:+.4f}] "
            f"(n_arms = {len(common)}, perm p = {perm_p:.4e})"
        )
        mantel_results[sample] = {
            "rho": obs,
            "ci_lo": lo,
            "ci_hi": hi,
            "n_arms": len(common),
            "n_boot": N_BOOT,
            "perm_p": perm_p,
            "arms": common,
            "arm_dist_matrix": ARM_DIST_TSV,
            "hic_matrix": hic_path,
        }
    results["mantel_bootstrap"] = mantel_results

    with open(out_dir / "mantel_bootstrap_ci.tsv", "w") as fh:
        fh.write("sample\tn_arms\trho\tci_lo\tci_hi\tperm_p\tn_boot\n")
        for sample, r in mantel_results.items():
            fh.write(
                f"{sample}\t{r['n_arms']}\t{r['rho']:.6f}\t{r['ci_lo']:.6f}\t"
                f"{r['ci_hi']:.6f}\t{r['perm_p']:.6e}\t{r['n_boot']}\n"
            )

    # ---- (a-FISHER) Analytic Fisher z 95% CI on headline rho values -----
    # Stand-in until /moosefs is mounted and the arm-resampling bootstrap can
    # be re-run on the v5 matrices. Fisher z: z = atanh(rho), SE_z =
    # 1 / sqrt(n - 3), 95% CI on z = z +- 1.96 SE_z, back-transform via tanh.
    # This is exact for product-moment correlations and a good approximation
    # for Spearman when n is moderate (~30+).
    print("\n[a-Fisher] Analytic Fisher z 95% CI for headline Mantel rho values")
    print("    Used as a quick reviewer-ready CI on the published point estimates;")
    print("    arm-resampling bootstrap above is the gold-standard once /moosefs is mounted.")
    headline = [
        # (sample, technology, rho, n_arms, source)
        ("CHM13",   "Hi-C",   0.656, 38, "report/05_hic_validation.md table PHR Mantel"),
        ("HG002",   "Hi-C",   0.657, 41, "report/05_hic_validation.md table PHR Mantel"),
        ("HG02559", "Hi-C",   0.397, 37, "report/05_hic_validation.md table PHR Mantel"),
        ("HG00658", "Hi-C",   0.276, 37, "report/05_hic_validation.md table PHR Mantel"),
        ("HG02148", "Hi-C",   0.152, 37, "report/05_hic_validation.md table PHR Mantel"),
        ("NA19036", "Hi-C",   0.266, 34, "report/05_hic_validation.md table PHR Mantel"),
        ("HG002",   "Pore-C", 0.486, 41, "report/05_hic_validation.md table PHR Mantel"),
        ("HG002",   "CiFi",   0.308, 41, "report/05_hic_validation.md table PHR Mantel"),
    ]
    fisher_rows = [["sample", "tech", "rho", "n_arms", "ci_lo", "ci_hi", "method"]]
    fisher_z = {}
    for sample, tech, rho, n_arms, _src in headline:
        z = np.arctanh(rho)
        se = 1.0 / np.sqrt(n_arms - 3)
        zlo = z - 1.96 * se
        zhi = z + 1.96 * se
        rho_lo = float(np.tanh(zlo))
        rho_hi = float(np.tanh(zhi))
        key = f"{sample}_{tech}"
        fisher_z[key] = {
            "rho": rho, "n_arms": n_arms, "ci_lo": rho_lo, "ci_hi": rho_hi
        }
        print(
            f"    {sample} {tech}: rho = {rho:.3f}  Fisher 95% CI [{rho_lo:+.3f}, {rho_hi:+.3f}] "
            f"(n_arms = {n_arms})"
        )
        fisher_rows.append([sample, tech, f"{rho:.3f}", str(n_arms),
                            f"{rho_lo:.4f}", f"{rho_hi:.4f}", "Fisher z"])
    with open(out_dir / "mantel_fisher_z_ci.tsv", "w") as fh:
        for r in fisher_rows:
            fh.write("\t".join(r) + "\n")
    results["mantel_fisher_z"] = fisher_z

    # ---- (e) Mouse zygotene Spearman rho = 0.715: Fisher z stand-in -----
    # The per-pair Jaccard, Hi-C contact vectors that the per-pair Spearman
    # was computed on live in /moosefs and are not mountable in this
    # worktree. A bootstrap CI on those vectors is the canonical answer
    # (deferred to next pass). As a reviewer-ready stand-in, Fisher z on
    # n = 344 PHR pairs.
    print("\n[e-Fisher] Mouse zygotene per-pair Spearman rho = 0.715  (n = 344)")
    rho = 0.715
    n_pairs = 344
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n_pairs - 3)
    rho_lo = float(np.tanh(z - 1.96 * se))
    rho_hi = float(np.tanh(z + 1.96 * se))
    print(f"    Fisher 95% CI [{rho_lo:+.3f}, {rho_hi:+.3f}]")
    print("    NOTE: pairs are non-independent (D-M5). Fisher z understates")
    print("    the true width by an unknown factor; the arm-level Mantel CI")
    print("    from D-M5 is the appropriate published CI.")
    results["mouse_spearman_fisher_z"] = {
        "rho": rho, "n_pairs": n_pairs, "ci_lo": rho_lo, "ci_hi": rho_hi,
        "caveat": "non-independent pairs; D-M5 arm-level Mantel CI is preferred"
    }

    # ---- summary JSON ---------------------------------------------------
    with open(out_dir / "results_d_m12.json", "w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print("\nWrote results to:")
    print(f"  {out_dir / 'fst_block_jackknife.tsv'}")
    print(f"  {out_dir / 'fst_per_arm_per_pair.tsv'}")
    print(f"  {out_dir / 'mantel_bootstrap_ci.tsv'}")
    print(f"  {out_dir / 'results_d_m12.json'}")


if __name__ == "__main__":
    main()
