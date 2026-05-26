#!/usr/bin/env python3
"""D-M5: Arm-level Mantel test on mouse PHR Jaccard distance vs Zuo 2021
meiotic Hi-C contact, per stage (leptotene, zygotene, pachytene, diplotene).

Replaces the peer-review-flagged per-arm-pair Spearman (n=344 non-independent
pairs) with the statistically appropriate Mantel test plus a permutation
p-value and an arm-block bootstrap 95% CI on rho.

Inputs (must be readable on moosefs):
  --arm-dist /moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/subtelo_1Mb/similarity/mouse.dist_matrix.tsv
  --hic-dir  /moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/community_analysis_1Mb/50000bp
  (--hic-dir holds zuo2021_<stage>_oe_matrix.tsv for each stage)

Output:
  per-stage rows printed and written to --out-tsv:
    stage  n_arms  mantel_rho  perm_p  ci_lo  ci_hi  per_pair_spearman

Mantel test:
  - Spearman correlation between the upper-triangle (off-diagonal) entries
    of the two arm-level distance matrices, aligned by shared arm labels.
  - Permutation p: 10,000 row+column permutations of the Hi-C matrix.

Bootstrap 95% CI:
  - Arm-block resampling: with replacement, draw n_arms armnames; build the
    induced sub-matrices; compute Mantel rho; repeat n_boot=10,000 times;
    take 2.5/97.5 percentiles of the bootstrap distribution.
  - Arm-block (not pair-block) is the correct unit because the
    non-independence among the n_arms*(n_arms-1)/2 pairs is induced by
    sharing arms across pairs.

Sign convention. The mouse Hi-C contact matrix is in contact-strength units
(higher = closer in 3D). The mouse PHR Jaccard *distance* matrix (the file
`mouse.dist_matrix.tsv` produced by the pipeline) is `1 - jaccard` (higher =
more diverged). Feeding the distance matrix directly to a Spearman
correlation against contact strength yields a NEGATIVE rho. The upstream
pipeline (`analyze_hic_communities.py:1788`) converts the distance to
similarity via `sim = 1 - dist` before running Mantel so that the reported
rho is positive (similar arms contact more). This script does the same when
`--invert-arm-dist` is set (the default, matching the upstream convention).
"""

from __future__ import annotations

import argparse
import itertools
import math
import os
import sys
from typing import Iterable, Sequence

import numpy as np
from scipy.stats import spearmanr


STAGES = ("leptotene", "zygotene", "pachytene", "diplotene")


def log(msg: str) -> None:
    print(f"[mantel_d_m5] {msg}", file=sys.stderr, flush=True)


def load_labeled_matrix(path: str) -> tuple[list[str], np.ndarray]:
    """Read a TSV whose first column is the row label and first row is the
    column header. Returns (labels, NxN matrix) with diagonal preserved as
    written."""
    log(f"load_labeled_matrix path={path}")
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col_labels = header[1:]
        rows = []
        row_labels = []
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            row_labels.append(parts[0])
            rows.append([float(x) if x not in ("", "NA", "nan") else float("nan")
                         for x in parts[1:]])
    if row_labels != col_labels:
        raise ValueError(f"row/col labels disagree in {path}: "
                         f"rows={row_labels[:3]}... cols={col_labels[:3]}...")
    arr = np.asarray(rows, dtype=float)
    if arr.shape[0] != arr.shape[1]:
        raise ValueError(f"non-square matrix in {path}: {arr.shape}")
    log(f"  matrix shape={arr.shape}, labels={len(row_labels)}")
    return row_labels, arr


def align_shared(labels_a: Sequence[str], mat_a: np.ndarray,
                 labels_b: Sequence[str], mat_b: np.ndarray
                 ) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Reduce both matrices to the same arm set (intersection), reordered
    consistently."""
    shared = [a for a in labels_a if a in labels_b]
    log(f"align_shared: |A|={len(labels_a)} |B|={len(labels_b)} shared={len(shared)}")
    idx_a = [labels_a.index(a) for a in shared]
    idx_b = [labels_b.index(a) for a in shared]
    sub_a = mat_a[np.ix_(idx_a, idx_a)]
    sub_b = mat_b[np.ix_(idx_b, idx_b)]
    return shared, sub_a, sub_b


def upper_offdiag(mat: np.ndarray) -> np.ndarray:
    """Upper-triangle (strict, off-diagonal) entries of an NxN matrix as a
    flat 1-D array."""
    n = mat.shape[0]
    iu, ju = np.triu_indices(n, k=1)
    return mat[iu, ju]


def mantel_rho(mat_a: np.ndarray, mat_b: np.ndarray) -> float:
    va = upper_offdiag(mat_a)
    vb = upper_offdiag(mat_b)
    # drop pairs where either side is NaN
    keep = ~(np.isnan(va) | np.isnan(vb))
    if keep.sum() < 4:
        return float("nan")
    rho, _ = spearmanr(va[keep], vb[keep])
    return float(rho)


def mantel_perm_p(mat_a: np.ndarray, mat_b: np.ndarray,
                  n_perm: int, seed: int) -> tuple[float, float, int]:
    """Mantel with row+column permutation of mat_b. Returns
    (observed_rho, two-sided p, n_pairs_compared)."""
    log(f"mantel_perm_p: n_perm={n_perm} seed={seed}")
    rng = np.random.default_rng(seed)
    n = mat_a.shape[0]
    if n != mat_b.shape[0]:
        raise ValueError("matrices must be aligned to the same arm order")
    rho_obs = mantel_rho(mat_a, mat_b)
    if not np.isfinite(rho_obs):
        return float("nan"), float("nan"), 0
    count_extreme = 0
    n_finite_pairs = int(np.sum(~(np.isnan(upper_offdiag(mat_a))
                                  | np.isnan(upper_offdiag(mat_b)))))
    for _ in range(n_perm):
        perm = rng.permutation(n)
        mat_b_perm = mat_b[np.ix_(perm, perm)]
        rho_perm = mantel_rho(mat_a, mat_b_perm)
        if not np.isfinite(rho_perm):
            continue
        # two-sided: |perm| >= |obs|
        if abs(rho_perm) >= abs(rho_obs):
            count_extreme += 1
    p_val = (count_extreme + 1) / (n_perm + 1)  # add-one convention
    return rho_obs, p_val, n_finite_pairs


def bootstrap_ci(mat_a: np.ndarray, mat_b: np.ndarray, n_boot: int,
                 seed: int, ci_level: float = 0.95
                 ) -> tuple[float, float, int, np.ndarray]:
    """Arm-block bootstrap: resample n arms with replacement, build induced
    submatrices, compute Mantel rho per resample. Returns (lo, hi, n_valid,
    distribution)."""
    log(f"bootstrap_ci: n_boot={n_boot} seed={seed} ci_level={ci_level}")
    rng = np.random.default_rng(seed)
    n = mat_a.shape[0]
    dist = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        # When the same arm is drawn twice, the induced submatrix has its
        # diagonal repeated. We exclude self-pairs (i==j in resample space)
        # by relying on upper-triangle indexing later — that means duplicate
        # rows still contribute pair entries among them. Use Manly-style
        # resample (all indices kept).
        sub_a = mat_a[np.ix_(idx, idx)]
        sub_b = mat_b[np.ix_(idx, idx)]
        dist[b] = mantel_rho(sub_a, sub_b)
    valid = dist[np.isfinite(dist)]
    alpha = (1 - ci_level) / 2
    lo = float(np.quantile(valid, alpha))
    hi = float(np.quantile(valid, 1 - alpha))
    return lo, hi, int(valid.size), dist


def per_pair_spearman(mat_a: np.ndarray, mat_b: np.ndarray
                      ) -> tuple[float, int]:
    """The naive per-arm-pair Spearman the v5 paper reports (treats the
    n*(n-1)/2 pairs as independent observations). Same number as Mantel
    rho but reported separately for the per-stage table."""
    va = upper_offdiag(mat_a)
    vb = upper_offdiag(mat_b)
    keep = ~(np.isnan(va) | np.isnan(vb))
    if keep.sum() < 4:
        return float("nan"), 0
    rho, _ = spearmanr(va[keep], vb[keep])
    return float(rho), int(keep.sum())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm-dist", required=True,
                    help="mouse arm-level Jaccard distance matrix TSV "
                         "(named `mouse.dist_matrix.tsv` in the pipeline)")
    ap.add_argument("--invert-arm-dist", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="convert distance to similarity via sim = 1 - dist "
                         "before Mantel (default: True; matches the upstream "
                         "pipeline convention so rho is positive)")
    ap.add_argument("--hic-dir", required=True,
                    help="directory holding zuo2021_<stage>_oe_matrix.tsv")
    ap.add_argument("--hic-pattern", default="zuo2021_{stage}_oe_matrix.tsv",
                    help="filename pattern with {stage} placeholder")
    ap.add_argument("--stages", nargs="+", default=list(STAGES))
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-tsv", required=True,
                    help="per-stage Mantel result TSV")
    ap.add_argument("--boot-dist-dir", default=None,
                    help="if set, write per-stage bootstrap distribution to "
                         "<dir>/mantel_d_m5_<stage>_boot.tsv")
    args = ap.parse_args()

    log(f"args={vars(args)}")
    arm_labels, arm_dist = load_labeled_matrix(args.arm_dist)
    log(f"arm distance matrix loaded: n_arms={len(arm_labels)}")
    if args.invert_arm_dist:
        log("converting arm distance -> similarity via sim = 1 - dist")
        arm_dist = 1.0 - arm_dist
        # the variable name keeps "arm_dist" for code reuse, but it now
        # holds the similarity matrix.

    out_rows = []
    for stage in args.stages:
        hic_path = os.path.join(args.hic_dir,
                                args.hic_pattern.format(stage=stage))
        log(f"--- stage={stage} hic_path={hic_path}")
        if not os.path.exists(hic_path):
            log(f"  MISSING: {hic_path} — skipping stage")
            out_rows.append((stage, "NA", "NA", "NA", "NA", "NA", "NA"))
            continue
        hic_labels, hic_mat = load_labeled_matrix(hic_path)

        shared, sub_dist, sub_hic = align_shared(arm_labels, arm_dist,
                                                 hic_labels, hic_mat)
        n_shared = len(shared)
        log(f"  n_shared_arms={n_shared}")

        rho_obs, perm_p, n_pairs = mantel_perm_p(
            sub_dist, sub_hic, n_perm=args.n_perm, seed=args.seed)
        log(f"  mantel_rho={rho_obs:.4f}, perm_p={perm_p:.6f}, "
            f"n_pairs={n_pairs}")

        ci_lo, ci_hi, n_valid, dist = bootstrap_ci(
            sub_dist, sub_hic, n_boot=args.n_boot, seed=args.seed + 1)
        log(f"  bootstrap CI 95%: ({ci_lo:.4f}, {ci_hi:.4f}) "
            f"n_valid={n_valid}/{args.n_boot}")

        per_rho, per_n = per_pair_spearman(sub_dist, sub_hic)
        log(f"  per-pair spearman rho={per_rho:.4f} (n={per_n} pairs)")

        out_rows.append((stage, str(n_shared), f"{rho_obs:.4f}",
                         f"{perm_p:.6f}",
                         f"({ci_lo:.4f}, {ci_hi:.4f})",
                         f"{per_rho:.4f}",
                         str(per_n)))

        if args.boot_dist_dir:
            os.makedirs(args.boot_dist_dir, exist_ok=True)
            boot_path = os.path.join(args.boot_dist_dir,
                                     f"mantel_d_m5_{stage}_boot.tsv")
            with open(boot_path, "w") as fh:
                fh.write("rep\trho\n")
                for i, r in enumerate(dist):
                    fh.write(f"{i}\t{r}\n")
            log(f"  wrote bootstrap distribution -> {boot_path}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)) or ".",
                exist_ok=True)
    header = ["stage", "n_arms", "mantel_rho", "perm_p", "ci95",
              "per_pair_spearman", "per_pair_n"]
    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in out_rows:
            fh.write("\t".join(row) + "\n")
    log(f"wrote {args.out_tsv} ({len(out_rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
