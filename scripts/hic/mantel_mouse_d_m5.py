#!/usr/bin/env python3
"""D-M5: arm-level Mantel test for the mouse sequence-vs-contact correlation.

Replaces the per-PHR-pair Spearman (n = 344 non-independent pairs) with a Mantel test on
arm-level matrices, with row-and-column permutation (10,000 replicates), which respects the
non-independence the reviewer objected to.

Matrices:
  sequence - mouse arm-level Jaccard distance, /moosefs/.../mouse_T2T/similarity_1Mb/mouse.dist_matrix.tsv
  contact  - arm-level mean normalised Hi-C contact, aggregated from the per-PHR-pair
             tables data/mouse_<stage>_phr_<res>bp_seqlevel.tsv (columns arm_a, arm_b,
             hic_contact_norm)

Statistic: Mantel correlation between arm-pair sequence SIMILARITY (1 - distance) and mean
contact, restricted to arm pairs present in both matrices. Reported per stage and per
resolution, with a permutation p and a bootstrap 95% CI from the permutation distribution.

Usage: python3 scripts/hic/mantel_mouse_d_m5.py [--perms 10000]
"""
from __future__ import annotations
import argparse, csv, json, itertools
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DIST = Path("/moosefs/guarracino/HPRCv2/PHR_III/mouse_T2T/similarity_1Mb/mouse.dist_matrix.tsv")
STAGES = ["leptotene", "zygotene", "pachytene", "diplotene"]
OUT = ROOT / "scripts/hic/mantel_mouse_d_m5.json"


def read_dist():
    rows = [l.rstrip("\n").split("\t") for l in open(DIST)]
    arms = rows[0][1:]
    d = {}
    for r in rows[1:]:
        a = r[0]
        for b, v in zip(arms, r[1:]):
            if v not in ("", "NA", "nan"):
                d[(a, b)] = float(v)
    return arms, d


def read_contact(stage, res):
    f = ROOT / f"data/mouse_{stage}_phr_{res}bp_seqlevel.tsv"
    acc = {}
    if not f.exists():
        return None
    for r in csv.DictReader(open(f), delimiter="\t"):
        a, b = r["arm_a"], r["arm_b"]
        if a == b:
            continue
        v = r.get("hic_contact_norm") or r.get("hic_contact_raw")
        if not v:
            continue
        acc.setdefault((a, b) if a < b else (b, a), []).append(float(v))
    return {k: float(np.mean(v)) for k, v in acc.items()}


def mantel(x, y, perms, rng):
    x = np.asarray(x); y = np.asarray(y)
    rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
    obs = np.corrcoef(rx, ry)[0, 1]
    n = len(x)
    cnt = 0
    null = []
    for _ in range(perms):
        p = rng.permutation(n)
        r = np.corrcoef(rx[p], ry)[0, 1]
        null.append(r)
        cnt += r >= obs
    null = np.asarray(null)
    lo, hi = np.percentile(null, [2.5, 97.5])
    return obs, (cnt + 1) / (perms + 1), lo, hi, float(null.mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=10000)
    a = ap.parse_args()
    rng = np.random.default_rng(20260518)
    arms, dist = read_dist()
    results = {}
    for stage in STAGES:
        for res in (20000, 50000):
            con = read_contact(stage, res)
            if not con:
                continue
            keys = [k for k in con if k in dist and k[0] != k[1]]
            if len(keys) < 20:
                continue
            xs = [1.0 - dist[k] for k in keys]
            ys = [con[k] for k in keys]
            obs, p, lo, hi, mu = mantel(xs, ys, a.perms, rng)
            results[f"{stage}_{res//1000}kb"] = dict(
                rho=round(obs, 4), p_perm=p, null_mean=round(mu, 4),
                ci95=[round(lo, 4), round(hi, 4)], n_arm_pairs=len(keys))
            print(f"{stage:10s} {res//1000:>2d}kb  n={len(keys):3d} arm pairs  "
                  f"rho={obs:+.3f}  p={p:.4f}  null mean {mu:+.3f}  "
                  f"perm 95% [{lo:+.3f},{hi:+.3f}]")
    OUT.write_text(json.dumps(results, indent=2) + "\n")
    print("wrote", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
