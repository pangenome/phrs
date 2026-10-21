#!/usr/bin/env python3
"""Offline generator for the D-M6 matched F_ST comparison table.

`matched_fst_d_m6.py` recomputes per-arm F_ST from an upstream enrichment table that lives on
moosefs; that mount intermittently fails, which makes the manuscript's F_ST numbers hard to
reproduce. This script produces the same comparison table from the repository's own per-block
F_ST table (`scripts/ci/fst_per_arm_per_pair.tsv`) plus the published genome-wide continental
baselines, reusing the constants and the block-jackknife from `matched_fst_d_m6.py` so there is
only one implementation of either.

Output: scripts/popgen/matched_fst_d_m6.tsv (consumed by scripts/ci/audit_response_numbers.py)
"""
from __future__ import annotations
import collections, csv, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/popgen"))
from matched_fst_d_m6 import (  # noqa: E402  (single source of truth)
    PUBLISHED_FST_1000G_HGDP,
    block_jackknife_ci,
)

BLOCKS = ROOT / "scripts/ci/fst_per_arm_per_pair.tsv"
OUT = ROOT / "scripts/popgen/matched_fst_d_m6.tsv"


def main():
    by = collections.defaultdict(list)
    with open(BLOCKS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by[(r["pop1"], r["pop2"])].append(float(r["fst"]))
    rows = []
    for (p1, p2), vals in sorted(by.items()):
        mean, lo, hi = block_jackknife_ci(vals)
        key = (p1, p2) if (p1, p2) in PUBLISHED_FST_1000G_HGDP else (p2, p1)
        matched = PUBLISHED_FST_1000G_HGDP.get(key, float("nan"))
        if matched != matched:
            verdict = "no_published_baseline"
        elif lo <= matched <= hi:
            verdict = "equivalent"
        elif matched > hi:
            verdict = "depressed"
        else:
            verdict = "elevated"
        rows.append([f"{p1}-{p2}", f"{mean:+.4f}", f"{lo:+.4f}", f"{hi:+.4f}",
                     f"{matched:.3f}", f"{mean - matched:+.4f}",
                     f"{mean / matched:.3f}" if matched else "nan",
                     len(vals), verdict])
    with open(OUT, "w") as fh:
        fh.write("superpop_pair\tsubtelo_fst\tjk95_lo\tjk95_hi\tmatched_genomewide_fst\t"
                 "delta_subtelo_minus_matched\tratio_subtelo_over_matched\tn_blocks_used\tverdict\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)} ({len(rows)} superpopulation pairs)")
    for r in rows:
        if r[0].startswith("AFR"):
            print(f"  {r[0]:9s} subtelo {r[1]} [{r[2]},{r[3]}]  matched {r[4]}  "
                  f"ratio {r[6]}  -> {r[8]}")


if __name__ == "__main__":
    main()
