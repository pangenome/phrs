#!/usr/bin/env python3
"""Synteny (module ORDER) concordance between chromosome ends.

Modules shared between subtelomeric arm ends are orientation-consistent (see
arm_module_orientation.py). The remaining question is ORDER: when two arm ends share the
same partner arms, do those partners sit at the same rank along both arms?

Metric. For every ordered pair (X = query arm, P = partner arm) take the bp-weighted mean
midpoint of X's alignments to P, as a fraction of the 500 kb flank. For a pair of arms
(X1, X2) sharing >= 4 partners (each with >= 5 kb), compute the Spearman correlation across
shared partners between the two arms' partner midpoints. A pair of arms with the same module
order scores near +1; a rearranged pair scores near 0 or negative. Repertoire overlap
(Jaccard of partner sets) is reported alongside so order and content can be told apart.

Data: HPRCv2 all-vs-all alignments of the telomeric 500 kb flanks
(/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95).

Usage: python3 scripts/cladistics/arm_module_order.py [--min-kb 5]
Writes scripts/cladistics/arm_module_order.tsv and prints the same-community vs
different-community comparison.
"""
from __future__ import annotations
import argparse, collections, glob, gzip, itertools, json, re, statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PAF_DIR = Path("/moosefs/guarracino/HPRCv2/PHR_III/all-vs-all.1Mb.p95")
ASSIGN = ROOT / "data/hprcv2.1Mb.subtelo.arm-leiden-k15.assignments.tsv"
OUT = ROOT / "scripts/cladistics/arm_module_order.tsv"
NAMERE = re.compile(r"_(chr[0-9XY]+_[pq])arm")


def spearman(x, y):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-kb", type=float, default=5.0)
    a = ap.parse_args()
    minbp = a.min_kb * 1000
    # (query arm, partner arm) -> [sum(bp*mid_frac), sum(bp)]
    agg = collections.defaultdict(lambda: [0.0, 0.0])
    files = sorted(glob.glob(str(PAF_DIR / "hprcv2.1Mb.telo_500kb_trimmed.all-vs-CHM13_0_chr*.paf.gz")))
    for f in files:
        m = re.search(r"_(chr[0-9XY]+_[pq])arm\.paf\.gz$", f)
        if not m:
            continue
        partner = m.group(1)
        with gzip.open(f, "rt") as fh:
            for line in fh:
                c = line.split("\t")
                qm = NAMERE.search(c[0])
                if not qm:
                    continue
                x = qm.group(1)
                if x == partner:
                    continue
                qs, qe = int(c[2]), int(c[3])
                L = qe - qs
                if L <= 0:
                    continue
                mid = (qs + qe) / 2.0 / float(c[1])
                e = agg[(x, partner)]
                e[0] += L * mid
                e[1] += L
    # partner midpoint profiles per arm
    prof = collections.defaultdict(dict)
    for (x, p), (w, tot) in agg.items():
        if tot >= minbp:
            prof[x][p] = w / tot
    comm = {r["ChromArm"]: r["Community"] for r in
            __import__("csv").DictReader(open(ASSIGN), delimiter="\t")}
    rows = []
    for x1, x2 in itertools.combinations(sorted(prof), 2):
        if x1 not in comm or x2 not in comm:
            continue
        shared = sorted(set(prof[x1]) & set(prof[x2]))
        if len(shared) < 4:
            continue
        rho = spearman([prof[x1][p] for p in shared], [prof[x2][p] for p in shared])
        if rho != rho:
            continue
        jac = len(shared) / len(set(prof[x1]) | set(prof[x2]))
        same = comm[x1] == comm[x2]
        rows.append([x1, x2, "same" if same else "diff", len(shared), f"{jac:.3f}", f"{rho:.3f}"])
    with open(OUT, "w") as fh:
        fh.write("arm1\tarm2\tcommunity_relation\tn_shared_partners\trepertoire_jaccard\torder_rho\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    same = [float(r[5]) for r in rows if r[2] == "same"]
    diff = [float(r[5]) for r in rows if r[2] == "diff"]
    sj = [float(r[4]) for r in rows if r[2] == "same"]
    dj = [float(r[4]) for r in rows if r[2] == "diff"]
    print(f"arm pairs with >=4 shared partners: {len(rows)} (same-community {len(same)}, different {len(diff)})")
    print(f"order concordance (Spearman rho): same-community median {st.median(same):+.3f} "
          f"(n={len(same)}), different median {st.median(diff):+.3f} (n={len(diff)})")
    print(f"repertoire Jaccard:               same-community median {st.median(sj):.3f}, "
          f"different median {st.median(dj):.3f}")
    print(f"fraction of same-community pairs with rho>0.5: "
          f"{sum(1 for x in same if x>0.5)/len(same):.2f}; different: "
          f"{sum(1 for x in diff if x>0.5)/len(diff):.2f}")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
