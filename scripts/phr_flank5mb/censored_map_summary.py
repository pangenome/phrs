#!/usr/bin/env python3
"""Per 100 kb from the telomere of each query (censored_map.sbatch output): aligned bp per CHM13 chromosome (query
bases covered by records to that chromosome, PAF columns 1-4 and 6), the top chromosomes, and the CHM13 position of the
longest record. A query p arm starts at its telomere; a q arm ends at it."""
import argparse
from collections import defaultdict

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--paf", required=True)
ap.add_argument("--fai", required=True)
ap.add_argument("--bin", type=int, default=100000)
a = ap.parse_args()
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai)}
cov = defaultdict(lambda: defaultdict(int))
best = {}
for l in open(a.paf):
    x = l.split("\t")
    q, qs, qe, t, ts = x[0], int(x[2]), int(x[3]), x[5].split("#")[-1], int(x[7])
    p = q.endswith("_parm")
    for b0 in range(qs // a.bin * a.bin, qe, a.bin):
        ov = min(qe, b0 + a.bin) - max(qs, b0)
        d = b0 if p else L[q] - b0 - a.bin  # bin start as distance from the telomere
        cov[(q, d)][t] += ov
        if ov > best.get((q, d), (0, ""))[0]:
            best[(q, d)] = (ov, f"{t}:{ts / 1e6:.2f}Mb")
for q in L:
    print(f"== {q} ({L[q]} bp)")
    for d in range(0, L[q] - a.bin + 1, a.bin):
        c = cov.get((q, d), {})
        top = ", ".join(f"{t} {v // 1000} kb" for t, v in sorted(c.items(), key=lambda kv: -kv[1])[:4]) or "unaligned"
        print(f"   {d / 1e6:.1f}-{(d + a.bin) / 1e6:.1f} Mb: {top}; longest record at {best.get((q, d), (0, '.'))[1]}")
