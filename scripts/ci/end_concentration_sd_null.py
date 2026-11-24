#!/usr/bin/env python3
"""Are high-identity inter-chromosomal windows concentrated at chromosome ends beyond what segmental
duplication density predicts? (R1-m6)

Windows: the 100 kb CHM13 windows of the genome-wide survey (Fig. 1A input, IMPG queries with up to five
hops at >= 98% identity). A window is SHARED when at least one chromosome other than its own (chrUn
excluded) is reached by >= --min-haps haplotypes. A window is at an END when it lies within --end-mb of
either end of its chromosome. SD density: the fraction of the window covered by segmental duplications
of the CHM13v2.0 SD annotation whose partner is on another chromosome (inter-chromosomal SDs), and,
separately, by any SD.
Statistic: the share of shared windows that are at an end. Null: shared labels permuted among windows
with the same SD density class (--classes edges), the number of shared windows per class kept fixed;
one-sided p = (1 + permutations with a share >= observed) / (1 + permutations).
"""
import argparse, bisect, gzip, random, re, sys
from collections import defaultdict

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--windows", default="data/hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz")
ap.add_argument("--sd", required=True, help="chm13v2.0_SD.bed (partner location in column 4)")
ap.add_argument("--min-haps", type=int, default=47, help="10%% of 466 haplotypes")
ap.add_argument("--end-mb", type=float, default=1.0)
ap.add_argument("--classes", default="0,0.001,0.1,0.25,0.5,0.75,1.01")
ap.add_argument("--perm", type=int, default=10000)
ap.add_argument("--seed", type=int, default=20261124)
a = ap.parse_args()
sd_all, sd_inter = defaultdict(list), defaultdict(list)
for l in open(a.sd):
    f = l.split("\t")
    sd_all[f[0]].append((int(f[1]), int(f[2])))
    if f[3].split(":")[0] != f[0]:
        sd_inter[f[0]].append((int(f[1]), int(f[2])))


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def cover(m, s, e):
    """bp of merged intervals m within [s, e)"""
    i = bisect.bisect_left(m, [s, -1]) - 1
    n = 0
    for x, y in m[max(i, 0):]:
        if x >= e:
            break
        n += max(0, min(e, y) - max(s, x))
    return n


ma = {c: merge(v) for c, v in sd_all.items()}
mi = {c: merge(v) for c, v in sd_inter.items()}
rows = [l.rstrip("\n").split("\t") for l in gzip.open(a.windows, "rt") if not l.startswith("chrom")]
size = defaultdict(int)  # chromosome length = end of its last window (the windows tile each chromosome)
for f in rows:
    size[f[0].split("#")[-1]] = max(size[f[0].split("#")[-1]], int(f[2]))
W = []
for f in rows:
    c = f[0].split("#")[-1]
    if c == "chrM":
        continue
    s, e = int(f[1]), int(f[2])
    parts = [] if len(f) < 8 or not f[7] else [p.rsplit("-", 1) for p in f[7].split(",")]
    partners = {x for x, n in parts if x not in (c, "chrUn") and int(n) >= a.min_haps}
    end = s < a.end_mb * 1e6 or e > size[c] - a.end_mb * 1e6
    W.append(dict(c=c, shared=bool(partners), end=end, fa=cover(ma.get(c, []), s, e) / (e - s),
                  fi=cover(mi.get(c, []), s, e) / (e - s)))
edges = [float(x) for x in a.classes.split(",")]
cls = lambda x: bisect.bisect_right(edges, x) - 1
n_sh = sum(w["shared"] for w in W)
obs = sum(w["shared"] and w["end"] for w in W) / n_sh
print(f"windows {len(W)}, at an end {sum(w['end'] for w in W)} ({sum(w['end'] for w in W) / len(W):.4f}); "
      f"shared (>= 1 other chromosome, >= {a.min_haps} haplotypes) {n_sh}; shared at an end "
      f"{sum(w['shared'] and w['end'] for w in W)} = {obs:.4f} of shared windows")
rng = random.Random(a.seed)
for key, lab in (("fi", "inter-chromosomal SD"), ("fa", "any SD")):
    by = defaultdict(list)
    for w in W:
        by[cls(w[key])].append(w)
    print(f"-- null within {lab} density classes {edges}")
    for k in sorted(by):
        v = by[k]
        sh = sum(w["shared"] for w in v)
        print(f"   class {edges[k]}-{edges[k + 1]}: windows {len(v)}, shared {sh}, at an end {sum(w['end'] for w in v)}, "
              f"shared at an end {sum(w['shared'] and w['end'] for w in v)}")
    null = []
    for _ in range(a.perm):
        hit = 0
        for v in by.values():
            sh = sum(w["shared"] for w in v)
            hit += sum(w["end"] for w in rng.sample(v, sh))
        null.append(hit / n_sh)
    null.sort()
    p = (1 + sum(x >= obs for x in null)) / (1 + a.perm)
    print(f"   expected share of shared windows at an end {sum(null) / len(null):.4f} "
          f"(95% {null[int(0.025 * a.perm)]:.4f}-{null[int(0.975 * a.perm)]:.4f}); observed {obs:.4f}; "
          f"enrichment {obs / (sum(null) / len(null)):.2f}x; one-sided p {p:.2g} ({a.perm} permutations)")
