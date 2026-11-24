#!/usr/bin/env python3
"""Do chromosome ends share the positional profile of the acrocentric and pseudoautosomal ends? (R3-M2)

Profile of an end: for the 100 kb CHM13 windows of the genome-wide survey (Fig. 1A input), ordered from the
telomere inward over --span-mb, the number of chromosomes other than its own (chrUn excluded) reached by
>= --min-haps haplotypes. Reference profile: the mean over the acrocentric short arms (13p, 14p, 15p, 21p,
22p) and the X pseudoautosomal ends (Xp, Xq); the CHM13v2.0 Y has PAR1 masked and is left out.
For every other end: Pearson r between its profile and the reference. Null: the end's profile circularly
shifted by 1..n-1 windows (keeps its own autocorrelation, breaks the alignment to the telomere);
p = (1 + shifts with r >= observed) / n. An end with no sharing in the span has no profile (r undefined).
"""
import argparse, gzip, statistics as st
from collections import defaultdict

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--windows", default="data/hprc25272.CHM13.w100kb-xm5-id098-l5k.tsv.gz")
ap.add_argument("--min-haps", type=int, default=47)
ap.add_argument("--span-mb", type=float, default=5.0)
ap.add_argument("--alpha", type=float, default=0.05)
a = ap.parse_args()
rows = [l.rstrip("\n").split("\t") for l in gzip.open(a.windows, "rt") if not l.startswith("chrom")]
size, val = defaultdict(int), {}
for f in rows:
    c = f[0].split("#")[-1]
    s, e = int(f[1]), int(f[2])
    size[c] = max(size[c], e)
    parts = [] if len(f) < 8 or not f[7] else [p.rsplit("-", 1) for p in f[7].split(",")]
    val[(c, s)] = len({x for x, n in parts if x not in (c, "chrUn") and int(n) >= a.min_haps})
n = int(a.span_mb * 10)
prof = {}
for c in size:
    if c == "chrM":
        continue
    p = [val.get((c, i * 100000), 0) for i in range(n)]
    last = (size[c] - 1) // 100000 * 100000
    q = [val.get((c, last - i * 100000), 0) for i in range(n)]
    prof[c.replace("chr", "") + "p"], prof[c.replace("chr", "") + "q"] = p, q
ref_ends = ["13p", "14p", "15p", "21p", "22p", "Xp", "Xq"]
ref = [st.mean(prof[e][i] for e in ref_ends) for i in range(n)]


def pearson(x, y):
    mx, my = st.mean(x), st.mean(y)
    sx = sum((v - mx) ** 2 for v in x) ** 0.5
    sy = sum((v - my) ** 2 for v in y) ** 0.5
    return None if sx == 0 or sy == 0 else sum((u - mx) * (v - my) for u, v in zip(x, y)) / (sx * sy)


print(f"span {a.span_mb} Mb = {n} windows; partner chromosome counted at >= {a.min_haps} haplotypes; "
      f"reference ends {ref_ends}; reference profile first 10 windows " + ",".join(f"{v:.1f}" for v in ref[:10]))
print("end\tr\tp_circular\tterminal_window_partners\tmax_partners_in_span\tverdict")
res = []
for e in sorted(prof, key=lambda x: (x[-1] == "q", int(x[:-1]) if x[:-1].isdigit() else 99, x)):
    if e in ref_ends or e in ("Yp", "Yq"):
        continue
    x = prof[e]
    r = pearson(x, ref)
    if r is None:
        res.append((e, None, None))
        print(f"{e}\tNA\tNA\t{x[0]}\t{max(x)}\tno sharing in span")
        continue
    null = [pearson(x[k:] + x[:k], ref) for k in range(1, n)]
    p = (1 + sum(v is not None and v >= r for v in null)) / n
    res.append((e, r, p))
    print(f"{e}\t{r:.3f}\t{p:.3f}\t{x[0]}\t{max(x)}\t{'same profile' if p < a.alpha else 'not significant'}")
ok = [x for x in res if x[2] is not None and x[2] < a.alpha]
print(f"SUMMARY: {len(res)} ends tested; {len(ok)} share the reference profile (circular-shift p < {a.alpha}); "
      f"{sum(x[1] is None for x in res)} without sharing in the span; "
      f"{sum(x[2] is not None and x[2] >= a.alpha for x in res)} with sharing but not significant")
