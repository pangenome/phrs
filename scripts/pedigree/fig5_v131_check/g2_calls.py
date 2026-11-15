#!/usr/bin/env python3
"""Exchange calls from G2c window tables.

usage: g2_calls.py WIN_DIR ANNOT_DIR WINDOW_KB OUT_TSV
A window is "other" when hom_id > 0 and other_id - hom_id >= delta, for delta in 0.0005, 0.001, 0.002 (other_id >= 0.99).
Runs <= 20 kb within 20 kb of a chromosome end are flagged "edge" (telomere repeats); longer terminal runs stay.
Consecutive "other" windows (gap <= 1 window) form a run. Per run: span, windows, max delta, targets,
query cenSat and released problematic overlap (v1.1 annotations).
"""
import collections, glob, os, sys

wd, ad, wk, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
DELTAS = (0.0005, 0.001, 0.002)
KIDM = set("chr1 chr3 chr4 chr5 chr7 chr9 chr12 chr14 chr15 chr18 chr21 chr22 chrX".split())
def qname(pair, c):
    if pair == "pat": return "PAN027", f"PAN027.{c}.paternal"
    if pair == "mat": return "PAN027", f"PAN027.{c}.maternal"
    return "PAN028", f"PAN028.{c}.haplotype{2 if c in KIDM else 1}"
def load(f):
    d = collections.defaultdict(list)
    for l in open(f):
        t = l.rstrip("\n").split("\t")
        if len(t) >= 3 and t[1].isdigit(): d[t[0]].append((int(t[1]), int(t[2]), t[3] if len(t) > 3 else ""))
    return d
ann = {s: (load(f"{ad}/{s}.polished.diploid.cenSat.bed"), load(f"{ad}/problematic.{s}.bed")) for s in ("PAN027", "PAN028")}
def ov(iv, a, b):
    x = [(min(b, e) - max(a, s), n) for s, e, n in iv if min(b, e) > max(a, s)]
    return sum(v for v, _ in x), ",".join(sorted({n.split("(")[0] for _, n in x if n}))[:50]
rows = []
for f in sorted(glob.glob(f"{wd}/*.w{wk}.tsv")):
    tag, pair, c = os.path.basename(f).split(".")[:3]
    w = [l.rstrip("\n").split("\t") for l in open(f)][1:]
    s, an = qname(pair, c)
    for dl in DELTAS:
        hit = [i for i, r in enumerate(w) if float(r[2]) > 0 and float(r[4]) >= 0.99 and float(r[4]) - float(r[2]) >= dl]
        runs = []
        for i in hit:
            if runs and i - runs[-1][-1] <= 2: runs[-1].append(i)
            else: runs.append([i])
        for r in runs:
            a, b = int(w[r[0]][0]), int(w[r[-1]][1])
            md = max(float(w[i][4]) - float(w[i][2]) for i in r)
            tg = collections.Counter(w[i][5].split("#")[1] + ":" + w[i][5].split("#")[-1] for i in r)
            cb, cn = ov(ann[s][0][an], a, b); pb, pn = ov(ann[s][1][an], a, b)
            edge = "edge" if (a < 20000 or int(w[-1][1]) - b < 20000) and b - a <= 20000 else "inner"  # long terminal runs kept
            rows.append((tag, pair, c, dl, a, b, b - a, len(r), f"{md:.5f}", ",".join(f"{k}x{v}" for k, v in tg.most_common(3)), cb, cn, pb, pn, edge))
with open(out, "w") as o:
    o.write("tag\tpair\tchrom\tdelta\tstart\tend\tspan\twindows\tmax_delta\ttargets\tcensat_bp\tcensat\tproblematic_bp\tproblematic\tposition\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
cnt = collections.Counter((r[0], r[1], r[3]) for r in rows)
for k in sorted(cnt): print(f"{k[0]} {k[1]} delta {k[2]}: runs {cnt[k]}, outside cenSat and problematic "
                            f"{sum(1 for r in rows if (r[0], r[1], r[3]) == k and r[10] == 0 and r[12] == 0)}, kept (not short end runs) "
                            f"{sum(1 for r in rows if (r[0], r[1], r[3]) == k and r[14] == 'inner')}")
