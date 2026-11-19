#!/usr/bin/env python3
"""One table of candidate assembly errors, v1.1.

usage: g1_evidence.py WASHU OUT_TSV
Candidates = union of (a) ONT clusters (check/wg_scan/clusters/<s>.v11.clusters.tsv, not satellite_block, >= 20 sites),
(b) merqury error k-mer clusters (merqury/v11/<s>/<s>_hap{1,2}_only.bed.gz, k-mers <= 10 kb apart, >= 20 k-mers),
(c) G2 residual gaps >= 1 kb (g2/homolog/residual.mm.tsv, v11; query = PAN027 or PAN028). Merged when <= 10 kb apart.
Per candidate: ONT sites and fraction, merqury k-mers, cenSat bp and classes, released problematic labels, distance to end.
"""
import collections, glob, gzip, sys

D, out = sys.argv[1:3]
SM = ("PAN010", "PAN011", "PAN027", "PAN028")
def aname(c):  # PAN027#1#chr1 -> PAN027.chr1.maternal, PAN011#2#chr9 -> PAN011.chr9.haplotype2
    s, h, ch = c.split("#")
    return f"{s}.{ch}." + ({"1": "maternal", "2": "paternal"}[h] if s == "PAN027" else f"haplotype{h}")
def load_bed(f):
    d = collections.defaultdict(list)
    for l in open(f):
        t = l.rstrip("\n").split("\t")
        if len(t) >= 3 and t[1].isdigit(): d[t[0]].append((int(t[1]), int(t[2]), t[3] if len(t) > 3 else ""))
    return d
def ov(iv, a, b):
    x = [(min(b, e) - max(a, s), n) for s, e, n in iv if min(b, e) > max(a, s)]
    return sum(v for v, _ in x), ",".join(sorted({n.split("(")[0] for _, n in x if n}))[:60]
cand = collections.defaultdict(list)  # contig -> [start, end, source]
ont = collections.defaultdict(list)
for s in SM:
    for l in list(open(f"{D}/check/wg_scan/clusters/{s}.v11.clusters.tsv"))[1:]:
        t = l.rstrip("\n").split("\t")
        if t[13] == "satellite_block" or int(t[9]) < 20: continue
        a, b = int(t[1]), int(t[2]); cand[t[0]].append((a, b, "ont")); ont[t[0]].append((a, b, t))
mq = collections.defaultdict(list)
for s in SM:
    for h in (1, 2):
        pos = collections.defaultdict(list)
        for l in gzip.open(f"{D}/merqury/v11/{s}/{s}_hap{h}_only.bed.gz", "rt"):
            t = l.split("\t"); pos[t[0]].append(int(t[1]))
        for c, p in pos.items():
            p.sort(); mq[c] = p; g = [p[0]]
            for x in p[1:] + [None]:
                if x is not None and x - g[-1] <= 10000: g.append(x); continue
                if len(g) >= 20: cand[c].append((g[0], g[-1] + 21, "merqury"))
                if x is not None: g = [x]
for l in list(open(f"{D}/g2/homolog/residual.mm.tsv"))[1:]:
    t = l.rstrip("\n").split("\t")
    if t[0] != "v11": continue
    q = {"pat": "PAN027#2", "mat": "PAN027#1"}.get(t[1])
    if q is None:
        q = "PAN028#2" if t[2] in "chr1 chr3 chr4 chr5 chr7 chr9 chr12 chr14 chr15 chr18 chr21 chr22 chrX".split() else "PAN028#1"
    cand[f"{q}#{t[2]}"].append((int(t[3]), int(t[4]), f"g2_{t[1]}"))
fai = {}
for s in SM:
    for h in (1, 2):
        for l in open(f"{D}/norm_v11/{s}#{h}.fa.gz.fai"): t = l.split("\t"); fai[t[0]] = int(t[1])
cs = {s: load_bed(f"{D}/annot_v11/{s}.polished.diploid.cenSat.bed") for s in SM}
pb = {s: load_bed(f"{D}/annot_v11/problematic.{s}.bed") for s in SM}
rows = []
for c, iv in cand.items():
    iv.sort(); m = []
    for a, b, src in iv:
        if m and a - m[-1][1] <= 10000: m[-1][1] = max(m[-1][1], b); m[-1][2].add(src)
        else: m.append([a, b, {src}])
    s = c.split("#")[0]; an = aname(c); L = fai.get(c, 0)
    for a, b, src in m:
        o = [x for x in ont[c] if x[1] > a and x[0] < b]
        sites = sum(int(x[2][9]) for x in o); fr = max((float(x[2][10]) for x in o if x[2][10] != "nan"), default=float("nan"))
        cls = ",".join(sorted({x[2][13] for x in o}))
        import bisect
        k = bisect.bisect_left(mq[c], b) - bisect.bisect_left(mq[c], a) if c in mq else 0
        cb, cn = ov(cs[s][an], a, b); pbp, pn = ov(pb[s][an], a, b)
        rows.append((c, a, b, b - a, min(a, L - b), "+".join(sorted(src)), sites, f"{fr:.3f}", cls, k, cb, cn, pbp, pn))
with open(out, "w") as o:
    o.write("contig\tstart\tend\tlen\tdist_end\tsources\tont_sites\tont_frac_max\tont_class\tmerqury_kmers\tcensat_bp\tcensat\tproblematic_bp\tproblematic\n")
    for r in sorted(rows): o.write("\t".join(map(str, r)) + "\n")
n = collections.Counter(r[5] for r in rows)
print(f"candidates {len(rows)}, bp {sum(r[3] for r in rows)}", file=sys.stderr)
for k, v in n.most_common(): print(f"  {k}: {v}", file=sys.stderr)
both = [r for r in rows if "ont" in r[5] and "merqury" in r[5]]
print(f"ONT and merqury: {len(both)}, of which in cenSat {sum(1 for r in both if r[10] > 0)}, in problematic {sum(1 for r in both if r[12] > 0)}", file=sys.stderr)
