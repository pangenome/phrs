#!/usr/bin/env python3
"""Copy-number check of a block with k-mers: read count vs copies in two assemblies (one with, one without the block).
usage: kmer_cn.py READS.wig ASM_WITH.wig ASM_WITHOUT.wig CTRL_READS.wig CTRL_ASM.wig
wig = meryl-lookup -wig-count output (k-mer count at each position). Per-copy read depth d = median read count of control
k-mers with 1 copy in the assembly. For block k-mers: est = reads / d; a k-mer favours the assembly whose copy count is
closer to est.
"""
import statistics, sys
def wig(f):  # position -> count; variableStep lines "pos count", fixedStep lines "count"; missing positions = 0
    v, pos, step = {}, None, 1
    for l in open(f):
        t = l.split()
        if not t or t[0].startswith(("track", "#")): continue
        if t[0] in ("variableStep", "fixedStep"):
            kv = dict(x.split("=") for x in t[1:]); pos = int(kv["start"]) if "start" in kv else None; step = int(kv.get("step", 1)); continue
        if len(t) == 2: v[int(t[0])] = float(t[1])
        else: v[pos] = float(t[0]); pos += step
    return v
R, A1, A0, CR, CA = (wig(f) for f in sys.argv[1:6])
P = sorted(set(R) | set(A1) | set(A0)); CP = sorted(set(CR) | set(CA))
print(f"positions: block {len(P)} (reads {len(R)}, with {len(A1)}, without {len(A0)}); control {len(CP)}")
r, a1, a0 = [R.get(i, 0) for i in P], [A1.get(i, 0) for i in P], [A0.get(i, 0) for i in P]
cr, ca = [CR.get(i, 0) for i in CP], [CA.get(i, 0) for i in CP]
d = statistics.median([x for x, y in zip(cr, ca) if y == 1 and x > 0])
print(f"per-copy read depth (control k-mers with 1 copy): {d}")
w = n0 = t = mis = 0; est = []
for x, y1, y0 in zip(r, a1, a0):
    if y1 == 0: continue
    t += 1; e = x / d; est.append(e)
    if x == 0: mis += 1
    if abs(e - y1) < abs(e - y0): w += 1
    elif abs(e - y0) < abs(e - y1): n0 += 1
print(f"block k-mers {t}; absent from reads {mis}; closer to the assembly WITH the block {w}, WITHOUT {n0}, tie {t - w - n0}")
print(f"median read copies {statistics.median(est):.2f}; median copies in assembly with block {statistics.median([y for y in a1 if y > 0]):.0f}, without {statistics.median([y for x, y in zip(a1, a0) if x > 0]):.0f}")
