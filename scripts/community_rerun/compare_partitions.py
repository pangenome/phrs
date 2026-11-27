#!/usr/bin/env python3
"""Compare two community assignments (TSV with a header; the first column is the item, the community is the column named
Community or community, else the second column). On the items in both: number of communities, whether the partitions
are identical up to labels, the adjusted Rand index, and the items whose group of co-members changed.
usage: compare_partitions.py A.tsv B.tsv [label]"""
import collections, csv, math, sys

def load(f):
    r = list(csv.reader(open(f), delimiter="\t")); h = r[0]
    c = next((i for i, x in enumerate(h) if x.lower() == "community"), 1)
    return {x[0]: x[c] for x in r[1:] if x}

a, b = load(sys.argv[1]), load(sys.argv[2]); lab = sys.argv[3] if len(sys.argv) > 3 else ""
common = sorted(set(a) & set(b))
ga, gb = collections.defaultdict(set), collections.defaultdict(set)
for x in common:
    ga[a[x]].add(x); gb[b[x]].add(x)
same = {frozenset(v) for v in ga.values()} == {frozenset(v) for v in gb.values()}
n = len(common); c2 = lambda k: k * (k - 1) / 2
ct = collections.Counter((a[x], b[x]) for x in common)
sij = sum(c2(v) for v in ct.values()); si = sum(c2(len(v)) for v in ga.values()); sj = sum(c2(len(v)) for v in gb.values())
exp = si * sj / c2(n) if n > 1 else 0; mx = (si + sj) / 2
ari = (sij - exp) / (mx - exp) if mx != exp else 1.0
moved = [x for x in common if ga[a[x]] != gb[b[x]]]
print(f"{lab}: items A {len(a)}, B {len(b)}, common {n}; communities A {len(ga)}, B {len(gb)}; identical partition {same}; "
      f"ARI {ari:.4f}; items whose co-members changed {len(moved)}")
if moved and len(ga) <= 60 and n <= 100:
    for x in moved:
        print(f"  {x}: A with {sorted(ga[a[x]] - {x})} | B with {sorted(gb[b[x]] - {x})}")
