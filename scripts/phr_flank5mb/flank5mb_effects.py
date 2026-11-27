#!/usr/bin/env python3
"""R1-M8: split the two ways the 5 Mb design differs from the paper's 500 kb design, per flank PHR length.
  partner set   : paper / 500 kb rerun with all 18,827 flanks as partners  vs  500 kb with only the acro p + X/Y flanks
  partner length: 500 kb with only the acro p + X/Y flanks  vs  5 Mb run replayed within 500 kb (flank5mb_compare.py derive TSV)
Lengths are region_end - region_start of the caller (0 = no PHR)."""
import argparse, csv, statistics as st
from collections import Counter

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--all", required=True, help="500 kb calls, all partners (rerun)")
ap.add_argument("--paper", required=True)
ap.add_argument("--acroxy", required=True, help="500 kb calls, acro p + X/Y partners only")
ap.add_argument("--derive", required=True, help="flank5mb_compare.py derive TSV of the 5 Mb run vs --acroxy")
ap.add_argument("--acroxy-windows", help="per-window rows of --acroxy; with --w5mb-windows and --names: where the windows that")
ap.add_argument("--w5mb-windows", help="pass only at 5 Mb lie (distance from the telomere) and which chromosomes they reach")
ap.add_argument("--names", help="5 Mb name -> 500 kb name")
a = ap.parse_args()
rd = lambda f: {r["seq"]: r for r in csv.DictReader(open(f), delimiter="\t")}
L = lambda r: 0 if r["region_start"] == "." else int(r["region_end"]) - int(r["region_start"])
arm = lambda s: s.rsplit("_", 2)[-2] + ("p" if s.endswith("parm") else "q")
al, pa, ax = rd(a.all), rd(a.paper), rd(a.acroxy)
ks = sorted(ax)


def report(lab, x, y):
    d = [(s, y[s] - x[s]) for s in ks if x[s] != y[s]]
    print(f"{lab}: same {len(ks) - len(d)} of {len(ks)}, shorter {sum(v < 0 for _, v in d)}, longer {sum(v > 0 for _, v in d)}, "
          f"median change {st.median(v for _, v in d) if d else 0} bp; differing by arm {Counter(arm(s) for s, _ in d).most_common(9)}")
    print(f"   at 500 kb: {sum(v >= 500000 for v in x.values())} -> {sum(v >= 500000 for v in y.values())}; "
          f"no PHR: {sum(v == 0 for v in x.values())} -> {sum(v == 0 for v in y.values())}")


lall, lpap, lax = ({s: L(t[s]) for s in ks} for t in (al, pa, ax))
r = {x["seq500"]: int(x["len_replay500"]) for x in csv.DictReader(open(a.derive), delimiter="\t")}
assert set(r) == set(ks), "derive TSV and acro+XY calls hold different flanks"
report("paper -> rerun, all partners (caller version)", lpap, lall)
report("rerun all partners -> acro+XY partners only (partner set)", lall, lax)
report("acro+XY partners 500 kb -> 5 Mb replayed within 500 kb (partner length)", lax, r)
report("paper -> 5 Mb replayed within 500 kb (both)", lpap, r)
if a.acroxy_windows:
    m5 = dict(l.rstrip("\n").split("\t") for l in open(a.names))
    def win(f, keep):
        d = {}
        for x in csv.DictReader(open(f), delimiter="\t"):
            if x["seq"] in keep:
                d[(x["seq"], int(x["dist_from_telomere"]))] = x
        return d
    w5, w500 = win(a.w5mb_windows, set(m5)), win(a.acroxy_windows, set(m5.values()))
    gained = [(s5, k, w5[(s5, k)]["chrs_involved"], w500[(m5[s5], k)]["chrs_involved"]) for (s5, k) in w5
              if k + 5000 <= 500000 and (m5[s5], k) in w500 and w5[(s5, k)]["passes_threshold"] == "True" and w500[(m5[s5], k)]["passes_threshold"] == "False"]
    print(f"windows within 500 kb that pass at 5 Mb but not at 500 kb (same partners): {len(gained)}")
    for lo in range(0, 500000, 100000):
        print(f"   distance {lo // 1000}-{lo // 1000 + 100} kb: {sum(lo <= k < lo + 100000 for _, k, _, _ in gained)}")
    by = Counter()
    for s5, k, c5, c500 in gained:
        by[(arm(s5), c5)] += 1
    print("   top (arm, chromosomes reached at 5 Mb):", by.most_common(8))
    near = Counter(arm(s5) for s5, k, _, _ in gained if k < 100000)
    print("   within 100 kb of the telomere, by arm:", near.most_common(6))
    for s5, k, c5, c500 in [g for g in gained if g[1] < 100000][:4]:
        print(f"   e.g. {s5} at {k}: 5 Mb {c5} | 500 kb {c500}")
