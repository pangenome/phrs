#!/usr/bin/env python3
"""Profile one copy of a duplicon alignment (duplicon_monophyly.py --msa-dir) against the consensus of
each arm, per bin of reference columns (R1-M1, test 3a). A converted copy matches its own arm in some
bins and another arm in others; a divergent outlier differs from every arm everywhere. Consensus of an
arm = majority base over that arm's copies (the profiled copy excluded) with >= 90% non-N columns.
"""
import argparse, collections

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--msa", required=True)
ap.add_argument("--copy", required=True)
ap.add_argument("--bin", type=int, default=250)
a = ap.parse_args()
seqs, name = {}, None
for line in open(a.msa):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:]
        seqs[name] = []
    else:
        seqs[name].append(line)
seqs = {k: "".join(v) for k, v in seqs.items()}
if a.copy not in seqs:
    raise SystemExit(f"FATAL: {a.copy} not in {a.msa}")
u, L = seqs[a.copy], len(seqs[a.copy])
arms = sorted({k.split("_")[0] for k in seqs})
cons, n = {}, {}
for arm in arms:
    vs = [v for k, v in seqs.items() if k.split("_")[0] == arm and k != a.copy and sum(c != "N" for c in v) >= 0.9 * L]
    n[arm] = len(vs)
    cons[arm] = "".join((collections.Counter(c for c in (s[j] for s in vs) if c in "ACGT").most_common(1)
                         or [("N", 0)])[0][0] for j in range(L)) if vs else "N" * L
print(f"copy {a.copy}; arms (copies in consensus): " + ", ".join(f"{x} ({n[x]})" for x in arms))
print("cols\t" + "\t".join(f"diff_vs_{x}" for x in arms) + "\tcompared")
for s in range(0, L, a.bin):
    e = min(L, s + a.bin)
    ok = [j for j in range(s, e) if u[j] in "ACGT" and all(cons[x][j] in "ACGT" for x in arms)]
    print(f"{s + 1}-{e}\t" + "\t".join(str(sum(u[j] != cons[x][j] for j in ok)) for x in arms) + f"\t{len(ok)}")
