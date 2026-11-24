#!/usr/bin/env python3
"""Known-answer test for geneconv_families.py: two arms of 12 copies each (3% fixed between-arm
differences, 0.3% private differences), one planted 600 bp tract copied from a 10q copy into a 4q copy,
and the same alignment without the plant. Pass = a cross-arm GI fragment involving the planted copy
and overlapping the tract, and no cross-arm GI fragment without the plant.
Run: python3 scripts/conversion/test_geneconv_families.py GENECONV_BINARY WORKDIR"""
import os, random, subprocess, sys

geneconv, wd = sys.argv[1], sys.argv[2]
here = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(7)
L, N = 3000, 12
base = [rng.choice("ACGT") for _ in range(L)]
other = {"A": "C", "C": "G", "G": "T", "T": "A"}
arm10 = list(base)
for j in rng.sample(range(L), int(0.03 * L)):
    arm10[j] = other[arm10[j]]


def copies(anc, arm):
    out = {}
    for i in range(N):
        s = list(anc)
        for j in rng.sample(range(L), int(0.003 * L)):
            s[j] = other[s[j]]
        out[f"{arm}_HG{i:03d}_1"] = s
    return out


seqs = {**copies(base, "4q"), **copies(arm10, "10q")}
plant = {k: list(v) for k, v in seqs.items()}
plant["4q_HG005_1"][1000:1600] = seqs["10q_HG007_1"][1000:1600]
ok = True
for name, s in (("plant", plant), ("noplant", seqs)):
    fa = os.path.join(wd, f"TEST_{name}_arrays.fa")
    os.makedirs(wd, exist_ok=True)
    with open(fa, "w") as fh:
        for k, v in s.items():
            fh.write(f">{k}\n{''.join(v)}\n")
    out = os.path.join(wd, f"out_{name}")
    subprocess.run([sys.executable, os.path.join(here, "geneconv_families.py"), fa, "--geneconv", geneconv,
                    "--out", out], check=True)
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, "fragments.tsv"))][1:]
    cross = [r for r in rows if r[1] == "cross_arm"]
    if name == "plant":
        hit = [r for r in cross if "4q_HG005_1" in r[3] and int(r[6]) < 1600 and int(r[7]) > 1000]
        good = len(hit) >= 1
        print(f"{'PASS' if good else 'FAIL'} planted tract found: {[(r[3], r[6], r[7], r[4]) for r in hit]}")
    else:
        good = cross == []
        print(f"{'PASS' if good else 'FAIL'} no cross-arm fragment without plant: {len(cross)}")
    ok &= good
sys.exit(0 if ok else 1)
