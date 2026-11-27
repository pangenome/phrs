#!/usr/bin/env python3
"""CPU time of the all-vs-one wfmash runs, per target class (acrocentric p, X/Y, other arms).
500 kb run: /usr/bin/time -v in slurm-<job>_<task>.out, task i = line i of sequences.txt.
(The 5 Mb run has no CPU record: sacct TotalCPU is 0 for its tasks.)"""
import argparse, glob, os, re

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--dir500", required=True)
a = ap.parse_args()
ACRO = {"chr13_parm", "chr14_parm", "chr15_parm", "chr21_parm", "chr22_parm"}
cls = lambda s: "acro p" if s.rsplit("_", 2)[-2] + "_" + s.rsplit("_", 1)[-1] in ACRO else ("X/Y" if re.search(r"_chr[XY]_[pq]arm$", s) else "other")
seqs = [l.strip() for l in open(os.path.join(a.dir500, "sequences.txt"))]
tot, n, wall = {}, {}, {}
for f in glob.glob(os.path.join(a.dir500, "slurm-*_*.out")):
    i = int(f.rsplit("_", 1)[1].split(".")[0])
    t = open(f).read()
    u = re.search(r"User time \(seconds\): ([\d.]+)", t); s = re.search(r"System time \(seconds\): ([\d.]+)", t)
    w = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)", t)
    if not (u and s and w):
        raise SystemExit(f"FATAL: no time -v block in {f}")
    c = cls(seqs[i - 1])
    tot[c] = tot.get(c, 0) + float(u.group(1)) + float(s.group(1)); n[c] = n.get(c, 0) + 1
    p = [float(x) for x in w.group(1).split(":")]
    wall[c] = wall.get(c, 0) + (p[0] * 3600 + p[1] * 60 + p[2] if len(p) == 3 else p[0] * 60 + p[1])
print("500 kb all-vs-all (18,827 targets x 18,827 queries, 9.4 Gb):")
for c in sorted(tot):
    print(f"   {c}: {n[c]} targets, CPU {tot[c] / 3600:.0f} h ({tot[c] / n[c]:.0f} s per target), wall sum {wall[c] / 3600:.0f} h")
print(f"   total CPU {sum(tot.values()) / 3600:.0f} h")
