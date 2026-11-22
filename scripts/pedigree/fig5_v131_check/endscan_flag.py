#!/usr/bin/env python3
"""Suspect ends from the whole-genome ONT end scan (end_scan.sbatch on ont_full), and their end-sweep task lines.

usage: endscan_flag.py SCAN_TSV TASKS_TSV[,TASKS_TSV...] OUT_TSV
Per end (contig, start|end): wrong = sum of wrong positions over its 10 kb bins; low = bins with MAPQ>=20 depth < 0.25 x
the median MAPQ>=20 depth of all non-zero bins. Flag = wrong >= 20, or wrong >= 5 and low >= 5. The known DAD chr9 h2 end
error has wrong 97, low 9 on v1.3.1. Writes the task lines (endsweep_tasks.py lists) of the flagged ends, renumbered.
"""
import collections, statistics, sys
scan, tasks, out = sys.argv[1], sys.argv[2].split(","), sys.argv[3]
rows = [l.rstrip("\n").split("\t") for l in open(scan)][1:]
md = statistics.median(float(r[4]) for r in rows if float(r[4]) > 0)
w, lo = collections.Counter(), collections.Counter()
for r in rows:
    k = (r[0], r[1]); w[k] += int(r[3]); lo[k] += float(r[4]) < 0.25 * md
flag = {k for k in w if w[k] >= 20 or (w[k] >= 5 and lo[k] >= 5)}
print(f"{scan}: {len(w)} ends, median q20 depth {md}, flagged {len(flag)}", file=sys.stderr)
for k in sorted(flag): print(f"  {k[0]} {k[1]}: wrong {w[k]}, low-depth bins {lo[k]}", file=sys.stderr)
names = {f"{c.split('#')[2]}h{c.split('#')[1]}{'p' if e == 'start' else 'q'}" for c, e in flag}
n = 0
with open(out, "w") as o:
    for tf in tasks:
        for l in open(tf):
            t = l.rstrip("\n").split("\t")
            if t[1].split("_w")[0] in names: n += 1; o.write("\t".join([str(n)] + t[1:]) + "\n")
print(f"task lines written: {n} (expected {len(flag)} x {len(tasks)})", file=sys.stderr)
if n != len(flag) * len(tasks): sys.exit("ERROR task lines missing for some flagged ends")
