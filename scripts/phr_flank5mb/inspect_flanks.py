#!/usr/bin/env python3
"""For given 5 Mb flanks: alignment records in their all-vs-one PAF (as target), how many start in the first 5 kb
from the telomere, and the first per-window rows of the PHR caller. Used to explain flanks whose PHR call differs
between the 500 kb and the 5 Mb runs (R1-M8)."""
import argparse, csv, gzip, os, sys
from collections import Counter

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--paf-dir", required=True)
ap.add_argument("--per-window", required=True)
ap.add_argument("--fai", required=True)
ap.add_argument("--prefix", default="", help="PAF file name prefix before the sanitized target name")
ap.add_argument("seqs", nargs="+")
a = ap.parse_args()
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai)}
want = set(a.seqs)
win = {s: [] for s in want}
for r in csv.DictReader(open(a.per_window), delimiter="\t"):
    if r["seq"] in want:
        win[r["seq"]].append(r)
for s in a.seqs:
    f = os.path.join(a.paf_dir, a.prefix + s.translate(str.maketrans(":@#-", "____")) + ".paf.gz")
    if not os.path.exists(f):
        sys.exit(f"FATAL: no PAF {f}")
    n = other = near = 0
    arms = Counter()
    p = s.endswith("_parm")
    for l in gzip.open(f, "rt"):
        x = l.split("\t")
        n += 1
        if x[0] == x[5]:
            continue
        other += 1
        ts, te = int(x[7]), int(x[8])
        d = ts if p else L[s] - te  # distance of the record from the telomere end of the target
        if d < 5000:
            near += 1
            arms[x[0].rsplit("_", 2)[-2]] += 1
    w = sorted(win[s], key=lambda r: int(r["dist_from_telomere"]))[:5]
    print(f"{s} (flank {L[s]} bp): PAF records {n}, from other flanks {other}, within 5 kb of the telomere {near} "
          f"by chromosome {dict(arms.most_common(6))}")
    for r in w:
        print(f"   window at {r['dist_from_telomere']}: chromosomes {r['n_chrs']} ({r['chrs_involved']}) pass {r['passes_threshold']}")
