#!/usr/bin/env python3
"""For given flanks: the windows whose pass/fail differs between two per-window tables of the PHR caller
(paper vs rerun), with the chromosomes each table reached there; then per-arm PHR length summaries
(n, median, mean, max) from two call tables, for the arms of those flanks."""
import argparse, csv, statistics as st

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--a-windows", required=True)
ap.add_argument("--b-windows", required=True)
ap.add_argument("--a-calls", required=True)
ap.add_argument("--b-calls", required=True)
ap.add_argument("seqs", nargs="+")
x = ap.parse_args()
want = set(x.seqs)


def win(f):
    d = {}
    for r in csv.DictReader(open(f), delimiter="\t"):
        if r["seq"] in want:
            d[(r["seq"], int(r["dist_from_telomere"]))] = r
    return d


A, B = win(x.a_windows), win(x.b_windows)
for k in sorted(B):
    if k in A and A[k]["passes_threshold"] != B[k]["passes_threshold"]:
        print(f"{k[0]} window at {k[1]} bp from the telomere: A pass {A[k]['passes_threshold']} chromosomes {A[k]['chrs_involved']} "
              f"| B pass {B[k]['passes_threshold']} chromosomes {B[k]['chrs_involved']}")
    elif k not in A:
        print(f"{k[0]} window at {k[1]} bp: only in B (A stopped before), B pass {B[k]['passes_threshold']} chromosomes {B[k]['chrs_involved']}")
arm = lambda s: s.rsplit("_", 2)[-2] + ("_p" if s.endswith("parm") else "_q")
L = lambda r: 0 if r["region_start"] == "." else int(r["region_end"]) - int(r["region_start"])
for lab, f in (("A", x.a_calls), ("B", x.b_calls)):
    rows = [r for r in csv.DictReader(open(f), delimiter="\t") if "NA18982#1#JBKABS010000018.1" not in r["seq"]]
    for a in sorted({arm(s) for s in want}):
        v = [L(r) for r in rows if arm(r["seq"]) == a and L(r) > 0]
        print(f"{lab} {a}: PHRs {len(v)}, median {st.median(v)}, mean {st.mean(v):.1f}, max {max(v)}")
