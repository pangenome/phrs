#!/usr/bin/env python3
"""Sensitivity of the PHR calls (Methods: stopping rule, scan length, chromosomes per window), replayed exactly from the
caller's per-window rows (find-multichr-regions-incremental.py --per-window-output). A window's pass/fail does not
depend on these settings, and every variant scans a subset of the windows the default run scanned (an earlier stop, a
shorter scan, or a stricter chromosome count that fails at least where the default fails), so the replay gives the
caller's own result. Rule: windows in telomere order; a window passes with >= K chromosomes (n_chrs); the region runs
from the telomere to the last passing window before the first run of STOP failing windows; windows beyond MAX are not
scanned. Check: the default (K 2, STOP 4, MAX 500 kb) must reproduce the calls file. Prints PHR and arm counts per variant.
usage: phr_sensitivity.py PER_WINDOW.tsv CALLS.tsv [EXCLUDE.txt]"""
import csv, collections, sys

pw, calls = sys.argv[1], sys.argv[2]
excl = set(l.strip() for l in open(sys.argv[3])) if len(sys.argv) > 3 else set()
win = collections.defaultdict(list)
for r in csv.DictReader(open(pw), delimiter="\t"):
    if r["seq"] not in excl:
        win[r["seq"]].append((int(r["dist_from_telomere"]), int(r["n_chrs"])))
ref = {}
for r in csv.DictReader(open(calls), delimiter="\t"):
    if r["seq"] not in excl:
        ref[r["seq"]] = 0 if r["region_start"] == "." else int(r["region_end"]) - int(r["region_start"])
arm = lambda s: s.rsplit("_", 2)[-2] + s.rsplit("_", 1)[-1][0]


def call(ws, k=2, stop=4, mx=500000, w=5000):
    last, fails = -1, 0
    for d, n in sorted(ws):
        if d + w > mx:
            break
        if n >= k:
            last, fails = d, 0
        else:
            fails += 1
            if fails >= stop:
                break
    return 0 if last < 0 else last + w


def summary(lab, **kw):
    L = {s: call(win.get(s, []), **kw) for s in ref}
    arms = {arm(s) for s, v in L.items() if v > 0}
    print(f"{lab}: PHRs {sum(v > 0 for v in L.values())} of {len(L)} flanks, arms {len(arms)}")
    return L


d = summary("default (2 chromosomes, stop 4, 500 kb)")
bad = [s for s in ref if d[s] != ref[s]]
print(f"CHECK default vs calls file: {len(ref) - len(bad)} of {len(ref)} identical" + (f", e.g. {bad[:3]}" if bad else ""))
if bad:
    sys.exit("FATAL: the replay does not reproduce the calls")
for stop in (1, 2, 3):
    summary(f"stop after {stop} failed windows", stop=stop)
for mx in (25000, 50000, 100000, 200000, 300000, 400000):
    summary(f"scan limited to {mx // 1000} kb", mx=mx)
for k in (3, 6, 15):
    summary(f"{k} chromosomes per window", k=k)
