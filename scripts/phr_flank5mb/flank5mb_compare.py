#!/usr/bin/env python3
"""Check the PHR caller rerun against the paper's calls, then measure PHR lengths on 5 Mb flanks (R1-M8).

exact  : NEW calls on the original 500 kb PAFs vs the paper's calls (all-vs-all.1Mb.p95.id95.len.tsv) and per-window
         rows; every field must match for the sequences NEW processed.
derive : 5 Mb run. From its per-window rows we replay the caller's rule on the windows of the first 500 kb from the
         telomere (pass = >= 2 chromosomes; the region runs from the telomere to the last passing window before the
         first run of 4 failing windows) and compare that length with the paper's call of the same flank (names map
         5 Mb -> 500 kb). Per window we also compare pass/fail with the paper's per-window rows. Then we report the
         5 Mb lengths of the PHRs that reached the 500 kb edge.
"""
import argparse, csv, statistics as st, sys
from collections import defaultdict

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("mode", choices=["exact", "derive"])
ap.add_argument("--orig-calls", required=True)
ap.add_argument("--orig-windows", required=True)
ap.add_argument("--new-calls", required=True)
ap.add_argument("--new-windows", required=True)
ap.add_argument("--names", help="5 Mb name -> 500 kb name (derive)")
ap.add_argument("--fai5", help="fai of the 5 Mb flanks (derive)")
ap.add_argument("--out", required=True, help="per-flank TSV")
a = ap.parse_args()
W, STOP = 5000, 4
rows = lambda f: list(csv.DictReader(open(f), delimiter="\t"))
log = lambda m: sys.stderr.write(m + "\n")


def windows(f, keep):
    """per-window rows of the sequences in keep, as {seq: {dist_from_telomere: (passes, chrs)}}"""
    d = defaultdict(dict)
    for r in csv.DictReader(open(f), delimiter="\t"):
        if r["seq"] in keep:
            d[r["seq"]][int(r["dist_from_telomere"])] = (r["passes_threshold"] == "True", r["chrs_involved"])
    return d


oc = {r["seq"]: r for r in rows(a.orig_calls)}
nc = {r["seq"]: r for r in rows(a.new_calls)}
length = lambda r: 0 if r["region_start"] == "." else int(r["region_end"]) - int(r["region_start"])
log(f"original calls {len(oc)} ({sum(length(r) > 0 for r in oc.values())} with a region); new calls {len(nc)}")
if a.mode == "exact":
    miss = sorted(set(nc) - set(oc))
    log(f"new sequences absent from the original calls: {len(miss)} {miss[:3]}")
    same = [s for s in nc if s in oc and all(nc[s][k] == oc[s][k] for k in oc[s])]
    diff = [s for s in nc if s in oc and s not in same]
    ow, nw = windows(a.orig_windows, set(nc)), windows(a.new_windows, set(nc))
    wsame = sum(ow[s] == nw[s] for s in nc)
    with open(a.out, "w") as fh:
        fh.write("seq\torig_start\torig_end\tnew_start\tnew_end\tsame_row\tsame_windows\n")
        for s in sorted(nc):
            o = oc.get(s, {"region_start": "NA", "region_end": "NA"})
            fh.write(f"{s}\t{o['region_start']}\t{o['region_end']}\t{nc[s]['region_start']}\t{nc[s]['region_end']}\t{s in same}\t{ow[s] == nw[s]}\n")
    print(f"EXACT: {len(same)} of {len(nc)} rows identical in every field (region, chromosomes, arms), {len(diff)} differ; "
          f"per-window rows identical for {wsame} of {len(nc)} sequences")
    both = [s for s in nc if s in oc]
    rdiff = [s for s in both if length(nc[s]) != length(oc[s]) or nc[s]["region_start"] != oc[s]["region_start"]]
    arm = lambda s: s.rsplit("_", 2)[-2] + ("p" if s.endswith("parm") else "q")
    print(f"REGIONS: identical for {len(both) - len(rdiff)} of {len(both)}; arms_involved identical for "
          f"{sum(nc[s]['arms_involved'] == oc[s]['arms_involved'] for s in both)}; PHRs paper {sum(length(oc[s]) > 0 for s in both)}, "
          f"new {sum(length(nc[s]) > 0 for s in both)}; PHR bp paper {sum(length(oc[s]) for s in both)}, new {sum(length(nc[s]) for s in both)}")
    flips = defaultdict(int)
    for s in nc:
        for k, (p, _) in ow[s].items():
            if k in nw[s] and p != nw[s][k][0]:
                flips["paper pass, new fail" if p else "paper fail, new pass"] += 1
                flips["inside the paper PHR" if k + W <= length(oc[s]) else "past the paper PHR end"] += 1
    print(f"WINDOWS: pass/fail differs in {dict(flips)} of {sum(len(nw[s]) for s in nc)} windows")
    for s in rdiff:
        print(f"  region differs: {s} ({arm(s)}) paper {length(oc[s])} new {length(nc[s])}")
    for s in diff[:3]:
        print("  e.g. field differs:", s, {k: (oc[s][k], nc[s][k]) for k in oc[s] if oc[s][k] != nc[s][k]})
    sys.exit(0)
m5 = dict(l.rstrip("\n").split("\t") for l in open(a.names))
L5 = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai5)}
nw = windows(a.new_windows, set(m5))
ow = windows(a.orig_windows, set(m5.values()))
log(f"names {len(m5)}; 5 Mb sequences with per-window rows {len(nw)}; their 500 kb flanks with per-window rows {len(ow)}")


def replay(win, limit):
    """caller rule on the windows that end within limit bp of the telomere, in telomere order; region length"""
    last, fail = -1, 0
    for i, k in enumerate(sorted(x for x in win if x + W <= limit)):
        if win[k][0]:
            last, fail = i, 0
        else:
            fail += 1
            if fail >= STOP:
                break
    return (last + 1) * W


out, cmp, wdiff, wtot, unaligned, selfok = [], defaultdict(int), defaultdict(int), 0, 0, 0
for s5, s500 in sorted(m5.items()):
    arm = s5.rsplit("_", 2)[-2] + ("p" if s5.endswith("parm") else "q")
    if s5 not in nc:
        sys.exit(f"FATAL: {s5} has no call in the 5 Mb run")
    d = nw.get(s5, {})
    if d and min(d) != 0:
        unaligned += 1  # q-arm flank whose length is not a multiple of 5 kb: windows are offset from the telomere
    l500, lorig, l5 = replay(d, 500000), length(oc[s500]), length(nc[s5])
    selfok += replay(ow.get(s500, {}), 500000) == lorig  # known answer: the replay on the paper's own windows
    cmp["same" if l500 == lorig else ("longer" if l500 > lorig else "shorter")] += 1
    for k, (p, ch) in ow.get(s500, {}).items():
        if k in d:
            wtot += 1
            if p != d[k][0]:
                wdiff["orig pass, 5 Mb fail" if p else "orig fail, 5 Mb pass"] += 1
    out.append(dict(seq5=s5, seq500=s500, arm=arm, len_paper=lorig, len_replay500=l500, len_5mb=l5,
                    flank_bp=L5[s5], at_edge=l5 > ((L5[s5] - W) // W + 1 - STOP) * W))
with open(a.out, "w") as fh:
    cols = list(out[0])
    fh.write("\t".join(cols) + "\n")
    for r in out:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
print(f"SELF-CHECK: the replay on the paper's per-window rows gives the paper's length for {selfok} of {len(out)} flanks")
print(f"DERIVE 500 kb: replayed length vs paper for {len(out)} flanks: same {cmp['same']}, longer {cmp['longer']}, "
      f"shorter {cmp['shorter']}; windows compared {wtot}, pass/fail differs {dict(wdiff)}; q flanks with offset windows {unaligned}")
cen = [r for r in out if r["len_paper"] >= 500000]
print(f"PHRs at the 500 kb edge in the paper: {len(cen)}")
by = defaultdict(list)
for r in cen:
    by[r["arm"]].append(r)
print("arm\tn\tmedian_5mb_Mb\tmin\tmax\tat_flank_edge (region ends in the last 4 windows: length is a floor)")
for arm in sorted(by, key=lambda x: (x[3:-1].isdigit() and int(x[3:-1]) or 99, x)):
    v = [r["len_5mb"] / 1e6 for r in by[arm]]
    print(f"{arm}\t{len(v)}\t{st.median(v):.2f}\t{min(v):.2f}\t{max(v):.2f}\t{sum(r['at_edge'] for r in by[arm])}")
for t in (1e6, 2e6, 3e6, 4e6):
    print(f"   censored in the paper with 5 Mb length >= {t / 1e6:.0f} Mb: {sum(r['len_5mb'] >= t for r in cen)} of {len(cen)}")
print(f"   censored in the paper still at the 5 Mb flank edge: {sum(r['at_edge'] for r in cen)} of {len(cen)}")
